import os
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import cumulative_trapezoid

import struct
import numpy as np

# from mpi4py import MPI

# comm = MPI.COMM_WORLD
rank = 0 # comm.Get_rank()
proc = 1 # comm.Get_size()


from lib_plot import *
from lib_xy_rphi import *

# from nx_140_np_200.parameters import *
# from nx_140_np_1000.parameters import *
from nx_140_np_1000_glinskiy.parameters import *

# Data layout in fields files
fields = [ "Ex", "Ey", "Ez", "Bx", "By", "Bz" ]

# Pressures remap
pressures = {
    "Prr": "Pxx",
    "Paa": "Pyy"
}

sorts = [ "Electrons", "Ions" ]

boundaries = {
    'X': (-0.5 * Ny * dy, +0.5 * Ny * dy, 0, Nz * dz),
    'Y': (-0.5 * Nx * dx, +0.5 * Nx * dx, 0, Nz * dz),
    'Z': (-0.5 * Nx * dx, +0.5 * Nx * dx, -0.5 * Ny * dy, +0.5 * Ny * dy),
}

data_shape = {
    'X': (Ny + 3, Nz + 3),
    'Y': (Nx + 3, Nz + 3),
    'Z': (Nx + 3, Ny + 3),
}


COS, SIN, R_MAP = init_COS_SIN_RMAP((0, data_shape['Z'][0], 0, data_shape['Z'][1]))

def mkdir(dirname):
    if not os.path.exists(dirname) and rank == 0:
        os.mkdir(dirname)

mkdir(f"./{params_path}/Video")

# Berendeev's data reading utilities
def get_prefix(timestep, restarts=restart_timesteps, prefixes=prefixes):
    i = 0
    for restart in restarts:
        if (timestep > restart):
            i += 1
    return prefixes[i]

def get_particles_file(sort, diag_name, timestep, prefix=None):
    p = get_prefix(timestep) if prefix == None else prefix
    t_str = str(int(timestep)).zfill(4)
    return f"{p}/Particles/{sort}/Diag2D/{diag_name}{t_str}"

def get_fields_file(timestep, prefix=None):
    p = get_prefix(timestep) if prefix == None else prefix
    t_str = str(int(timestep)).zfill(4)
    return f"{p}/Fields/Diag2D/FieldPlaneZ_{t_str}"

def parse_file(path, offset=0):
    file = open(path, "rb")

    float_bytesize = 4
    buffer = file.read(2 * float_bytesize)
    nx_ny = struct.unpack("ff", buffer[: 2 * float_bytesize])
    nx = int(nx_ny[0])
    ny = int(nx_ny[1])

    raw = np.fromfile(
        file,
        dtype=np.float32,
        offset=(nx * ny) * float_bytesize * offset,
        count=(nx * ny),
    )
    return raw.reshape((nx, ny)).transpose()

def is_correct_timestep(t):
    plane = "Z"
    fields_file = f"{get_prefix(t)}/Fields/Diag2D/FieldPlane{plane}_{slices[plane][-1]}_{str(t).zfill(4)}"
    fields_file_bytesize = 4 * (2 + data_shape[plane][0] * data_shape[plane][1] * len(fields))
    return os.path.isfile(fields_file) and os.path.getsize(fields_file) == fields_file_bytesize

def check_consistency(tmin: int, tmax: int):
    for t in range(tmin, tmax):
        if not is_correct_timestep(t):
            print(f"Data is incosistent. Valid data range is: ({tmin}, {t}) [dts] or ({tmin * dts / tau}, {t * dts / tau}) [tau].")
            return

    print(f"Data range ({tmin}, {tmax}) [dts] or ({tmin * dts / tau}, {t * dts / tau}) [tau] is consistent.")
    return

# Figure utilities
def subplot(fig, gs, x, y):
    return fig.add_subplot(gs[x + y * gs.ncols])

bbox = dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.25')

def annotate_x(axis, annotation, x=0.5, y=1, size=big, ha="center"):
    axis.annotate(
        annotation,
        xy=(x, y),
        xytext=(0, 1),
        xycoords="axes fraction",
        textcoords="offset points",
        ha=ha,
        va="baseline",
        size=size,
        bbox=bbox
    )

def annotate_y(axis, annotation):
    axis.annotate(
        annotation,
        xy=(0, 0.5),
        xytext=(-axis.yaxis.labelpad - 1, 0),
        xycoords=axis.yaxis.label,
        textcoords="offset points",
        ha="right",
        va="center",
        rotation=90,
        size=big,
    )

def sliding_average(linear, w=5):
    return np.convolve(linear, np.ones(w) / w, mode='valid')


# MPI utilities
def create_t_range(tmin, tmax, offset):
    return np.arange(tmin + rank * offset, tmax + 1, proc * offset)

def reduce_array(a):
    length = len(a)
    chunk = length // proc
    remain = length % proc
    prev = 0

    if rank < remain:
        prev = (chunk + 1) * rank
        chunk += 1
    else:
        prev = chunk * rank + remain

    return a[prev: prev + chunk]

def aggregate_array(a):
    result = np.array([])
    for s in a:
        result = np.append(result, s)

    return result


# Fourier transform utilities. We seek for decomposition by exp(-iwt + ikx)
def fourier_transform(data):
    f_data = np.fft.ifftshift(data)
    f_data = np.fft.fft2(f_data)
    f_data = np.fft.fftshift(f_data)

    # np.fft uses linear frequenses, so we multiply the arrays by (2 * np.pi)
    shape = data.shape
    w = np.fft.fftfreq(shape[0], d=dts) * (2 * np.pi)
    k = np.fft.fftfreq(shape[1], d=(2 * np.pi /  shape[1])) * (2 * np.pi)

    w = np.fft.fftshift(w)
    k = np.fft.fftshift(k)

    return f_data, w, k

def inverse_fourier_transform(f_data):
    data = np.fft.ifftshift(f_data)
    data = np.fft.ifft2(data)
    data = np.fft.fftshift(data)

    r_data = np.real(data)
    i_data = np.imag(data)

    return r_data, i_data


# 3D-specific part
def generate_info(diag, plane, title):
    axes_args = {
        'X': [ "$(y, z, x=0)$", 'y', 'z' ],
        'Y': [ "$(x, z, y=0)$", 'x', 'z' ],
        'Z': [ "$(x, y, z=100)$", 'x', 'y' ]
    }

    diag.boundaries = boundaries[plane]
    bx = boundaries[plane][0] + buff * dx
    ex = boundaries[plane][1] - buff * dx
    by = boundaries[plane][2] + (buff * dy if plane == 'Z' else 0)
    ey = boundaries[plane][3] - (buff * dy if plane == 'Z' else 0)

    diag.set_axes_args(
        title=title + axes_args[plane][0],
        xlim=(bx, ex),
        ylim=(by, ey),
        xlabel=f"${axes_args[plane][1]},""~c/\\omega_{pe}$",
        ylabel=f"${axes_args[plane][2]},""~c/\\omega_{pe}$",
        xticks=np.linspace(bx, ex, 5),
        yticks=np.linspace(by, ey, 5),
    )

def timestep_should_be_processed(t, filename):
    if not is_correct_timestep1(t):
        print(f"Timestep {t} [dts], {t * dts / tau} is incorrect, it would be skipped.")
        return False
    # if os.path.exists(filename):
    #     print(f"Timestep {t} [dts], {t * dts / tau} was already processed, it would be skipped.")
    #     return False
    return True

def get_parsed_field(field, name, plane, comp, t):
    filename = f"{get_prefix(t)}/{field.path_to_file}_{str(t).zfill(4)}"
    if comp == 'z':
        return parse_file(filename, fields.index(name + comp))
    elif (plane == "X" and comp == "x"):
        # we return A_phi, thus we should invert the second half in y
        data = parse_file(filename, fields.index(name + comp))
        data[:, (data_shape[plane][0] // 2 + 1):] *= -1
        return data
    elif (plane == "X" and comp == "y") or \
         (plane == "Y" and comp in "xy"):
        data = parse_file(filename, fields.index(name + comp))
        data[:, :(data_shape[plane][0] // 2)] *= -1
        return data
    elif plane == "Z":
        e_x = parse_file(filename, fields.index(f"{name}x"))
        e_y = parse_file(filename, fields.index(f"{name}y"))
        return vx_vy_to_vr_va(e_x, e_y, COS, SIN)

def get_parsed_scalar(field, t):
    return parse_file(f"{get_prefix(t)}/{field.path_to_file}_{str(t).zfill(4)}")
