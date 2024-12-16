import os
import sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import cumulative_trapezoid

import struct
import numpy as np

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc = comm.Get_size()


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

def agg(to_agg, data):
    return data + (to_agg if np.any(to_agg) else np.zeros_like(data))

def dump(name, units, t, data):
    prefix = f"{params_path}/Dump"
    mkdir(prefix)
    with open(f"{prefix}/{name}_{t}.txt", "w") as f:
        f.write(f"Grid {name} [{units}]\n")
        for i, d in enumerate(data):
            f.write(f"{i} {d}\n")

def sliding_average(linear, w=5):
    return np.convolve(linear, np.ones(w) / w, mode='valid')

# Berendeev's data reading utilities
def get_prefix(t, restarts=restart_timesteps, prefixes=prefixes):
    i = 0
    for restart in restarts:
        if (t > restart):
            i += 1
    return prefixes[i]

def get_parsed_file(t, path, prefix=None):
    p = get_prefix(t) if prefix == None else prefix
    t_str = str(int(t)).zfill(4)
    return f"{p}/{path}_{t_str}"

def get_fields_path(plane):
    return f"Fields/Diag2D/FieldPlane{plane}_{slices[plane][-1]}"

def get_fields_file(t, plane="Z", prefix=None):
    return get_parsed_file(t, get_fields_path(plane), prefix)

def get_particles_path(sort, diag_name, plane):
    return f"Particles/{sort}/Diag2D/{diag_name}Plane{plane}_{slices[plane][-1]}"

def get_particles_file(sort, diag_name, plane, t, prefix=None):
    return get_parsed_file(t, get_particles_path(sort, diag_name, plane), prefix)

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

# Timestep consistency utils
def is_correct_timestep(t):
    plane = "Z"
    fields_file = get_fields_file(t, plane)
    fields_file_bytesize = 4 * (2 + data_shape[plane][0] * data_shape[plane][1] * len(fields))
    return os.path.isfile(fields_file) and os.path.getsize(fields_file) == fields_file_bytesize

def check_consistency(tmin, tmax):
    for t in range(tmin, tmax):
        if not is_correct_timestep(t):
            print(f"Data is inconsistent. Valid data range is: ({tmin}, {t}) [dts] or ({tmin * dts / tau}, {t * dts / tau}) [tau].")
            return
    print(f"Data range ({tmin}, {tmax}) [dts] or ({tmin * dts / tau}, {t * dts / tau}) [tau] is consistent.")
    return

def timestep_should_be_processed(t, filename, skip_processed=True):
    msg = f"{filename} {t} [dts] {t * dts / tau:.3f} [tau]"
    if not is_correct_timestep(t):
        print(msg, "Data is incorrect, skipping")
        return False
    if skip_processed and os.path.exists(filename):
        print(msg, "Timestep was already processed, skipping.")
        return False
    print("Processing", msg)
    return True

def find_correct_timestep(t, t_range):
    for t_c in range(t_range[0], t + 1)[::-1]:
        if not is_correct_timestep(t_c):
            print(f"Warning! Timestep {t_c} is incorrect, first previous correct step will be used.")
            continue
        print(f"{t_c:4d} [dts]", f"{t_c * dts / tau:6.3f}", "[tau]")
        return t_c
    return t_c

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
def get_parsed_field(path, name, plane, comp, t, prefix=None):
    file = get_parsed_file(t, path, prefix)
    if comp == 'z':
        return parse_file(file, fields.index(name + comp))
    elif (plane == "X" and comp == "x"):
        # we return A_phi, thus we should invert the second half in y
        data = parse_file(file, fields.index(name + comp))
        data[:, (data_shape[plane][0] // 2 + 1):] *= -1
        return data
    elif (plane == "X" and comp == "y") or \
         (plane == "Y" and comp in "xy"):
        data = parse_file(file, fields.index(name + comp))
        data[:, :(data_shape[plane][0] // 2)] *= -1
        return data
    elif plane == "Z":
        fx = parse_file(file, fields.index(name + "x"))
        fy = parse_file(file, fields.index(name + "y"))
        return vx_vy_to_vr_va(fx, fy, COS, SIN)

def get_parsed_scalar(field, t):
    return parse_file(f"{get_prefix(t)}/{field.path_to_file}_{str(t).zfill(4)}")

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

def electric_field(plane, subplot=None, title=None):
    vmap = (-2e-2, +2e-2)
    field = Field(get_fields_path(plane), subplot, None, signed_cmap, vmap)
    if title != None:
        generate_info(field, plane, title)
    return field

def magnetic_field(plane, subplot=None, title=None):
    DB = B0
    cmap = unsigned_cmap if plane == 'Z' else signed_cmap
    vmap = (0, B0) if plane == 'Z' else (B0 - DB, B0 + DB)
    field = Field(get_fields_path(plane), subplot, None, cmap, vmap)
    if title != None:
        generate_info(field, plane, title)
    return field

def particles_field(sort, diag_name, plane, subplot=None, title=None, v=None, cmap=signed_cmap):
    field = Field(get_particles_path(sort, diag_name, plane), subplot, None, cmap, v)
    if title != None:
        generate_info(field, plane, title)
    return field

# magnetic-tubes
def select_magnetic_line(bz, xl):
    xs = data_shape["Y"][0]
    zs = data_shape["Y"][1]
    xc = xs // 2
    zc = zs // 2

    # 2 * np.pi can be removed since comparison is relative
    b_f0 = np.sum(bz[zc,xc:xl] * np.arange(0, xl - xc))

    xmap = np.zeros(zs, dtype=np.int32)
    for z in np.arange(0, zs):
        b_fz = 0
        for x in np.arange(xc, xs):
            b_fz += bz[z, x] * (x - xc)
            xmap[z] = x

            if (b_fz >= b_f0):
                break
    return xmap

def curvature_coefficients(bz_normalized, xl, xs1, xs2):
    cos_map = np.ones_like(bz_normalized[:,xl])
    for z, x in enumerate(0.5 * (xs1 + xs2)):
        cos_map[z] = bz_normalized[z,int(x)]
    return cos_map

def average_over_magnetic_tube(data, zs, xs1, xs2, xc):
    f_l = np.zeros(data_shape["Y"][1])
    for z in zs:
        rs = (np.arange(xs1[z], xs2[z] + 1) - xc) * dx
        area = cumulative_trapezoid(2 * np.pi * rs, rs, initial=0)[-1]
        f_l[z] += cumulative_trapezoid(data[z, xs1[z]:(xs2[z] + 1)] * (2 * np.pi * rs), rs, initial=0)[-1] / area
        f_l[z] += cumulative_trapezoid(data[z, (2 * xc - (xs2[z] + 1)):(2 * xc - xs1[z])] * (2 * np.pi * rs), rs, initial=0)[-1] / area
    return f_l / 2

def align(fr, fz, br, bz, b):
    dot = (fr * br + fz * bz)
    return np.divide(dot, b, where=(b > 1e-3), out=np.zeros_like(b))
