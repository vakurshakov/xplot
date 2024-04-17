import os
import sys
sys.path.append("../")

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import struct
import numpy as np

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc = comm.Get_size()


from lib_plot import *
from lib_xy_rphi import *

# from np_2000.parameters import * 
from np_5000.parameters import * 

# Data layout in fields files
fields = [ "Ex", "Ey", "Ez", "Bx", "By", "Bz" ]

# Pressures remap
pressures = {
    "Prr": "Pxx",
    "Paa": "Pyy"
}

sorts = [ "Ions", "Electrons" ]

# remap into c / w_pi units
dx /= np.sqrt(mi_me)
dy /= np.sqrt(mi_me)
dz /= np.sqrt(mi_me)
buff = 10

boundaries = (-0.5 * Nx * dx, +0.5 * Nx * dx, -0.5 * Ny * dy, +0.5 * Ny * dy)

data_shape = (Ny + 3, Nx + 3)
COS, SIN, R_MAP = init_COS_SIN_RMAP((0, data_shape[0], 0, data_shape[1]))

def mkdir(dirname):
    if not os.path.exists(dirname) and rank == 0:
        os.mkdir(dirname)

mkdir(f"./{params_path}/Video")

# Berendeev's data reading utilities
def get_prefix(timestep):
    i = 0
    for restart_timestep in restart_timesteps:
        if (timestep > restart_timestep):
            i += 1
    return prefixes[i]

def get_particles_file(sort, diag_name, timestep):
    t_str = str(int(timestep)).zfill(4)
    return f"{get_prefix(timestep)}/Particles/{sort}/Diag2D/{diag_name}{t_str}"

def get_fields_file(timestep):
    t_str = str(int(timestep)).zfill(4)
    return f"{get_prefix(timestep)}/Fields/Diag2D/FieldAvgPlaneZ_{t_str}"

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
    fields_file = get_fields_file(t)
    fields_file_bytesize = 4 * (2 + data_shape[0] * data_shape[1] * len(fields))
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

def annotate_x(axis, annotation, x=0.5, y=1):
    axis.annotate(
        annotation,
        xy=(x, y),
        xytext=(0, 1),
        xycoords="axes fraction",
        textcoords="offset points",
        ha="center",
        va="baseline",
        size=big,
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
def reduce_array(a, rank, proc):
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
        