import os
import sys
import numpy as np

from scipy.integrate import cumulative_trapezoid

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from plot import *
from xy_rphi import *
from mpi_utils import *
from file_utils import *

from tools.configuration import *

# Boundaries to all possible graphs (bx, ex, by, ey) [c / wpe]
boundaries = {
    'X': (-0.5 * Ny * dy, +0.5 * Ny * dy, 0, Nz * dz),
    'Y': (-0.5 * Nx * dx, +0.5 * Nx * dx, 0, Nz * dz),
    'Z': (-0.5 * Nx * dx, +0.5 * Nx * dx, -0.5 * Ny * dy, +0.5 * Ny * dy),
}

mkdir(f"{params_path}/Video")

def agg(to_agg, data):
    return data + (to_agg if np.any(to_agg) else np.zeros_like(data))

def sliding_average(linear, w=5):
    return np.convolve(linear, np.ones(w) / w, mode='valid')

# `Fields` generation
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