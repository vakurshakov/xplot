import os
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))

from xplot.lib.plot import *
from xplot.lib.reader import *
from xplot.lib.mpi_utils import *

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

# `Fields` generation utils
def generate_info(field: Field, title: str = None):
    if title == None:
        return

    axes_args = {
        'X': [ "$(y, z, x=0)$", 'y', 'z' ],
        'Y': [ "$(x, z, y=0)$", 'x', 'z' ],
        'Z': [ "$(x, y, z=100)$", 'x', 'y' ]
    }

    plane = get(field.path_to_file, "region.plane")

    field.boundaries = boundaries[plane]
    bx = boundaries[plane][0] + buff * dx
    ex = boundaries[plane][1] - buff * dx
    by = boundaries[plane][2] + (buff * dy if plane == 'Z' else 0)
    ey = boundaries[plane][3] - (buff * dy if plane == 'Z' else 0)

    field.set_axes_args(
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
    field = Field(find_diag(f"FieldView.E.{plane}"), subplot, None, signed_cmap, vmap)
    generate_info(field, title)
    return field

def magnetic_field(plane, cmap, vmap, subplot=None, title=None):
    field = Field(find_diag(f"FieldView.B.{plane}"), subplot, None, cmap, vmap)
    generate_info(field, title)
    return field

def particles_field(name, cmap=None, v=None, subplot=None, title=None):
    field = Field(find_diag(name), subplot, None, signed_cmap if cmap == None else cmap, v)
    generate_info(field, title)
    return field