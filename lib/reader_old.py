import os
import sys
import struct
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))

from xy_rphi import *
from file_utils import get_prefix
from tools.configuration import *

# Data layout in fields files
fields = [ "Ex", "Ey", "Ez", "Bx", "By", "Bz" ]

slices = {
    "X": [ "0071" ],
    "Y": [ "0071" ],
    "Z": [ "0201" ],
}

# Pressures remap
pressures = {
    "Prr": "Pxx",
    "Paa": "Pyy"
}

data_shape = {
    'X': (Ny + 3, Nz + 3),
    'Y': (Nx + 3, Nz + 3),
    'Z': (Nx + 3, Ny + 3),
}

COS, SIN, R_MAP = init_COS_SIN_RMAP((0, data_shape['Z'][0], 0, data_shape['Z'][1]))


def get_filepath(t: int, path: str, prefix=None):
    p = get_prefix(t) if prefix == None else prefix
    t_str = str(t).zfill(4)
    return f"{p}/{path}_{t_str}"

def get_fields_path(plane):
    return f"Fields/Diag2D/FieldPlane{plane}_{slices[plane][-1]}"

def get_fields_file(t, plane="Z", prefix=None):
    return get_filepath(t, get_fields_path(plane), prefix)

def get_particles_path(sort, diag_name, plane):
    return f"Particles/{sort}/Diag2D/{diag_name}Plane{plane}_{slices[plane][-1]}"

def get_particles_file(sort, diag_name, plane, t, prefix=None):
    return get_filepath(t, get_particles_path(sort, diag_name, plane), prefix)

def parse_file(path, offset=0):
    float_byte_size = 4

    with open(path, "rb") as file:
        buffer = file.read(2 * float_byte_size)
        nx_ny = struct.unpack("ff", buffer[: 2 * float_byte_size])
        nx = int(nx_ny[0])
        ny = int(nx_ny[1])

        raw = np.fromfile(
            file,
            dtype=np.float32,
            count=(nx * ny),
            offset=(nx * ny) * float_byte_size * offset,
        )
    return raw.reshape((nx, ny)).transpose()

def get_parsed_field(path, name, plane, comp, t, prefix=None):
    file = get_filepath(t, path, prefix)
    if comp == 'z':
        return parse_file(file, fields.index(name + comp))
    elif (plane == "X" and comp == "x"):
        # we return A_phi, thus we should invert the second half in y
        data = parse_file(file, fields.index(name + comp))
        data[:, (data.shape[1] // 2 + 1):] *= -1
        return data
    elif (plane == "X" and comp == "y") or \
         (plane == "Y" and comp in "xy"):
        data = parse_file(file, fields.index(name + comp))
        return data
    elif plane == "Z":
        fx = parse_file(file, fields.index(name + "x"))
        fy = parse_file(file, fields.index(name + "y"))
        return vx_vy_to_vr_va(fx, fy, COS, SIN)

def get_parsed_scalar(field, t):
    return parse_file(f"{get_prefix(t)}/{field.path_to_file}_{str(t).zfill(4)}")
