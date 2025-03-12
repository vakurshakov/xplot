import os
import sys
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../../"))

from xy_rphi import *
from file_utils import get_prefix
from tools.configuration import *

# Common shape of the data stored in files [cells]
data_shape = {
    'X': (Ny, Nz),
    'Y': (Nx, Nz),
    'Z': (Nx, Ny),
}

COS, SIN, R_MAP = init_COS_SIN_RMAP((0, data_shape['Z'][0], 0, data_shape['Z'][1]))

def get_diag_path(diag: dict | None, t: int, prefix: str = None):
    if diag == None:
        return None

    def get_suffix_2D(diag):
        plane = get(diag, "region.plane")
        pos = get(diag, "region.position")
        fill = -1

        if plane == "X": fill = Nx
        elif plane == "Y": fill = Ny
        elif plane == "Z": fill = Nz

        s = ""
        s += f"_Plane{plane}"
        s += f"_{str(pos).zfill(len(str(fill)))}cwpe"
        return s

    suffix = ""

    if diag["diagnostic"] in fields_diagnostics:
        suffix += get(diag, "field")
    elif diag["diagnostic"] in particles_diagnostics:
        suffix += f"{get(diag, "particles")}/{get(diag, "moment")}"

    if get(diag, "region.type") == "2D":
        suffix += get_suffix_2D(diag)

    p = get_prefix(t, restarts, prefixes) if prefix == None else prefix
    t_str = str(t).zfill(len(str(Nt)))
    return f"{p}/{suffix}/{t_str}"

def parse_file(diag: dict, t: int, prefix: str = None):
    ds = [Nx, Ny, Nz]

    if get(diag, "region.size") != None:
        ds = get(diag, "region.size")
    elif get(diag, "region.type") == "2D":
        ds = list(data_shape[get(diag, "region.plane")])

    count = np.prod(ds)
    if diag["diagnostic"] in fields_diagnostics and "comp" not in diag:
        ds.append(3)
        count *= ds[-1]

    with open(get_diag_path(diag, t) + ".bin", "rb") as file:
        raw = np.fromfile(
            file,
            dtype=np.float32,
            count=count,
            offset=0,
        )
    return raw.reshape(ds) #.transpose()

def get_parsed_field_cyl(diag, comp, t, prefix: str = None):
    plane = get(diag, "region.plane")

    if comp == 'z':
        return parse_file(diag, t)
    elif (plane == "X" and comp == "x"):
        # we return A_phi, thus we should invert the second half in y
        data = parse_file(diag, t)
        data[:, (data.shape[1] // 2 + 1):] *= -1
        return data
    elif (plane == "X" and comp == "y") or \
         (plane == "Y" and comp in "xy"):
        data = parse_file(diag, t)
        return data
    elif plane == "Z":
        fx = parse_file(diag, t)
        fy = parse_file(diag, t)
        return vx_vy_to_vr_va(fx, fy, COS, SIN)

def get_parsed_field(diag, t, prefix: str = None):
    return get_parsed_field_cyl(diag, t, prefix)

def get_parsed_scalar(diag, t, prefix: str = None):
    return parse_file(diag, t, prefix)