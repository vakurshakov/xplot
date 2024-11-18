#!/usr/bin/env python3

from lib_common import *

import np_2000.parameters as p2000
import np_5000.parameters as p5000
import mi1.parameters     as pmi1

bx = -60
ex = +60
by = bx
ey = ex

arg_1d = {
    "xlabel": "$r,~c/\\omega_{pe}$",
    "xlim": (0, ex),
    "xticks": np.linspace(0, ex, 5),
}

arg_2d = {
    "xlabel": "$x,~c/\\omega_{pe}$",
    "ylabel": "$y,~c/\\omega_{pe}$",
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

arg_arrow = {
    "head_width": 0.01,
    "head_length": 5,
    "length_includes_head": True,
    "color": "black",
    "linewidth": 2,
    "zorder": 100
}

arg_anno = {
    "arrowprops": dict(
        facecolor="black",
        width=1,
        headwidth=6
    ),
    "fontsize": ssmol,
    "horizontalalignment": "left",
    "bbox": bbox,
}


res_dir = f"../Final"
mkdir(res_dir)

def get_current_data(t, sort, prefix):
    t_str = str(int(t)).zfill(4)
    filename = f"{prefix}/Particles/{sort}/Diag2D/CurrentPlaneAvgZ{t_str}"
    je_x = parse_file(filename, 0)
    je_y = parse_file(filename, 1)
    return vx_vy_to_vr_va(je_x, je_y, COS, SIN)

def get_magnetic_field_data(t, prefix):
    t_str = str(int(t)).zfill(4)
    filename = f"{prefix}/Fields/Diag2D/FieldAvgPlaneZ_{t_str}"
    return parse_file(filename, fields.index("Bz"))

def get_electric_fields_data(t, prefix):
    t_str = str(int(t)).zfill(4)
    filename = f"{prefix}/Fields/Diag2D/FieldAvgPlaneZ_{t_str}"
    e_x = parse_file(filename, fields.index("Ex"))
    e_y = parse_file(filename, fields.index("Ey"))
    return vx_vy_to_vr_va(e_x, e_y, COS, SIN)