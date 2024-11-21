#!/usr/bin/env python3

from lib_common import *

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


res_dir = f"{params_path}/Final"
mkdir(res_dir)
