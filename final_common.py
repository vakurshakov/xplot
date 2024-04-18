#!/usr/bin/python3

from lib_common import *

import np_2000.parameters as p2000
import np_5000.parameters as p5000

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

res_dir = f"../Final"
mkdir(res_dir)

# We take currents from np_5000 because they are less noisy. It is also
# valid because we want to compare the result before the m=3 instability
def get_current_data(t, sort, params=p5000):
    t_str = str(int(t)).zfill(4)
    prefix = get_prefix(t, params.restart_timesteps, params.prefixes)
    filename = f"{prefix}/Particles/{sort}/Diag2D/CurrentPlaneAvgZ{t_str}"
    je_x = parse_file(filename, 0)
    je_y = parse_file(filename, 1)
    return vx_vy_to_vr_va(je_x, je_y, COS, SIN)[1]
