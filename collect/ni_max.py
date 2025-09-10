#!/usr/bin/env python3

from collect import *
from scipy.ndimage import uniform_filter

xc = data_shape["Z"][0] // 2
yc = data_shape["Z"][1] // 2
xw = 2
yw = 2

def parse(t):
    if t == 0:
        return 0, 0, 0
        
    ni = parse_file(get_particles_file("Ions", "Density", "Z", t))
    ni = uniform_filter(ni, size=25, mode="constant")

    i = np.argmax(ni)
    ny, nx = np.unravel_index(i, ni.shape)

    v = np.mean(ni[ny-yw:ny+yw,nx-xw:nx+xw])
    x = (nx-xc)*dx
    y = (ny-yc)*dx

    return v, x, y

named_arrays = [
    ["ni_max_value", []],
    ["ni_max_x", []],
    ["ni_max_y", []],
]

def output(name):
    return f"{res_dir}/{name}_t"

process_collection(named_arrays, parse, output)
