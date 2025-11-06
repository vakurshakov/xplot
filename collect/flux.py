#!/usr/bin/env python3

from collect import *

X, Y, R, _ = init_XY_RA((0, data_shape['Z'][0], 0, data_shape['Z'][1]))

x_lim = 50 # data_shape["Z"][0]
y_lim = 50 # data_shape["Z"][1]
r_lim = 50

def parse(t):
    b = get_parsed_field(magnetic_field("Z"), "B", "Z", "z", t)
    return np.sum(b, where=np.logical_and(np.abs(X) < x_lim, np.abs(Y) < y_lim)) * (2 * x_lim * dx) * (2 * y_lim * dy),

named_arrays = [
    ["flux", []],
]

def output(name):
    return f"{res_dir}/{name}_t"

process_collection(named_arrays, parse, output)
