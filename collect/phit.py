#!/usr/bin/env python3

from collect import *

r0    = int(8 / dx)   # dx units
rmax  = int(16 / dx)  # dx units
rstep = 5             # dimensionless

rmap = []
for r, map in enumerate(R_MAP):
    if (r0 <= r and r < rmax) and (r % rstep == 0):
        rmap.append((r*dx, map))

def parse(t, map):
    b = get_parsed_field(get_fields_path("Z"), "B", "Z", "z", t)
    er, ea = get_parsed_field(get_fields_path("Z"), "E", "Z", "", t)
    jri, jai = get_parsed_field(get_particles_path("Ions", "Current", "Z"), "E", "Z", "", t)
    jre, jae = get_parsed_field(get_particles_path("Electrons", "Current", "Z"), "E", "Z", "", t)
    return b[map], er[map], ea[map], jri[map], jai[map], jre[map], jae[map]

named_arrays = [
    ["b", []],
    ["er", []],
    ["ea", []],
    ["jri", []],
    ["jai", []],
    ["jre", []],
    ["jae", []],
]

def output(name, r):
    return f"{res_dir}/{name}_phit_r={r:.2f}"

for r, map in rmap:
    process_collection(
        named_arrays,
        lambda t: parse(t, map),
        lambda name: output(name, r),
        shape=((tmax - tmin), len(map[0])))

    # Clearing the arrays, after collection
    for _, arr in named_arrays:
        arr = []