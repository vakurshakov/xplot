#!/usr/bin/env python3

from collect import *

r0    = int(1.0 / dx)
rmax  = int(4.0 / dx)
rstep = 5

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
    ("b", []),
    ("er", []),
    ("ea", []),
    ("jri", []),
    ("jai", []),
    ("jre", []),
    ("jae", []),
]

def output(name):
    return f"{res_dir}/{name}_phit_r={r:.2f}"

for r, map in rmap:
    process_collection(lambda t: parse(t, map), named_arrays, output, ((tmax - tmin), len(map[0])))
