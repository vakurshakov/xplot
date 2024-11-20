#!/usr/bin/env python3

from lib_collect import *

res_dir = f"{params_path}/Other"
mkdir(res_dir)

def parse_data(t):
    b = get_parsed_field("B", "Z", "z", t)

    ni = parse_file(get_particles_file("Ions", "Dens", "Z", t))
    pri = parse_file(get_particles_file("Ions", pressures["Prr"], "Z", t))
    pai = parse_file(get_particles_file("Ions", pressures["Paa"], "Z", t))

    ne = parse_file(get_particles_file("Electrons", "Dens", "Z", t))
    pre = parse_file(get_particles_file("Electrons", pressures["Prr"], "Z", t))
    pae = parse_file(get_particles_file("Electrons", pressures["Paa"], "Z", t))

    b = center_avg(b)
    pri = phi_averaged(pri, R_MAP)
    pre = phi_averaged(pre, R_MAP)

    ni = center_avg(ni)
    ne = center_avg(ne)
    pai = phi_averaged(pai, R_MAP)
    pae = phi_averaged(pae, R_MAP)

    pr = pri + pre
    pa = pai + pae

    rs = np.arange(0, data_shape["Z"][0] // 2) * dx
    pd = cumulative_trapezoid((pr - pa) / (rs + 0.1), dx=dx, initial=0)
    pd = -(pd - pd[-1])
    return b, ni, ne, pri[0], pai[0], pre[0], pae[0], pd[0]

named_arrays = [
    ("b", []),
    ("ni", []),
    ("ne", []),
    ("pri", []),
    ("pai", []),
    ("pre", []),
    ("pae", []),
    ("pd", []),
]

def output_file(name):
    return f"{res_dir}/{name}_t"

process(parse_data, named_arrays, output_file)
