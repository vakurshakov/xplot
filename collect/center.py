#!/usr/bin/env python3

from lib_common import *

res_dir = f"{params_path}/Other"
mkdir(res_dir)

tmin = 0
tmax = int(time / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1))

def parse_data(t):
    b = get_parsed_field(magnetic_field("Z"), "B", "Z", "z", t)

    ni = get_parsed_scalar(particles_field("Ions", "Dens", "Z"), t)
    pri = get_parsed_scalar(particles_field("Ions", pressures["Prr"], "Z"), t)
    pai = get_parsed_scalar(particles_field("Ions", pressures["Paa"], "Z"), t)

    ne = get_parsed_scalar(particles_field("Electrons", "Dens", "Z"), t)
    pre = get_parsed_scalar(particles_field("Electrons", pressures["Prr"], "Z"), t)
    pae = get_parsed_scalar(particles_field("Electrons", pressures["Paa"], "Z"), t)

    def center(d, w=5):
        c0 = data_shape["Z"][0] // 2
        c1 = data_shape["Z"][1] // 2
        return np.mean(d[c1-w:c1+w, c0-w:c0+w])

    b = center(b)
    pri = phi_averaged(pri, R_MAP)
    pre = phi_averaged(pre, R_MAP)

    ni = center(ni)
    ne = center(ne)
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

for t in t_range:
    data = parse_data(find_correct_timestep(t))

    for d, (_, arr) in zip(data, named_arrays):
        arr.append(d)

for name, arr in named_arrays:
    gathered_list = comm.gather(arr, root=0)

    if (rank == 0):
        arr = aggregate_array(gathered_list)
        np.save(f"{res_dir}/{name}_t", arr)
