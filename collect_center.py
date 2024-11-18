#!/usr/bin/python3

from lib_common import *
from plot_fields import *
from plot_particles import *

res_dir = f"./{params_path}/Other"
mkdir(res_dir)

tmin = 0
tmax = int(time / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1))

b_arr = []

ni_arr = []
pri_arr = []
pai_arr = []

ne_arr = []
pre_arr = []
pae_arr = []

pd_arr = []

b_min_arr = []
beta_max_arr = []

B0 = 0.191
_, _, R, _ = init_XY_RA((0, data_shape['Z'][0], 0, data_shape['Z'][1]))

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

    _b = center(b)
    _pri = phi_averaged(pri, R_MAP)
    _pre = phi_averaged(pre, R_MAP)

    ni = center(ni)
    ne = center(ne)
    pai = phi_averaged(pai, R_MAP)
    pae = phi_averaged(pae, R_MAP)

    pr = _pri + _pre
    pa = pai + pae

    rs = np.arange(0, data_shape["Z"][0] // 2) * dx
    pd = cumulative_trapezoid((pr - pa) / (rs + 0.1), dx=dx, initial=0)
    pd = -(pd - pd[-1])

    b_min = None

    if t < 0.1 * tau / dts:
        b_min = np.where((R < 3 / dx), b, B0)
    elif t < 3 * tau / dts:
        b_min = np.where((R < 12 / dx), b, B0)
    else:
        b_min = np.where((R < 20 / dx), b, B0)

    b_min = np.min(b_min)

    beta_max = np.where((R < 20 / dx), 2 * (pri + pre) / (B0 * B0), 0)
    beta_max = np.max(beta_max)
    return _b, b_min, beta_max, ni, ne, _pri[0], pai[0], _pre[0], pae[0], pd[0]

def find_correct_timestep(t):
    for t_corr in range(t_range[0], t + 1, 1)[::-1]:
        if is_correct_timestep1(t_corr):
            print(f"{t_corr:5d} [dts]", f"{t_corr * dts / tau:6.3f}", "[tau]")
            break
        if t_corr == t_range[0]:
            print(f"Warning! Timestep is incorrect, last correct step will be used.")
    return t_corr

arr_name = [
    (b_arr, "b"),
    (b_min_arr, "b_min"),
    (beta_max_arr, "beta_max"),
    (ni_arr, "ni"),
    (ne_arr, "ne"),
    (pri_arr, "pri"),
    (pai_arr, "pai"),
    (pre_arr, "pre"),
    (pae_arr, "pae"),
    (pd_arr, "pd"),
]

for t in t_range:
    data = parse_data(find_correct_timestep(t))

    for d, (arr, _) in zip(data, arr_name):
        arr.append(d)

for arr, name in arr_name:
    gathered_list = comm.gather(arr, root=0)

    if (rank == 0):
        arr = aggregate_array(gathered_list)
        np.save(f"{res_dir}/{name}_t", arr)
