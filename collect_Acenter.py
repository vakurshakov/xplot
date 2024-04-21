#!/usr/bin/python3

from final_common import *

import scipy.integrate as I

# beta calculation
res_dir = f"./{params_path}/Data1D"
mkdir(res_dir)

tmin = 0
tmax = int(9 * tau / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1), rank, proc)

b_arr = []
ni_arr = []
pri_arr = []
pai_arr = []
pre_arr = []
pae_arr = []
pd_arr = []

def parse_data(t):
    b = parse_file(get_fields_file(t), fields.index("Bz"))
    ni = parse_file(get_particles_file("Ions", "DensPlaneAvgZ", t))

    pri = parse_file(get_particles_file("Ions", f"{pressures['Prr']}PlaneAvgZ", t))
    pai = parse_file(get_particles_file("Ions", f"{pressures['Paa']}PlaneAvgZ", t))
    pre = parse_file(get_particles_file("Electrons", f"{pressures['Prr']}PlaneAvgZ", t))
    pae = parse_file(get_particles_file("Electrons", f"{pressures['Paa']}PlaneAvgZ", t))

    b = phi_averaged(b, R_MAP) / B0
    ni = phi_averaged(ni, R_MAP)
    pri = phi_averaged(pri, R_MAP) / ((B0 * B0) / 2)
    pai = phi_averaged(pai, R_MAP) / ((B0 * B0) / 2)
    pre = phi_averaged(pre, R_MAP) / ((B0 * B0) / 2)
    pae = phi_averaged(pae, R_MAP) / ((B0 * B0) / 2)

    pr = pri + pre
    pa = pai + pae

    rs = np.arange(0, data_shape[0] // 2) * dx
    pd = I.cumtrapz((pr - pa) / (rs + 0.1), dx=dx, initial=0)
    pd = -(pd - pd[-1])

    return b, ni, pri, pai, pre, pae, pd

def find_correct_timestep(t):
    for t_corr in range(t_range[0], t + 1, 1)[::-1]:
        if is_correct_timestep(t_corr):
            print(f"{t_corr:5d} [dts]", f"{t_corr * dts / tau:6.3f}", "[tau]")
            break
        if t_corr == t_range[0]:
            print(f"Warning! Timestep is incorrect, last correct step will be used.")
    return t_corr

for t in t_range:
    b, ni, pri, pai, pre, pae, pd = parse_data(find_correct_timestep(t))

    b_arr.append(b[0])
    ni_arr.append(ni[0])
    pri_arr.append(pri[0])
    pai_arr.append(pai[0])
    pre_arr.append(pre[0])
    pae_arr.append(pae[0])
    pd_arr.append(pd[0])

arr_name = [
    (b_arr, "b"),
    (ni_arr, "ni"),
    (pri_arr, "pri"),
    (pai_arr, "pai"),
    (pre_arr, "pre"),
    (pae_arr, "pae"),
    (pd_arr, "pd"),
]

for arr, name in arr_name:
    gathered_list = comm.gather(arr, root=0)

    if (rank == 0):
        arr = aggregate_array(gathered_list)
        np.save(f"{res_dir}/{name}_t", arr)
