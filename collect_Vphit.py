#!/usr/bin/env python3

from lib_common import *

res_dir = f"./{params_path}/Data1D"
mkdir(res_dir)

tmin   = 0
tmax = int(time / dts) + 1
t_range = np.arange(tmin, tmax, 1)

Ja_e = []
Ja_i = []
n_e = []
n_i = []
v_e = []
v_i = []

for t in reduce_array(t_range, rank, proc):
    print(f"{t:5d} [dts]", f"{t * dts / tau:6.3f}", "[tau]")

    comp_Jx_e = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 0)
    comp_Jy_e = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 1)
    _, comp_Ja_e = vx_vy_to_vr_va(comp_Jx_e, comp_Jy_e, COS, SIN)

    comp_Jx_i = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 0)
    comp_Jy_i = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 1)
    _, comp_Ja_i = vx_vy_to_vr_va(comp_Jx_i, comp_Jy_i, COS, SIN)

    comp_n_e = parse_file(get_particles_file("Electrons", "DensPlaneZ", t))
    comp_n_i = parse_file(get_particles_file("Ions", "DensPlaneZ", t))

    comp_Ja_e = phi_averaged(comp_Ja_e, R_MAP)
    comp_Ja_i = phi_averaged(comp_Ja_i, R_MAP)
    comp_n_e = phi_averaged(comp_n_e, R_MAP)
    comp_n_i = phi_averaged(comp_n_i, R_MAP)

    Ja_e.append(comp_Ja_e)
    Ja_i.append(comp_Ja_i)
    n_e.append(comp_n_e)
    n_i.append(comp_n_i)

    v_e.append(np.divide(comp_Ja_e, comp_n_e, out=np.zeros_like(comp_Ja_e), where=(np.abs(comp_n_e) >= 1e-4)))
    v_i.append(np.divide(comp_Ja_i, comp_n_i, out=np.zeros_like(comp_Ja_i), where=(np.abs(comp_n_i) >= 1e-4)))

arr_name = [
    (Ja_e, "Ja_e"),
    (Ja_i, "Ja_i"),
    (n_e,  "n_e"),
    (n_i,  "n_i"),
    (v_e,  "v_e"),
    (v_i,  "v_i"),
]

for arr, name in arr_name:
    gathered_list = comm.gather(arr, root=0)

    if (rank == 0):
        arr = aggregate_array(gathered_list)
        arr = np.reshape(arr, (len(t_range), data_shape[0] // 2))
        np.save(f"{res_dir}/{name}_r_t", arr)

