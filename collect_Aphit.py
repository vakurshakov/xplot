#!/usr/bin/python3

from lib_common import *

res_dir = f"./{params_path}/Data1D"
mkdir(res_dir)

x0   = int(1.0 / dx)
xmax = int(4.0 / dx)
xoff = 5

rmap = []
for i, map in enumerate(R_MAP):
    if (i % xoff == 0) and (x0 <= i and i < xmax):
        rmap.append((i*dx, map))

tmin   = 0
tmax = int(time / dts) + 1
t_range = np.arange(tmin, tmax, 1)

for r, map in rmap:
    Er = []
    Ea = []
    Bz = []
    Jr_e = []
    Ja_e = []
    Jr_i = []
    Ja_i = []

    for t in reduce_array(t_range, rank, proc):
        print(f"r: {r}, phi_n: {len(map[0])}", f"{t:5d} [dts]", f"{t * dts / tau:6.3f}", "[tau]")

        # TODO: Create some fourier helper
        comp_Bz = parse_file(get_fields_file(t), fields.index("Bz"))

        comp_Ex = parse_file(get_fields_file(t), fields.index("Ex"))
        comp_Ey = parse_file(get_fields_file(t), fields.index("Ey"))
        comp_Er, comp_Ea = vx_vy_to_vr_va(comp_Ex, comp_Ey, COS, SIN)

        comp_Jx_e = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 0)
        comp_Jy_e = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 1)
        comp_Jr_e, comp_Ja_e = vx_vy_to_vr_va(comp_Jx_e, comp_Jy_e, COS, SIN)

        comp_Jx_i = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 0)
        comp_Jy_i = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 1)
        comp_Jr_i, comp_Ja_i = vx_vy_to_vr_va(comp_Jx_i, comp_Jy_i, COS, SIN)

        Bz.append(comp_Bz[map])
        Er.append(comp_Er[map])
        Ea.append(comp_Ea[map])
        Jr_e.append(comp_Jr_e[map])
        Ja_e.append(comp_Ja_e[map])
        Jr_i.append(comp_Jr_i[map])
        Ja_i.append(comp_Ja_i[map])

    arr_name = [
        (Bz,   "Bz"),
        (Er,   "Er"),
        (Ea,   "Ea"),
        (Jr_e, "Jr_e"),
        (Ja_e, "Ja_e"),
        (Jr_i, "Jr_i"),
        (Ja_i, "Ja_i"),
    ]

    for arr, name in arr_name:
        gathered_list = comm.gather(arr, root=0)

        if (rank == 0):
            arr = aggregate_array(gathered_list)
            arr = np.reshape(arr, (len(t_range), len(map[0])))
            np.save(f"{res_dir}/{name}_phi_t_r={r:.2f}", arr)

