#!/usr/bin/python3

from lib_common import *

res_dir = f"./{params_path}/Data2D"
mkdir(res_dir)

tmin = int(0 * tau / dts)
tmax = int(8 * tau / dts)

t_range = np.arange(tmin, tmax, 1)
t_range = reduce_array(t_range, rank, proc)

Er_xyt = np.ndarray((len(t_range), data_shape[0], data_shape[1]), dtype=np.float32)
Ea_xyt = np.ndarray((len(t_range), data_shape[0], data_shape[1]), dtype=np.float32)

for t in t_range:
    print(f"{t:5d} [dts]", f"{t * dts / tau:6.3f}", "[tau]")

    # TODO: Create some fourier helper
    comp_Ex = parse_file(get_fields_file(t), fields.index("Ex"))
    comp_Ey = parse_file(get_fields_file(t), fields.index("Ey"))
    comp_Er, comp_Ea = vx_vy_to_vr_va(comp_Ex, comp_Ey, COS, SIN)

    Er_xyt[t - t_range[0]] = comp_Er
    Ea_xyt[t - t_range[0]] = comp_Ea

arr_name = [
    (Er_xyt, "Er"),
    (Ea_xyt, "Ea"),
]

for arr, name in arr_name:
    for y in range(data_shape[0]):
        tmp = arr[:,y,:]
        gathered_list = comm.gather(tmp, root=0)

        if (rank == 0):
            tmp = aggregate_array(gathered_list)
            tmp = np.reshape(tmp, (tmax - tmin, data_shape[1]))
            np.save(f"{res_dir}/{name}_xt_y={y}", tmp)
