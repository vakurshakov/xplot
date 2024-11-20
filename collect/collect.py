import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from lib_common import *

tmin = 0
tmax = int(time / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1))

res_dir = f"{params_path}/Collection"
mkdir(res_dir)

def center_avg(d, w=5):
    c0 = data_shape["Z"][0] // 2
    c1 = data_shape["Z"][1] // 2
    return np.mean(d[c1-w:c1+w, c0-w:c0+w])

def process_collection(parse, named_arrays, output, shape=None):
    for t in t_range:
        data = parse(find_correct_timestep(t, t_range))

        for d, (_, arr) in zip(data, named_arrays):
            arr.append(d)

    for name, arr in named_arrays:
        gathered_list = comm.gather(arr, root=0)
        if (rank != 0):
            return

        arr = aggregate_array(gathered_list)
        if shape != None:
            arr = np.reshape(arr, shape)

        np.save(output(name), arr)
