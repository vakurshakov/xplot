from lib_common import *

tmin = 0
tmax = int(time / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1))

def center_avg(d, w=5):
    c0 = data_shape["Z"][0] // 2
    c1 = data_shape["Z"][1] // 2
    return np.mean(d[c1-w:c1+w, c0-w:c0+w])


def process_collection(parse_data, named_arrays, output_file, shape=None):
    for t in t_range:
        data = parse_data(find_correct_timestep(t))

        for d, (_, arr) in zip(data, named_arrays):
            arr.append(d)

    for name, arr in named_arrays:
        gathered_list = comm.gather(arr, root=0)

        if (rank == 0):
            arr = aggregate_array(gathered_list)
            if shape != None:
                arr = np.reshape(arr, shape)
            np.save(output_file(name), arr)
