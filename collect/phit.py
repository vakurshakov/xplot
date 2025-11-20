#!/usr/bin/env python3

from collect import *

# NOTE: This diagnostic can be really slow :9(

r0    = int(10 / dx)  # dx units
rmax  = int(15 / dx)  # dx units
rstep = 5             # dimensionless

rmap = []
for r, map in enumerate(R_MAP):
    if (r0 <= r and r < rmax) and (r % rstep == 0):
        rmap.append((r*dx, map))

def parse(t):
    b = get_parsed_field(magnetic_field("Z"), "B", "Z", "z", t)
    er, ea = get_parsed_field(electric_field("Z"), "E", "Z", "", t)
    ni = get_parsed_scalar(particles_field("Ions", "Density", "Z"), t)
    ne = get_parsed_scalar(particles_field("Electrons", "Density", "Z"), t)
    jri, jai = get_parsed_field(particles_field("Ions", "Current", "Z"), "E", "Z", "", t)
    jre, jae = get_parsed_field(particles_field("Electrons", "Current", "Z"), "E", "Z", "", t)
    
    arr = [[] for _ in named_arrays]
    for (_, map) in rmap:
        arr[0].append(b[map])
        arr[1].append(er[map])
        arr[2].append(ea[map])
        arr[3].append(ni[map])
        arr[4].append(ne[map])
        arr[5].append(jri[map])
        arr[6].append(jai[map])
        arr[7].append(jre[map])
        arr[8].append(jae[map])
    return arr

named_arrays = [
    ["b", [[] for _ in rmap]],
    ["er", [[] for _ in rmap]],
    ["ea", [[] for _ in rmap]],
    ["ni", [[] for _ in rmap]],
    ["ne", [[] for _ in rmap]],
    ["jri", [[] for _ in rmap]],
    ["jai", [[] for _ in rmap]],
    ["jre", [[] for _ in rmap]],
    ["jae", [[] for _ in rmap]],
]

def output(name, r):
    return f"{res_dir}/{name}_phit_r={r:.2f}.npy"

tmin = np.inf

if rank == 0:
    for name, _1 in named_arrays:
        for r, _2 in rmap:
            filename = output(name, r)
            
            if os.path.exists(filename):
                tmin = np.min([tmin, np.load(filename, allow_pickle=True, mmap_mode="r").shape[0]])

            print(f"Minimum timestep after {filename} is {tmin}")

# TODO: It would be better to start timers per collection, not globally
tmin = comm.bcast(tmin)
tmin = 0 if tmin == np.inf else int(tmin)
tmax = int(time / dts) + 1
t_range = reduce_array(np.arange(tmin, tmax, 1))

if len(t_range) == 0:
    print("There is nothing to be processed, exiting")
    exit(0)

for t in t_range:
    data = parse(find_correct_timestep(t, t_range))

    for dd, (_, arr) in zip(data, named_arrays):
        for d, a in zip(dd, arr):
            a.append(d)

for name, arr in named_arrays:
    for a, (r, map) in zip(arr, rmap):
        gathered_list = comm.gather(a, root=0)

        if (rank != 0):
            continue

        a = aggregate_array(gathered_list)
        a = np.reshape(a, ((tmax - tmin), len(map[0])))

        filename = output(name, r)
        print(f"Saving {filename}")

        if tmin > 0:
            la = np.load(filename, allow_pickle=True, mmap_mode="r")
            a = np.concat((la[:tmin,:], a), axis=0)

        np.save(filename, a)