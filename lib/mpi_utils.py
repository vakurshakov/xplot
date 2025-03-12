import numpy as np

from mpi4py import MPI

# MPI utilities

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
proc = comm.Get_size()

def create_t_range(tmin, tmax, offset):
    return np.arange(tmin + rank * offset, tmax + 1, proc * offset)

def reduce_array(a):
    length = len(a)
    chunk = length // proc
    remain = length % proc
    prev = 0

    if rank < remain:
        prev = (chunk + 1) * rank
        chunk += 1
    else:
        prev = chunk * rank + remain

    return a[prev: prev + chunk]

def aggregate_array(a):
    result = np.array([])
    for s in a:
        result = np.append(result, s)

    return result
