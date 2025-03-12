import numpy as np

# Generator of utility maps for cartesian and polar coordinates
def init_XY_RA(boundaries):
    bx = boundaries[0]
    ex = boundaries[1]
    by = boundaries[2]
    ey = boundaries[3]

    if (ex-bx) != (ey-by):
        raise RuntimeError("Boundaries must be of equal size!")

    x0 = (ex-bx)/2 - 0.5

    X = np.zeros((ey-by, ex-bx), dtype=np.float32) +\
        np.fromiter((x-x0 for x in range(bx, ex, 1)), dtype=np.float32)

    # No upside down turn here
    Y = X.transpose()

    return X, Y, np.sqrt(X * X + Y * Y), np.arctan2(X, Y)


# Creates COS, SIN matricies and a radius - (x,y) cell coordinates map
# It's important that cells in RMAP do follow counter-clockwise order
def init_COS_SIN_RMAP(boundaries):
    X, Y, R, A = init_XY_RA(boundaries)

    zero_tolerance = 0.5

    rmax = (boundaries[1] - boundaries[0]) // 2
    rmap = []
    for r in range(0, rmax, 1):
        mask = np.where(np.abs(R - r) <= zero_tolerance)
        sorted_i = np.argsort(A[mask])
        sorted_mask = (mask[0][sorted_i], mask[1][sorted_i])
        rmap.append(sorted_mask)

    R = np.where(np.abs(R) > zero_tolerance, R, 1.0)

    return X/R, Y/R, rmap


def vx_vy_to_vr_va(vx, vy, COS, SIN):
    vr = + COS * vx + SIN * vy
    va = - SIN * vx + COS * vy

    return vr, va


def phi_averaged(data, RMAP):
    result = np.zeros(len(RMAP))

    for i, r in enumerate(RMAP):
        result[i] = data[r].mean()

    return result
