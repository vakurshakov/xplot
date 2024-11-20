#!/usr/bin/env python3

from plot import *

ncols=2
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.3, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

Prr = Field(None, subplot(fig, gs, 0, 0), boundaries, unsigned_cmap)
Paa = Field(None, subplot(fig, gs, 0, 1), boundaries, unsigned_cmap)

Prr_avg = Field(None, subplot(fig, gs, 1, 0))
Paa_avg = Field(None, subplot(fig, gs, 1, 1))

Prr_l = Field(None, Prr_avg.axes_position)
Paa_l = Field(None, Paa_avg.axes_position)

diag_2d = [Prr, Paa]
diag_1d = [Prr_avg, Paa_avg, Prr_l, Paa_l]

bx = boundaries[0] + buff * dx
ex = boundaries[1] - buff * dx
by = boundaries[2] + buff * dy
ey = boundaries[3] - buff * dy

arg_2d = {
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xlabel": "$x,~c/\\omega_{pi}$",
    "ylabel": "$y,~c/\\omega_{pi}$",
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

arg_1d = {
    "xlabel": "$r,~c/\\omega_{pi}$",
    "xlim": (0, (ex-bx)/2),
    "ylim": (0, 1),
    "xticks": np.linspace(0, (ex-bx)/2, 5),
    "yticks": np.linspace(0, 1, 5),
    "yticklabels": [],
}

res_dir = f"./{params_path}/Other"
mkdir(res_dir)

tmax = int(8 * tau / dts)

for sort in sorts:
    def prr_title(sort):
        return "\\Pi_{rr}^" + str.lower(sort[0])

    def paa_title(sort):
        return "\\Pi_{\\phi \\phi}^" + str.lower(sort[0])

    Prr.set_axes_args(title=f"${prr_title(sort)}(x, y)$", **arg_2d)
    Paa.set_axes_args(title=f"${paa_title(sort)}(x, y)$", **arg_2d)

    Prr_avg.set_axes_args(title=f"$\\langle {prr_title(sort)}(x, y) \\rangle_\\phi$ and ${prr_title(sort)}(x, 0)$", **arg_1d)
    Paa_avg.set_axes_args(title=f"$\\langle {paa_title(sort)}(x, y) \\rangle_\\phi$ and ${paa_title(sort)}(x, 0)$", **arg_1d)

    Prr.data = parse_file(get_particles_file(sort, f"{pressures['Prr']}PlaneAvgZ", tmax))
    Paa.data = parse_file(get_particles_file(sort, f"{pressures['Paa']}PlaneAvgZ", tmax))

    Pmax = np.max(Prr.data)
    Prr.data /= Pmax
    Paa.data /= Pmax

    for diag in diag_2d:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    Prr_avg.data = phi_averaged(Prr.data, R_MAP)
    Paa_avg.data = phi_averaged(Paa.data, R_MAP)

    Prr_l.data = Prr.data[data_shape[1] // 2, data_shape[0] // 2:]
    Paa_l.data = Paa.data[data_shape[1] // 2, data_shape[0] // 2:]

    for diag in diag_1d:
        rs = np.arange(0, diag.data.shape[0]) * dx
        diag.axes_position.plot(rs, diag.data)

        if (diag.axes_args != None):
            diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{tmax * dts / tau:.3f}" "}$", x=0.515, size=big, bbox=bbox)

    fig.tight_layout()
    fig.tight_layout()

    filename = f"{res_dir}/pressures_" + str.lower(sort[0])
    print(filename)

    fig.savefig(f"{filename}.png")

    for diag in diag_1d + diag_2d:
        diag.clear()
