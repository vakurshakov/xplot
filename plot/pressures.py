#!/usr/bin/env python3

from plot import *

ncols=4
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

vPi = (0, +0.03)
vPe = (0, +0.004)

Prr_i = Field(None, subplot(fig, gs, 0, 0), boundaries, unsigned_cmap, vPi)
Prr_e = Field(None, subplot(fig, gs, 2, 0), boundaries, unsigned_cmap, vPe)
Paa_i = Field(None, subplot(fig, gs, 0, 1), boundaries, unsigned_cmap, vPi)
Paa_e = Field(None, subplot(fig, gs, 2, 1), boundaries, unsigned_cmap, vPe)

Prr_i_avg = Field(None, subplot(fig, gs, 1, 0), None, None, vPi)
Prr_e_avg = Field(None, subplot(fig, gs, 3, 0), None, None, vPe)
Paa_i_avg = Field(None, subplot(fig, gs, 1, 1), None, None, vPi)
Paa_e_avg = Field(None, subplot(fig, gs, 3, 1), None, None, vPe)

Prr_i_l = Field(None, Prr_i_avg.axes_position, None, None, vPi)
Prr_e_l = Field(None, Prr_e_avg.axes_position, None, None, vPe)
Paa_i_l = Field(None, Paa_i_avg.axes_position, None, None, vPi)
Paa_e_l = Field(None, Paa_e_avg.axes_position, None, None, vPe)

diag_2d = [Prr_i, Prr_e, Paa_i, Paa_e]
diag_1d = [Prr_i_avg, Prr_e_avg, Paa_i_avg, Paa_e_avg, Prr_i_l, Prr_e_l, Paa_i_l, Paa_e_l]

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
    "xticks": np.linspace(0, (ex-bx)/2, 5),
    "yticklabels": [],
}

arg_e = {
    **arg_1d,
    "ylim": (vPe[0], vPe[1]),
    "yticks": np.linspace(vPe[0], vPe[1], 5),
}

arg_i = {
    **arg_1d,
    "ylim": (vPi[0], vPi[1]),
    "yticks": np.linspace(vPi[0], vPi[1], 5),
}

Prr_i.set_axes_args(title="$\\Pi_{rr}^i(x, y)$", **arg_2d)
Prr_e.set_axes_args(title="$\\Pi_{rr}^e(x, y)$", **arg_2d)
Paa_i.set_axes_args(title="$\\Pi_{\\phi\\phi}^i(x, y)$", **arg_2d)
Paa_e.set_axes_args(title="$\\Pi_{\\phi\\phi}^e(x, y)$", **arg_2d)

Prr_i_avg.set_axes_args(title="$\\langle \\Pi_{rr}^i(x, y) \\rangle_\\phi$ and $\\Pi_{rr}^i(x, 0)$", **arg_i)
Prr_e_avg.set_axes_args(title="$\\langle \\Pi_{rr}^e(x, y) \\rangle_\\phi$ and $\\Pi_{rr}^e(x, 0)$", **arg_e)
Paa_i_avg.set_axes_args(title="$\\langle \\Pi_{\\phi\\phi}^i(x, y) \\rangle_\\phi$ and $\\Pi_{\\phi\\phi}^i(x, 0)$", **arg_i)
Paa_e_avg.set_axes_args(title="$\\langle \\Pi_{\\phi\\phi}^e(x, y) \\rangle_\\phi$ and $\\Pi_{\\phi\\phi}^e(x, 0)$", **arg_e)

res_dir = f"./{params_path}/Pressures"
mkdir(res_dir)

offset = 15
tmin = int(3.245 * tau / dts)
tmax = int(time / dts) + 1
t_range = np.arange(tmin + rank * offset, tmax, proc * offset)


for t in t_range:
    if not is_correct_timestep(t):
        print(f"Warning! Timestep {t} [dts], {t * dts / tau} is incorrect, it would be skipped.")
        continue

    filename = f"{res_dir}/{str(t // offset).zfill(4)}"
    if os.path.exists(filename):
        continue

    Prr_i.data = parse_file(get_particles_file("Ions", f"{pressures['Prr']}PlaneAvgZ", t))
    Paa_i.data = parse_file(get_particles_file("Ions", f"{pressures['Paa']}PlaneAvgZ", t))
    Prr_e.data = parse_file(get_particles_file("Electrons", f"{pressures['Prr']}PlaneAvgZ", t))
    Paa_e.data = parse_file(get_particles_file("Electrons", f"{pressures['Paa']}PlaneAvgZ", t))

    for diag in diag_2d:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    Prr_i_avg.data = phi_averaged(Prr_i.data, R_MAP)
    Paa_i_avg.data = phi_averaged(Paa_i.data, R_MAP)
    Prr_e_avg.data = phi_averaged(Prr_e.data, R_MAP)
    Paa_e_avg.data = phi_averaged(Paa_e.data, R_MAP)

    Prr_i_l.data = Prr_i.data[data_shape[1] // 2, data_shape[0] // 2:]
    Paa_i_l.data = Paa_i.data[data_shape[1] // 2, data_shape[0] // 2:]
    Prr_e_l.data = Prr_e.data[data_shape[1] // 2, data_shape[0] // 2:]
    Paa_e_l.data = Paa_e.data[data_shape[1] // 2, data_shape[0] // 2:]

    for diag in diag_1d:
        rs = np.arange(0, diag.data.shape[0]) * dx
        diag.axes_position.plot(rs, diag.data)

        if (diag.axes_args != None):
            diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", x=0.515, size=big, bbox=bbox)

    fig.tight_layout()
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(f"{filename}.png")

    for diag in diag_1d + diag_2d:
        diag.clear()
