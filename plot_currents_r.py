#!/usr/bin/python3

from lib_common import *

### Figure parameters.

ncols=3
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.3, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

vJe = (-0.04, +0.04)
vJi = (-0.004, +0.004)
vB = (0, B0)

Je = Field(None, subplot(fig, gs, 0, 0), boundaries, signed_cmap, vJe)
Ji = Field(None, subplot(fig, gs, 1, 0), boundaries, signed_cmap, vJi)
Bz = Field(None, subplot(fig, gs, 2, 0), boundaries, unsigned_cmap.reversed(), vB)

Je_avg = Field(None, subplot(fig, gs, 0, 1), None, None, vJe)
Ji_avg = Field(None, subplot(fig, gs, 1, 1), None, None, vJi)
Bz_avg = Field(None, subplot(fig, gs, 2, 1), None, None, vB)

Je_l = Field(None, Je_avg.axes_position, None, None, vJe)
Ji_l = Field(None, Ji_avg.axes_position, None, None, vJi)
Bz_l = Field(None, Bz_avg.axes_position, None, None, vB)

bx = boundaries[0] + buff * dx
ex = boundaries[1] - buff * dx
by = boundaries[2] + buff * dy
ey = boundaries[3] - buff * dy

titles = []

arg_2d = {
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xlabel": "$x,~c/\\omega_{pi}$",
    "ylabel": "$y,~c/\\omega_{pi}$",
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

arg_comm = {
    "xlabel": "$r,~c/\\omega_{pi}$",
    "xlim": (0, (ex-bx)/2),
    "xticks": np.linspace(0, (ex-bx)/2, 5),
}

arg_B = {
    **arg_comm,
    "ylim": (vB[0] - 0.01, vB[1] + 0.01),
    "yticks": np.linspace(vB[0], vB[1], 5),
}

Je.set_axes_args(**arg_2d, title="$J_r^e(x, y)$")
Ji.set_axes_args(**arg_2d, title="$J_r^i(x, y)$")
Bz.set_axes_args(**arg_2d, title="$B_z(x, y)$")

Je_avg.set_axes_args(title="$\\langle J_r^e(x, y) \\rangle_\\phi$ and $J_r^e(x, 0)$",
    ylim=(vJe[0], vJe[1]),
    yticks=np.linspace(vJe[0], vJe[1], 5),
    **arg_comm,
)

Ji_avg.set_axes_args(title="$\\langle J_r^i(x, y) \\rangle_\\phi$ and $J_r^i(x, 0)$",
    ylim=(vJi[0], vJi[1]),
    yticks=np.linspace(vJi[0], vJi[1], 5),
    **arg_comm,
)

Bz_avg.set_axes_args(title="$\\langle B_z(x, y) \\rangle_\\phi$ and $B_z(x, 0)$", **arg_B)

# Je_l.set_axes_args(**arg_E, title="$J_{\\phi}^e(x, 0)$")
# Ji_l.set_axes_args(**arg_E, title="$J_{\\phi}^i(x, 0)$")
# Bz_l.set_axes_args(**arg_B, title="$B_z(x, 0)$")


res_dir = f"./{params_path}/Currents_r"
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

    comp_x = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 0)
    comp_y = parse_file(get_particles_file("Electrons", "CurrentPlaneAvgZ", t), 1)
    Je.data, _ = vx_vy_to_vr_va(comp_x, comp_y, COS, SIN)

    comp_x = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 0)
    comp_y = parse_file(get_particles_file("Ions", "CurrentPlaneAvgZ", t), 1)
    Ji.data, _ = vx_vy_to_vr_va(comp_x, comp_y, COS, SIN)

    Bz.data = parse_file(get_fields_file(t), fields.index("Bz"))

    for diag in [Je, Ji, Bz]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    Je_avg.data = phi_averaged(Je.data, R_MAP)
    Ji_avg.data = phi_averaged(Ji.data, R_MAP)
    Bz_avg.data = phi_averaged(Bz.data, R_MAP)

    rmax = data_shape[1] // 2
    Je_l.data = Je.data[rmax, rmax:]
    Ji_l.data = Ji.data[rmax, rmax:]
    Bz_l.data = Bz.data[rmax, rmax:]

    for diag in [Je_avg, Ji_avg, Bz_avg, Je_l, Ji_l, Bz_l]:
        rs = np.linspace(0, rmax * dx, diag.data.shape[0])
        diag.axes_position.plot(rs, diag.data)

    for diag in [Je_avg, Ji_avg, Bz_avg]:
        diag.draw_info()


    bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3")
    annotate_x(fig.axes[1], "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(f"{filename}.png")

    for diag in [Je, Ji, Bz, Je_avg, Ji_avg, Bz_avg, Je_l, Ji_l, Bz_l]:
        diag.clear()
