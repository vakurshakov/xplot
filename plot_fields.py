#!/usr/bin/python3

from lib_common import *

ncols=3
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.3, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

vE = (-0.02, +0.02)
vB = (0, B0)

Er = Field(None, subplot(fig, gs, 0, 0), boundaries, signed_cmap, vE)
Ea = Field(None, subplot(fig, gs, 1, 0), boundaries, signed_cmap, vE)
Bz = Field(None, subplot(fig, gs, 2, 0), boundaries, unsigned_cmap.reversed(), vB)

Er_avg = Field(None, subplot(fig, gs, 0, 1), None, None, vE)
Ea_avg = Field(None, subplot(fig, gs, 1, 1), None, None, vE)
Bz_avg = Field(None, subplot(fig, gs, 2, 1), None, None, vB)

Er_l = Field(None, Er_avg.axes_position, None, None, vE)
Ea_l = Field(None, Ea_avg.axes_position, None, None, vE)
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

arg_E = {
    **arg_comm,
    "ylim": (vE[0], vE[1]),
    "yticks": np.linspace(vE[0], vE[1], 5),
}

arg_B = {
    **arg_comm,
    "ylim": (vB[0] - 0.01, vB[1] + 0.01),
    "yticks": np.linspace(vB[0], vB[1], 5),
}

Er.set_axes_args(**arg_2d, title="$E_r(x, y)$")
Ea.set_axes_args(**arg_2d, title="$E_{\\phi}(x, y)$")
Bz.set_axes_args(**arg_2d, title="$B_z(x, y)$")

Er_avg.set_axes_args(**arg_E, title="$\\langle E_r(x, y) \\rangle_\\phi$ and $E_r(x, 0)$")
Ea_avg.set_axes_args(**arg_E, title="$\\langle E_{\\phi}(x, y) \\rangle_\\phi$ and $E_{\\phi}(x, 0)$")
Bz_avg.set_axes_args(**arg_B, title="$\\langle B_z(x, y) \\rangle_\\phi$ and $B_z(x, 0)$")

# Er_l.set_axes_args(**arg_E, title="$E_r(x, 0)$")
# Ea_l.set_axes_args(**arg_E, title="$E_{\\phi}(x, 0)$")
# Bz_l.set_axes_args(**arg_B, title="$B_z(x, 0)$")


res_dir = f"./{params_path}/Fields"
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

    # Er, Ea
    comp_x = parse_file(get_fields_file(t), fields.index("Ex"))
    comp_y = parse_file(get_fields_file(t), fields.index("Ey"))
    Er.data, Ea.data = vx_vy_to_vr_va(comp_x, comp_y, COS, SIN)

    Bz.data = parse_file(get_fields_file(t), fields.index("Bz"))

    for diag in [Er, Ea, Bz]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

        # r_range = np.arange(2, 3, 0.1)
        # diag.axes_position.scatter(r_range, np.zeros_like(r_range), s=10, color='black')
    
    Er_avg.data = phi_averaged(Er.data, R_MAP)
    Ea_avg.data = phi_averaged(Ea.data, R_MAP)
    Bz_avg.data = phi_averaged(Bz.data, R_MAP)
    
    Er_l.data = Er.data[data_shape[1] // 2, (data_shape[0] // 2):]
    Ea_l.data = Ea.data[data_shape[1] // 2, (data_shape[0] // 2):]
    Bz_l.data = Bz.data[data_shape[1] // 2, (data_shape[0] // 2):]

    for diag in [Er_avg, Ea_avg, Bz_avg, Er_l, Ea_l, Bz_l]:
        rs = np.arange(0, diag.data.shape[0]) * dx
        diag.axes_position.plot(rs, diag.data)

    for diag in [Er_avg, Ea_avg, Bz_avg]:
        diag.draw_info()

    bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3")
    annotate_x(fig.axes[1], "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(f"{filename}.png")

    for diag in [Er, Ea, Bz, Er_avg, Ea_avg, Bz_avg, Er_l, Ea_l, Bz_l]:
        diag.clear()
