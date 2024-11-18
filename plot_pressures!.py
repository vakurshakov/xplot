#!/usr/bin/env python3

from lib_common import *

ncols=4
nrows=4

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

bx =  -60 #  boundaries[0] + buff * dx #
ex =  +60 #  boundaries[1] - buff * dx #
by =  bx  #  boundaries[2] + buff * dy #
ey =  ex  #  boundaries[3] - buff * dy #

arg_2d = {
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xlabel": "$x,~c/\\omega_{pe}$",
    "ylabel": "$y,~c/\\omega_{pe}$",
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

arg_1d = {
    "xlabel": "$r,~c/\\omega_{pe}$",
    "xlim": (0, (ex-bx)/2),
    "xticks": np.linspace(0, (ex-bx)/2, 5),
}

offset = 15
tmin = 0 # int(3.245 * tau / dts)
tmax = int(time / dts) + 1
t_range = np.arange(tmin + rank * offset, tmax, proc * offset)[::-1]

rs = np.arange(0, data_shape[0] // 2) * dx

def plot(t, s):
    if not is_correct_timestep(t):
        print(f"Warning! Timestep {t} [dts], {t * dts / tau} is incorrect, it would be skipped.")
        return

    filename = f"{res_dir}/{str(t // offset).zfill(4)}"
    if os.path.exists(filename):
        return

    v_n = 2.0
    v_j = (0.02 if s == "Electrons" else 0.002)
    v_mvv = (0.02 if s == "Ions" else 0.004)
    v_prr = (0.02 if s == "Ions" else 0.004)

    n  = Field(None, subplot(fig, gs, 0, 0), boundaries, unsigned_cmap, (0, v_n))
    n.set_axes_args(title=f"$n_{s[0].lower()}(x, y)$", **arg_2d)

    n_l = Field(None, subplot(fig, gs, 0, 1))
    n_l.set_axes_args(title=f"$\\langle n_{s[0].lower()}(x, y) \\rangle$", ylim=(0, v_n), **arg_1d)

    e = (-1.0) if s == "Electrons" else (+1.0)
    n.data = e * parse_file(get_particles_file(s, f"DensPlaneAvgZ", t))
    n_l.data = phi_averaged(n.data, R_MAP)

    n.draw(add_cbar=True)
    n.draw_info()

    n_l.axes_position.plot(rs, n_l.data)
    n_l.draw_info()

    jx = parse_file(get_particles_file(s, f"CurrentPlaneAvgZ", t), 0)
    jy = parse_file(get_particles_file(s, f"CurrentPlaneAvgZ", t), 1)
    jra = vx_vy_to_vr_va(jx, jy, COS, SIN)

    j_names = ["J_r", "J_{\\phi}"]
    mvv_names = ["\\Pi_{rr}", "\\Pi_{\\phi \\phi}"]
    prr_names = ["P_{rr}", "P_{\\phi \\phi}"]
    for i, (p_val, j_name, mvv_name, prr_name) in enumerate(zip(pressures.values(), j_names, mvv_names, prr_names)):
        j = Field(None, subplot(fig, gs, 1, i*2), boundaries, signed_cmap, (-v_j, +v_j))
        mvv = Field(None, subplot(fig, gs, 2, i*2), boundaries, unsigned_cmap, (0, v_mvv))
        prr = Field(None, subplot(fig, gs, 3, i*2), boundaries, unsigned_cmap, (0, v_prr))

        j_l = Field(None, subplot(fig, gs, 1, i*2+1))
        mvv_l = Field(None, subplot(fig, gs, 2, i*2+1))
        prr_l = Field(None, subplot(fig, gs, 3, i*2+1))

        j.set_axes_args(title=f"${j_name}^{s[0].lower()}(x, y)$", **arg_2d)
        mvv.set_axes_args(title=f"${mvv_name}^{s[0].lower()}(x, y)$", **arg_2d)
        prr.set_axes_args(title=f"${prr_name}^{s[0].lower()}(x, y)$", **arg_2d)

        j_l.set_axes_args(title=f"$\\langle {j_name}^{s[0].lower()}(x, y) \\rangle$", ylim=(-v_j, v_j), **arg_1d)
        mvv_l.set_axes_args(title=f"$\\langle {mvv_name}^{s[0].lower()}(x, y) \\rangle$", ylim=(0, v_mvv), **arg_1d)
        prr_l.set_axes_args(title=f"$\\langle {prr_name}^{s[0].lower()}(x, y) \\rangle$", ylim=(0, v_prr), **arg_1d)

        j.data = jra[i]
        mvv.data = parse_file(get_particles_file(s, f"{p_val}PlaneAvgZ", t))

        v = np.divide(j.data, n.data, out=np.zeros_like(j.data), where=(np.abs(n.data) > 1e-4))
        m = mi_me if s == "Ions" else 1.0
        prr.data = mvv.data - n.data * m * v * v

        j_l.data = phi_averaged(j.data, R_MAP)
        mvv_l.data = phi_averaged(mvv.data, R_MAP)
        prr_l.data = phi_averaged(prr.data, R_MAP)

        for diag in [j, mvv, prr]:
            diag.draw(add_cbar=True)
            diag.draw_info()

        for diag in [j_l, mvv_l, prr_l]:
            diag.axes_position.plot(rs, diag.data)
            diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", x=0.5, y=0.99, size=big, bbox=bbox)

    fig.tight_layout()
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(f"{filename}.png")
    fig.clear()

for s in sorts[::-1]:
    res_dir = f"{params_path}/Pressures_{s[0].lower()}!"
    mkdir(res_dir)

    for t in t_range:
        plot(t, s)
