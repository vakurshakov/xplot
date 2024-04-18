#!/usr/bin/python3

from final_common import *

import scipy.integrate as I

ncols=4
nrows=3

fig = plt.figure(figsize=(8 * ncols * 1.18, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

t_range = np.array([1, 4, 6, 9]) * int(tau / dts)  # dts

for i, t in enumerate(t_range):
    prr = Field(None, subplot(fig, gs, i, 0), boundaries, unsigned_cmap, (0, +0.02))
    prr.set_axes_args(title="$\\Pi_{rr}^i(x,\,y)$", **arg_2d)
    prr.data = parse_file(get_particles_file("Ions", f"{pressures['Prr']}PlaneAvgZ", t))

    vmax = 0.002
    j_r = Field(None, subplot(fig, gs, i, 1), boundaries, signed_cmap, (-vmax, +vmax))
    j_r.set_axes_args(title="$J_r^i(x,\,y)$", **arg_2d)
    prefix = get_prefix(t, p2000.restart_timesteps, p2000.prefixes)
    j_r.data = get_current_data(t, "Ions", prefix)[0]

    for diag in [prr, j_r]:
        diag.draw(add_cbar=(t == t_range[-1]))
        diag.draw_info()
        diag.axes_position.set_aspect(1)
        diag.axes_position.grid(alpha=0.3)

    # annotate_x(prr.axes_position, f"$t = {t * dts / tau:.0f}\,\\tau$", 0.06, 0.9, smol, "left")
    annotate_x(prr.axes_position, f"$t = {t * dts / tau:.0f}\,\\tau$", 0.95, 0.9, smol, "right")

# beta calculation
t_range = np.array([1, 2, 4, 6]) * int(tau / dts)  # dts

for i, t in enumerate(t_range):
    prefix = get_prefix(t, p5000.restart_timesteps, p5000.prefixes) if t == t_range[0] \
        else get_prefix(t, p2000.restart_timesteps, p2000.prefixes)

    pr = parse_file(get_particles_file("Ions", f"{pressures['Prr']}PlaneAvgZ", t, prefix))
    pa = parse_file(get_particles_file("Ions", f"{pressures['Paa']}PlaneAvgZ", t, prefix))
    pr += parse_file(get_particles_file("Electrons", f"{pressures['Prr']}PlaneAvgZ", t, prefix))
    pa += parse_file(get_particles_file("Electrons", f"{pressures['Paa']}PlaneAvgZ", t, prefix))
    b   = parse_file(get_fields_file(t, prefix), fields.index("Bz"))

    pr = phi_averaged(pr, R_MAP) / ((B0 * B0) / 2)
    pa = phi_averaged(pa, R_MAP) / ((B0 * B0) / 2)
    pb = 1 - np.square(phi_averaged(b, R_MAP) / B0)

    rs = np.arange(0, data_shape[0] // 2) * dx
    pd = I.cumtrapz((pr - pa) / (rs + 0.1), dx=dx, initial=0)
    pd = -(pd - pd[-1])
    ps = pr - pd

    beta = Field(None, subplot(fig, gs, i, 2), boundaries, signed_cmap, None)
    beta.set_axes_args(
        title="${\\rm Plasma~beta}$",
        ylim=(-0.1, +1.25),
        **arg_1d
    )

    ax = beta.axes_position
    ax.plot(rs, pr, label="$\\beta_r$", linewidth=2)
    ax.plot(rs, pa, label="$\\beta_{\\phi}$", linewidth=2)
    ax.plot(rs, pd, label="$\\Delta \\beta$", linewidth=2)
    ax.plot(rs, ps, label="$\\beta_r - \\Delta \\beta$", linewidth=2)
    ax.plot(rs, pb, label="$\\beta_b$", linewidth=2)

    ax.grid(alpha=0.3)
    beta.draw_info()

    if (t == t_range[0]):
        ax.legend(loc="center right", fontsize=smol)

    annotate_x(ax, f"$t = {t * dts / tau:.0f}\,\\tau$", 0.95, 0.9, smol, "right")

fig.tight_layout()
fig.tight_layout()

# p_l.axes_position.legend(fontsize=ssmol, loc="upper right")
fig.savefig(f"{res_dir}/figure3.pdf")
