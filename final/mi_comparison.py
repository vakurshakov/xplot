#!/usr/bin/env python3

from final_common import *

import scipy.integrate as I

ncols=3
nrows=3

fig = plt.figure(figsize=(8 * ncols * 1.15, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_titlesize(42)
set_labelsize(40)
set_ticksize(36)

def draw_current_lines():
    ax = j_l.axes_position

    rho_i = np.sqrt(mi_me * T_i / T_e) # hw / 4
    x0 = 0
    y0 = -0.5
    text_offset = 0.05
    ax.arrow(x0, y0, +rho_i/2, 0, **arg_arrow)
    ax.arrow(x0, y0, -rho_i/2, 0, **arg_arrow)
    ax.text(x0, y0 + text_offset, "$32\,\\rho_e$", horizontalalignment="center", fontsize=ticksize)

    x0 = 70
    y0 = 0.04
    ax.arrow(x0, y0, +1 * rho_i, 0, **arg_arrow)
    ax.arrow(x0, y0, -1 * rho_i, 0, **arg_arrow)
    ax.text(x0, y0 + text_offset, "$2\,\\rho_i$", horizontalalignment="center", fontsize=ticksize)

    xlim = (40, -1)
    ax.fill_between(rc[xlim[0]:xlim[1]], np.zeros_like(j)[xlim[0]:xlim[1]], j[xlim[0]:xlim[1]], color="grey", alpha=0.2, zorder=-1)
    ax.annotate(
        "${\\rm drift~instability~effect}$",
        (72, -0.085),
        (46, -0.320),
        **arg_anno,
    )


def draw_magnetic_lines():
    ax = b_l.axes_position

    ylim = ax.get_ylim()
    ax.fill_betweenx(ylim, [0, 0], [15, 15], color="grey", alpha=0.2, zorder=-1)
    ax.annotate(
        "${\\rm injection~region}$",
        (16, 0.1),
        (30, 0.14),
        **arg_anno,
    )
    ax.set_ylim(ylim)

t0_mi_100 = 7_632  / dts
t0_mi_1   = 13_998 / dts

for i, (t, params) in enumerate(zip([t0_mi_1, t0_mi_100], [pmi1, p2000])):
    vmax = 0.01
    e_r = Field(None, subplot(fig, gs, 0, i), boundaries, signed_cmap, (-vmax, +vmax))
    e_r.set_axes_args(title=f"$E_r(x,\,y)$", **arg_2d)

    e_a = Field(None, subplot(fig, gs, 1, i), boundaries, signed_cmap, (-vmax, +vmax))
    e_a.set_axes_args(title="$E_{\\phi}(x,\,y)$", **arg_2d)

    j_a = Field(None, subplot(fig, gs, 2, i), boundaries, signed_cmap, (-vmax, +vmax))
    j_a.set_axes_args(title="$J_{\\phi}^e(x,\,y)$", **arg_2d)

    prefix = get_prefix(t, params.restart_timesteps, params.prefixes)
    e_r.data, e_a.data = get_electric_fields_data(t, prefix)
    j_a.data = get_current_data(t, "Electrons", prefix)[1] # J_phi

    for diag in [e_r, e_a, j_a]:
        diag.draw(add_cbar=(diag == j_a))
        diag.draw_info()
        diag.axes_position.set_aspect(1)
        diag.axes_position.grid(alpha=0.3)


b_l = Field(None, subplot(fig, gs, 0, 2))
p_l = Field(None, subplot(fig, gs, 1, 2))
j_l = Field(None, subplot(fig, gs, 2, 2))

b_l.set_axes_args(
    title="${\\rm Magnetic~field},~B_z$",
    **arg_1d,
)
p_l.set_axes_args(
    title="${\\rm Electric~potential},~e \\varphi / T_i$",
    xlabel="$r,~c/\\omega_{pe}$",
    xlim=(0, 100),
)
j_l.set_axes_args(
    title="${\\rm Electric~current},~J_{\\phi}$",
    xlabel="$(r - r_0) / \\rho_e$",
    xlim=(-60, 120),
    xticks=[0, -30, +30, -60, +60, +90, 120],
    yticks=[-1, -0.75, -0.5, -0.25, 0],
)

jpmi1 = None
jp2000 = None

for t0, params, label in zip([t0_mi_1, t0_mi_100], [pmi1, p2000], ["$\\rho_i = \\rho_e$", "$\\rho_i \\approx 32 \\rho_e$"]):
    prefix = get_prefix(t0, params.restart_timesteps, params.prefixes)
    b = get_magnetic_field_data(t0, prefix)
    e = get_electric_fields_data(t0, prefix)[0] # E_r
    j =  get_current_data(t0, "Electrons", prefix)[1] + get_current_data(t0, "Ions", prefix)[1] # J_phi, total!

    b = phi_averaged(b, R_MAP)
    e = phi_averaged(e, R_MAP)
    j = phi_averaged(j, R_MAP)

    rs = np.arange(0, data_shape[0] // 2) * dx

    p = - I.cumtrapz(e, rs, dx, initial=0) / T_i
    je_min = np.min(j)
    j /= np.abs(je_min)
    rho_e = np.sqrt(T_e) / B0

    r0 = np.argmin(j)
    rc = (np.arange(0, j.shape[0]) - r0) * dx / rho_e

    if (params == pmi1):
        jpmi1 = j
    elif (params == p2000):
        jp2000 = j

    lw = 3
    b_l.axes_position.plot(rs, b, linewidth=lw, label=label)
    p_l.axes_position.plot(rs, p, linewidth=lw, label=label)
    j_l.axes_position.plot(rc, j, linewidth=lw, label=label)

for diag in [b_l, p_l, j_l]:
    diag.draw_info()
    diag.axes_position.grid(alpha=0.3)

for ax, letter in zip([*fig.axes[:3], *fig.axes[4:7], *fig.axes[8:]], "abcdefghi"):
    annotate_x(ax, "${\\rm " + letter + "}$", 0.92, 0.06)

for ax, text in zip([fig.axes[0], fig.axes[4]], ["$m_i = 1,~\\rho_i = \\rho_e$", "$m_i = 100,~\\rho_i \\approx 32 \\rho_e$"]):
    annotate_x(ax, text, 0.06, 0.9, ticksize, "left")

draw_current_lines()
draw_magnetic_lines()

b_l.axes_position.legend(bbox_to_anchor=(0.47, 0.4), fontsize=ticksize, framealpha=1)
p_l.axes_position.legend(bbox_to_anchor=(0.47, 0.4), fontsize=ticksize, framealpha=1)
j_l.axes_position.legend(bbox_to_anchor=(0.47, 0.4), fontsize=ticksize, framealpha=1)

fig.tight_layout()
fig.tight_layout()

# p_l.axes_position.legend(fontsize=ticksize, loc="upper right")
fig.savefig(f"{res_dir}/figure2.pdf")
