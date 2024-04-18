#!/usr/bin/python3

from final_common import *

ncols=3
nrows=3

fig = plt.figure(figsize=(8 * ncols * 1.15, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

ni_l = Field(None, subplot(fig, gs, 0, 2))
pb_l = Field(None, subplot(fig, gs, 1, 2))
je_l = Field(None, subplot(fig, gs, 2, 2))

ni_l.set_axes_args(
    title="${\\rm Ion~density},~n_i$",
    yticks=[0, 0.5, 1.0, 1.5],
    **arg_1d,
)
pb_l.set_axes_args(
    title="${\\rm Pressure~displacement},~1 - b^2$",
    **arg_1d,
)
je_l.set_axes_args(
    title="${\\rm Electric~current},~j_e$",
    xlabel="$(r - r_0) / \\rho_e$",
    xlim=(-70, 120),
    xticks=[0, -30, +30, -60, +60, +90, 120],
    yticks=[-1, -0.75, -0.5, -0.25, 0],
)

t_range = np.arange(1, 5) * int(tau / dts)

hw = 0

for i, t in enumerate(t_range[1:]):
    for j, sort in enumerate(sorts):
        current = Field(None, subplot(fig, gs, i, j), boundaries, signed_cmap, (-0.02, +0.02))

        # We take currents from np_5000 because they are less noisy. It is also
        # valid because we want to compare the result before the m=3 instability
        prefix = get_prefix(t, p5000.restart_timesteps, p5000.prefixes)
        current.data = get_current_data(t, sort, prefix)[1]  # J_phi

        current.set_axes_args(title=f"$J_\\phi^{sort[0].lower()}(x,\,y,\,t = {int(t * dts / tau):d}\,\\tau)$", **arg_2d)
        current.draw(add_cbar=(t == t_range[-1]))
        current.draw_info()
        current.axes_position.set_aspect(1)
        current.axes_position.grid(alpha=0.3)

for t in t_range:
    ni = parse_file(get_particles_file("Ions", "DensPlaneAvgZ", t))
    b  = parse_file(get_fields_file(t), fields.index("Bz"))

    prefix = get_prefix(t, p5000.restart_timesteps, p5000.prefixes)
    je_a = get_current_data(t, "Electrons", prefix)[1]  # J_phi

    ni_l.data = phi_averaged(ni, R_MAP)
    pb_l.data = phi_averaged(b, R_MAP)
    je_l.data = phi_averaged(je_a, R_MAP)

    pb_l.data = 1 - np.square(pb_l.data / B0)

    rs = np.arange(0, data_shape[0] // 2) * dx
    lw = (3 if t == t_range[-1] else 1.5)

    ni_l.axes_position.plot(rs, ni_l.data, linewidth=lw)
    pb_l.axes_position.plot(rs, pb_l.data, linewidth=lw, label=f"$t = {int(t / tau * dts):d}\,\\tau$")
    
    je_min = np.min(je_l.data)
    r0 = np.argmin(je_l.data)
    je_l.data /= np.abs(je_min)
    rho_e = np.sqrt(T_e) / B0
    rs = (np.arange(0, je_l.data.shape[0]) - r0) * dx / rho_e
    je_l.axes_position.plot(rs, je_l.data, linewidth=lw)

    hw_i = np.argwhere(np.abs(je_l.data + 0.5) < 0.1)
    hw += (np.mean(hw_i[np.where(hw_i > r0)]) - np.mean(hw_i[np.where(hw_i < r0)])) * dx / rho_e

    if t == t_range[0]:
        for diag in [ni_l, pb_l, je_l]:
            diag.draw_info()
            diag.axes_position.grid()

    print(t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")

for ax, letter in zip([*fig.axes[:8], fig.axes[9]], "ghi" "ad" "be" "cf" ):
    # annotate_x(ax, "${\\rm " + letter + "}$", 0, 1.06)
    annotate_x(ax, "${\\rm " + letter + "}$", 0.07, 0.07)
            
hw = 32 # hw / 4
x0 = 3.5
je_l.axes_position.arrow(x0, -0.5, +hw/2, 0, **arg_arrow)
je_l.axes_position.arrow(x0, -0.5, -hw/2, 0, **arg_arrow)
je_l.axes_position.text(x0, -0.45, "$32\,\\rho_e$", horizontalalignment="center", fontsize=ssmol)

fig.tight_layout()
fig.tight_layout()

pb_l.axes_position.legend(fontsize=ssmol, loc="upper right")
fig.savefig(f"{res_dir}/figure1.pdf")
