#!/usr/bin/python3

from final_common import *
from final_resonance import *

ncols=3
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.01, 8 * nrows * 1.2))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

t_range = np.arange(3, 9) * int(tau / dts)
r_range = np.arange(0, data_shape[0] // 2) * dx

j_e = Field(None, subplot(fig, gs, 0, 1))
j_e.set_axes_args(
    title="${\\rm Electron~current}, J_{\\phi}^e$",
    ylim=(-0.02, 0.01),
    yticks=[-0.02, -0.01, 0, 0.01],
    **arg_1d
)

j_i = Field(None, subplot(fig, gs, 1, 1))
j_i.set_axes_args(
    title="${\\rm Ion~current}, J_{\\phi}^i$",
    ylim=(-0.002, 0.001),
    yticks=[-0.002, -0.001, 0, 0.001],
    **arg_1d
)

for t in t_range:
    prefix = get_prefix(t, p2000.restart_timesteps, p2000.prefixes)
    j_e.data = phi_averaged(get_current_data(t, "Electrons", prefix)[1], R_MAP)  # J_phi
    j_i.data = phi_averaged(get_current_data(t, "Ions", prefix)[1], R_MAP)  # J_phi

    for diag in [j_e, j_i]:
        ax = diag.axes_position
        ax.plot(r_range, diag.data, label=f"$t = {t * dts / tau:.0f}\,\\tau$")

        if t == t_range[-1]:
            diag.draw_info()
            ax.grid(alpha=0.3)

j_i.axes_position.legend(loc="lower right", fontsize=ssmol * 0.9)

# resonances
Res = get_resonances(subplot(fig, gs, 2, 1), t_range)


r0 = 2.5 # c / w_pi
name = "Ea"
title = names[name][0]

# map
F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"${title}(\\phi,\,t)$", r0, 0.02)

# filtered map
F_mw, w, m = prepare_field_mw(None, f"${title}(m,\,\\omega)$", r0, 50, F_at.data)
F_mw_data = F_mw.data
F_mw_data[:, np.where(np.abs(m) > m0)] = 0.0

F_at_filtered = prepare_field_phit(subplot(fig, gs, 1, 0), name, f"${title}(\\phi,\,t),~m = {m0:0d}$", r0, 0.002)
F_at_filtered.data = inverse_fourier_transform(F_mw_data)[0]

for diag in [F_at, F_at_filtered]:
    diag.draw(add_cbar=True) # , cbar_pad=1.6, cbar_orientation="horizontal")
    diag.draw_info()

# increment
Avg = Field(None, subplot(fig, gs, 2, 0))
Avg.set_axes_args(
    title="${\\rm Increment~of~} m = 3 {\\rm ~mode}$",

    xlabel="$t, \\tau$",
    xlim=(f_tmin, f_tmax),
    xticks=np.linspace(f_tmin, f_tmax, 7),

    # ylabel=f"$\\ln(\\langle |{title}(\\phi,\,t)|^2 \\rangle)$",
    ylim=(-21, -13),
    yticks=np.linspace(-21, -13, 5),
)

ts = np.arange(tmin, tmax) * dts / tau

ff_avg = np.log(np.mean(np.square(F_at_filtered.data), axis=1))
Avg.axes_position.plot(ts, ff_avg)

Avg.draw_info()

# filtered fourier lines
ax = F_at_filtered.axes_position
phi = np.linspace(0, 2 * np.pi)
const = -9.2
ax.plot(phi, (- m0 * phi - const) / w0, linestyle='--', c='r', linewidth=2)
ax.text(np.pi - 0.3, 3.5, f"$\\omega \\approx {w0 / (Omega_i * tau):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)

# increment lines
ax = Avg.axes_position
y0 = -20.2
x0 = f_tmin
gamma = 0.6
ax.plot(ts, -19.9 + 2 * gamma * (ts - f_tmin), linestyle='--', c='r', linewidth=2)
ax.text(5.0, -14, f"$\\Gamma \\approx {(gamma / (Omega_i * tau)):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)

for ax, letter in zip([*fig.axes[:2], fig.axes[2], fig.axes[4], fig.axes[3], fig.axes[-1]], "de" "fba" "c"):
    annotate_x(ax, "${\\rm " + letter + "}$", 0.09, 0.9)

fig.tight_layout()
fig.tight_layout()

fig.savefig(f"{res_dir}/figure4.pdf")
