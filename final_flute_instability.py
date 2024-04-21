#!/usr/bin/python3

from final_common import *
from plot_fourier_common import *

import scipy.integrate as I

ncols=3
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1, 8 * nrows * 1.2))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1, 0.8], figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

t_range = np.array([4, 5, 6, 7, 8]) * int(tau / dts)
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

j_i.axes_position.legend(loc="lower right", fontsize=ssmol)

r0 = 2.5 # c / w_pi
name = "Ea"
title = names[name][0]

# map
F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"${title}(\\phi,\,t)$", r0, 0.02)

# fourier
F_mw, w, m = prepare_field_mw(subplot(fig, gs, 2, 0), f"${title}(m,\,\\omega)$", r0, 50, F_at.data)
F_mw_data = F_mw.data
F_mw.data = np.abs(F_mw.data)

# filtered map
F_mw_data[:, np.where(np.abs(m) > m0)] = 0.0

F_at_filtered = prepare_field_phit(subplot(fig, gs, 1, 0), name, f"${title}(\\phi,\,t),~m = {m0:0d}$", r0, 0.002)
F_at_filtered.data = inverse_fourier_transform(F_mw_data)[0]

for diag in [F_at, F_mw, F_at_filtered]:
    diag.draw(add_cbar=True, cbar_pad=1.6, cbar_orientation="horizontal")
    diag.draw_info()

    # ax = diag.axes_position
    # annotate_x(ax, f"$r = {r0 * np.sqrt(mi_me):.0f}~c/\\omega_""{pe}""$", 0.1, 0.06, ha="left")

# increment
Avg = Field(None, subplot(fig, gs, 2, 1))
Avg.set_axes_args(
    title="${\\rm Increment~of~} m = 3 {\\rm ~mode}$",

    xlabel="$t, \\tau$",
    xlim=(f_tmin, f_tmax),
    xticks=np.linspace(f_tmin, f_tmax, 7),

    ylabel=f"$\\ln(\\langle |{title}(\\phi,\,t)|^2 \\rangle)$",
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

# fourier lines
ax = F_mw.axes_position
ax.plot([-m0, -m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=1)
ax.plot([+m0, +m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=1)

# increment lines
ax = Avg.axes_position
y0 = -20.2
x0 = f_tmin
gamma = 0.6
ax.plot(ts, -19.9 + 2 * gamma * (ts - f_tmin), linestyle='--', c='r', linewidth=2)
ax.text(5.0, -14, f"$\\Gamma \\approx {(gamma / (Omega_i * tau)):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)

# modified Simon-Hoh instability lines
'''
t0 = int(4 * tau / dts)
tw = int(0.5 * tau / dts)
tmin_avg = t0 - tw
tmax_avg = t0 + tw

e_r = np.load(f"{params_path}/Data1D/Er_phi_t_r={r0:.2f}.npy", allow_pickle=True)
b_z = np.load(f"{params_path}/Data1D/Bz_phi_t_r={r0:.2f}.npy", allow_pickle=True)

n_i = parse_file(get_particles_file("Ions", f"DensPlaneAvgZ", t))
n_i = phi_averaged(n_i, R_MAP)

nr0 = int(r0 * np.sqrt(mi_me) / dx)
n0_i = n_i[nr0]
dn_i = np.gradient(n_i, dx)[nr0]

e_r = np.mean(e_r[tmin_avg:tmax_avg,:])
b_z = np.mean(b_z[tmin_avg:tmax_avg,:])

omega_E = -(m0 / (r0 * np.sqrt(mi_me))) * e_r / b_z
omega_s = -(m0 / (r0 * np.sqrt(mi_me))) * (T_e / b_z) * (dn_i / n0_i)
gamma_sh = np.sqrt(T_e / mi_me) * (m0 / (r0 * np.sqrt(mi_me))) * np.sqrt(omega_E / omega_s)
ax.plot(ts, -19.9 + 2 * gamma_sh * ts * tau, linestyle='--', c='r', linewidth=2)

print("Er:", e_r, "Bz:", b_z, "n0:", n0_i, "dn0:", dn_i, "gamma_sh:", gamma_sh)
'''

for ax, letter in zip([*fig.axes[:2], fig.axes[2], fig.axes[4], fig.axes[3], fig.axes[-1]], "de" "abc" "f"):
    annotate_x(ax, "${\\rm " + letter + "}$", 0.09, 0.06)

fig.tight_layout()
fig.tight_layout()

# p_l.axes_position.legend(fontsize=ssmol, loc="upper right")
fig.savefig(f"{res_dir}/figure4.pdf")
