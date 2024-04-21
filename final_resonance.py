#!/usr/bin/python3

from final_common import *
from plot_fourier_common import *

import scipy.integrate as I

ncols=1
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

t_range = np.arange(3, 9) * int(tau / dts)
r_range = np.arange(0, data_shape[0] // 2) * dx

plot = Field(None, subplot(fig, gs, 0, 0))

ex = 30
plot.set_axes_args(
    ylabel="$\\omega / \\Omega_i$",
    ylim=(-0.05, 0.5),
    yticks=np.linspace(0, 0.5, 6),

    xlabel="$r,~c/\\omega_{pe}$",
    xlim=(0, ex),
    xticks=np.linspace(0, ex, 4),
)

ax = plot.axes_position

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for i, t in enumerate(t_range):
    prefix = get_prefix(t, p2000.restart_timesteps, p2000.prefixes)
    b_z = phi_averaged(get_magnetic_field_data(t, prefix), R_MAP)
    j_i = phi_averaged(get_current_data(t, "Ions", prefix)[1], R_MAP)
    n_i = phi_averaged(parse_file(get_particles_file("Ions", "DensPlaneAvgZ", t)), R_MAP)

    j_i = np.divide(j_i, (r_range + 0.01))
    omega_i = np.divide(j_i, n_i, out=np.zeros_like(j_i), where=(np.abs(n_i) > 1e-3))

    ax.plot(r_range, -m0 * omega_i / Omega_i, c=colors[i], linewidth=2)
    ax.plot(r_range, ((b_z / mi_me) / Omega_i), c=colors[i], linewidth=2, label=f"$t = {t * dts / tau:.0f}\,\\tau$")

plot.draw_info()

ax.legend(loc="upper left", fontsize=ssmol * 0.73)
ax.grid(alpha=0.3)

ax.annotate(
    "$\\Omega_i$",
    (15, 0.42),  # point to annotate
    (12, 0.45),  # text location
    **arg_anno,
)

ax.annotate(
    "$k V_d^i$",
    (21, 0.03),  # point to annotate
    (26, 0.1),   # text location
    **arg_anno,
)

fig.tight_layout()
fig.tight_layout()

fig.savefig(f"{res_dir}/figure41.pdf")
