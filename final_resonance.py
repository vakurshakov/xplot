#!/usr/bin/python3

from final_common import *
from plot_fourier_common import *

arg_resonance = {
    "title": "${\\rm Drift~resonances}$",

    "ylabel": "$\\omega / \\Omega_v^i$",
    "ylim": (-0.05, 0.5),
    "yticks": np.linspace(0, 0.5, 6),

    "xlabel": "$r,~c/\\omega_{pe}$",
    "xlim": (0, 30),
    "xticks": np.linspace(0, 30, 4),
}

def get_resonances(ax, t_range):
    plot = Field(None, ax)
    plot.set_axes_args(**arg_resonance)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    r_range = np.arange(0, data_shape[0] // 2) * dx
    for i, t in enumerate(t_range):
        prefix = get_prefix(t, p2000.restart_timesteps, p2000.prefixes)
        b_z = phi_averaged(get_magnetic_field_data(t, prefix), R_MAP)
        j_i = phi_averaged(get_current_data(t, "Ions", prefix)[1], R_MAP)
        n_i = phi_averaged(parse_file(get_particles_file("Ions", "DensPlaneAvgZ", t)), R_MAP)

        j_i = np.divide(j_i, (r_range + 0.01))
        omega_i = np.divide(j_i, n_i, out=np.zeros_like(j_i), where=(np.abs(n_i) > 1e-3))

        ax.plot(r_range, -m0 * omega_i / Omega_i, c=colors[i], linewidth=2)
        ax.plot(r_range, ((b_z / mi_me) / Omega_i), c=colors[i], linewidth=2, label=f"$t = {t * dts / tau:.0f}\,\\tau$")

    ax.annotate(
        "$\\Omega_i$",
        (15, 0.42),  # point to annotate
        (9, 0.45),  # text location
        **arg_anno,
    )
    ax.annotate(
        "$k V_d^i$",
        (21, 0.03),  # point to annotate
        (24, 0.1),   # text location
        **arg_anno,
    )
    plot.draw_info()
    return plot