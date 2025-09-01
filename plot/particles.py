#!/usr/bin/env python3

from plot import *


def plot_particles3(s, t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    jrs[0].data = get_parsed_field(jrs[0], "E", "X", "y", t)
    jas[0].data = get_parsed_field(jas[0], "E", "X", "x", t)

    jrs[1].data = get_parsed_field(jrs[1], "E", "Y", "x", t)
    jas[1].data = get_parsed_field(jas[1], "E", "Y", "y", t)

    jrs[2].data, jas[2].data = get_parsed_field(jrs[2], "E", "Z", "", t)

    for i, plane in enumerate(planes):
        e = (-1) if s == "Electrons" else (+1)
        ns[i].data = e * get_parsed_scalar(ns[i], t)
        jzs[i].data = get_parsed_field(jzs[i], "E", plane, "z", t)

    for i, diag in enumerate(jrs + jas + jzs + ns):
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=0.99, bbox=bbox, fontsize=big)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(filename)

    for diag in jrs + jas + ns:
        diag.clear()


if __name__ == "__main__":
    for s in sorts:
        ncols=4
        nrows=3

        fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
        gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.2, 1.2, 1], figure=fig)

        v = 4.0
        ns = []
        ns.append(particles_field(s, "Density", "X", subplot(fig, gs, 0, 0), f"$n_{s[0].lower()}$", (0, v), unsigned_cmap))
        ns.append(particles_field(s, "Density", "Y", subplot(fig, gs, 0, 1), f"$n_{s[0].lower()}$", (0, v), unsigned_cmap))
        ns.append(particles_field(s, "Density", "Z", subplot(fig, gs, 0, 2), f"$n_{s[0].lower()}$", (0, v), unsigned_cmap))

        v = 0.02
        jrs = []
        jrs.append(particles_field(s, "Current", "X", subplot(fig, gs, 1, 0), "$J_r^"f"{s[0].lower()}""$", (-v, +v)))
        jrs.append(particles_field(s, "Current", "Y", subplot(fig, gs, 1, 1), "$J_r^"f"{s[0].lower()}""$", (-v, +v)))
        jrs.append(particles_field(s, "Current", "Z", subplot(fig, gs, 1, 2), "$J_r^"f"{s[0].lower()}""$", (-v, +v)))

        jas = []
        jas.append(particles_field(s, "Current", "X", subplot(fig, gs, 2, 0), "$J_{\\phi}^"f"{s[0].lower()}""$", (-v, +v)))
        jas.append(particles_field(s, "Current", "Y", subplot(fig, gs, 2, 1), "$J_{\\phi}^"f"{s[0].lower()}""$", (-v, +v)))
        jas.append(particles_field(s, "Current", "Z", subplot(fig, gs, 2, 2), "$J_{\\phi}^"f"{s[0].lower()}""$", (-v, +v)))

        jzs = []
        jzs.append(particles_field(s, "Current", "X", subplot(fig, gs, 3, 0), "$J_z^"f"{s[0].lower()}""$", (-v, +v)))
        jzs.append(particles_field(s, "Current", "Y", subplot(fig, gs, 3, 1), "$J_z^"f"{s[0].lower()}""$", (-v, +v)))
        jzs.append(particles_field(s, "Current", "Z", subplot(fig, gs, 3, 2), "$J_z^"f"{s[0].lower()}""$", (-v, +v)))

        res_dir = f"{params_path}/Info_{s}"
        mkdir(res_dir)

        offset = 25
        t0 = 0
        t_range = create_t_range(t0, int(time / dts), offset)

        for t in t_range:
            plot_particles3(s, t)
