#!/usr/bin/env python3

from plot import *

def plot_currents(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    ja_xi.data = get_parsed_field(ja_xi, "E", "X", "x", t)
    ja_xe.data = get_parsed_field(ja_xe, "E", "X", "x", t)

    ja_zi.data = get_parsed_field(ja_zi, "E", "Z", "", t)[1]
    ja_ze.data = get_parsed_field(ja_ze, "E", "Z", "", t)[1]

    for diag in [ ja_xi, ja_xe, ja_zi, ja_ze ]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    fig.suptitle(f"$t / \\tau = {t * dts / tau:.3f}$", x=0.5, bbox=bbox, fontsize=labelsize)
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

    for diag in [ ja_xi, ja_xe, ja_zi, ja_ze ]:
        diag.clear()


if __name__ == "__main__":
    ncols=2
    nrows=2

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.2, 1], figure=fig)

    v = 0.02
    ja_xi = particles_field("Ions", "Current", "X", subplot(fig, gs, 0, 0), "$J_{\\phi}^i$", (-v, +v))
    ja_zi = particles_field("Ions", "Current", "Z", subplot(fig, gs, 0, 1), "$J_{\\phi}^i$", (-v, +v))

    ja_xe = particles_field("Electrons", "Current", "X", subplot(fig, gs, 1, 0), "$J_{\\phi}^e$", (-v, +v))
    ja_ze = particles_field("Electrons", "Current", "Z", subplot(fig, gs, 1, 1), "$J_{\\phi}^e$", (-v, +v))

    res_dir = f"{params_path}/Currents"
    mkdir(res_dir)

    offset = 5
    t0 = 0
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range:
        plot_currents(t)
