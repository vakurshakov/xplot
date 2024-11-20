#!/usr/bin/env python3

from plot import *

def plot_electric_fields(t):
    for tt in range(t, t + offset):
        update_data(tt)

    for diag in er + ea + ez:
        diag.data /= offset
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    annotate_x(er[0].axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(f"{res_dir}/{str(t // offset).zfill(4)}.png")

    for diag in er + ea + ez:
        diag.clear()

def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    er[0].data = agg(er[0].data, get_parsed_field(er[0], "E", "X", "y", t))
    ea[0].data = agg(ea[0].data, get_parsed_field(ea[0], "E", "X", "x", t))
    ez[0].data = agg(ez[0].data, get_parsed_field(ez[0], "E", "X", "z", t))
    er[1].data = agg(er[1].data, get_parsed_field(er[1], "E", "Y", "x", t))
    ea[1].data = agg(ea[1].data, get_parsed_field(ea[1], "E", "Y", "y", t))
    ez[1].data = agg(ez[1].data, get_parsed_field(ez[1], "E", "Y", "z", t))


if __name__ == "__main__":
    ncols=3
    nrows=2

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.2, 1.2], figure=fig)

    er = []
    ea = []
    ez = []

    ea.append(electric_field("X", subplot(fig, gs, 0, 0), "$E_{\\phi}$"))
    er.append(electric_field("X", subplot(fig, gs, 1, 0), "$E_r$"))
    ez.append(electric_field("X", subplot(fig, gs, 2, 0), "$E_z$"))

    ea.append(electric_field("Y", subplot(fig, gs, 0, 1), "$E_{\\phi}$"))
    er.append(electric_field("Y", subplot(fig, gs, 1, 1), "$E_r$"))
    ez.append(electric_field("Y", subplot(fig, gs, 2, 1), "$E_z$"))

    for diag in er + ea + ez:
        diag.vmin_vmax = (-0.6e-2, +0.6e-2)

    res_dir = f"{params_path}/Fields_Electric"
    mkdir(res_dir)

    offset = 100
    t0 = 500 * 5

    for t in  [t0]: # create_t_range(t0, int(time / dts), offset):
        plot_electric_fields(t)
