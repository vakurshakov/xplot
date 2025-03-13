#!/usr/bin/env python3

from plot import *

def plot_fields3(t: int):
    filename = f"{res_dir}/{get_formatted_time(t // offset)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    # TODO: Should be moved into `xplot.plot.Field` or abstracted in some form (cyl, dec problem)
    er[0].data = get_parsed_field(er[0].path_to_file, "x", t)
    ea[0].data = get_parsed_field(ea[0].path_to_file, "y", t)
    ez[0].data = get_parsed_field(ez[0].path_to_file, "z", t)

    er[1].data = get_parsed_field(er[1].path_to_file, "x", t)
    ea[1].data = get_parsed_field(ea[1].path_to_file, "y", t)
    ez[1].data = get_parsed_field(ez[1].path_to_file, "z", t)

    br[0].data = get_parsed_field(br[0].path_to_file, "x", t)
    ba[0].data = get_parsed_field(ba[0].path_to_file, "y", t)
    bz[0].data = get_parsed_field(bz[0].path_to_file, "z", t)

    br[1].data = get_parsed_field(br[1].path_to_file, "x", t)
    ba[1].data = get_parsed_field(ba[1].path_to_file, "y", t)
    bz[1].data = get_parsed_field(bz[1].path_to_file, "z", t)

    for diag in er + ea + ez + br + ba + bz:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    # annotate_x(er[0].axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)
    annotate_x(er[0].axes_position, "$t = {" f"{t * dts:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(filename)

    for diag in er + ea + ez + br + ba + bz:
        diag.clear()


if __name__ == "__main__":
    ncols=3
    nrows=4

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1, 1, 1, 1], figure=fig)

    er = []
    ea = []
    ez = []

    br = []
    ba = []
    bz = []

    ea.append(electric_field("Y", subplot(fig, gs, 0, 0), "$E_x$"))
    er.append(electric_field("Y", subplot(fig, gs, 1, 0), "$E_y$"))
    ez.append(electric_field("Y", subplot(fig, gs, 2, 0), "$E_z$"))

    ea.append(electric_field("Z", subplot(fig, gs, 0, 1), "$E_x$"))
    er.append(electric_field("Z", subplot(fig, gs, 1, 1), "$E_y$"))
    ez.append(electric_field("Z", subplot(fig, gs, 2, 1), "$E_z$"))

    ba.append(magnetic_field("Y", signed_cmap, (  -2e-2,   +2e-2), subplot(fig, gs, 0, 2), "$B_x$"))
    br.append(magnetic_field("Y", signed_cmap, (  -2e-2,   +2e-2), subplot(fig, gs, 1, 2), "$B_y$"))
    bz.append(magnetic_field("Y", signed_cmap, (B0-2e-2, B0+2e-2), subplot(fig, gs, 2, 2), "$B_z$"))

    ba.append(magnetic_field("Z", signed_cmap, (  -2e-2,   +2e-2), subplot(fig, gs, 0, 3), "$B_x$"))
    br.append(magnetic_field("Z", signed_cmap, (  -2e-2,   +2e-2), subplot(fig, gs, 1, 3), "$B_y$"))
    bz.append(magnetic_field("Z", signed_cmap, (B0-2e-2, B0+2e-2), subplot(fig, gs, 2, 3), "$B_z$"))

    res_dir = f"{params_path}/Fields"
    mkdir(res_dir)

    offset = int(dts / dt)
    t_range = create_t_range(0, Nt, offset)

    for t in t_range:
        plot_fields3(t)
