#!/usr/bin/env python3

from lib_common import *

def plot_fields3(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    def get_parsed_fields_xy(fr, fa, name):
        return (
            get_parsed_field(fr[0], name, "Y", "x", t),
            get_parsed_field(fa[0], name, "Y", "y", t),
            *get_parsed_field(fr[1], name, "Z", "", t))

    er[0].data, ea[0].data, er[1].data, ea[1].data = get_parsed_fields_xy(er, er, "E")
    br[0].data, ba[0].data, br[1].data, ba[1].data = get_parsed_fields_xy(br, br, "B")

    def get_parsed_fields_z(fz, name):
        return (
            get_parsed_field(fz[0], name, planes[0], "z", t),
            get_parsed_field(fz[1], name, planes[1], "z", t))

    ez[0].data, ez[1].data = get_parsed_fields_z(ez, "E")
    bz[0].data, bz[1].data = get_parsed_fields_z(bz, "B")

    bz[0].data = np.sqrt(np.square(bz[0].data) + np.square(br[0].data))

    for diag in er + ea + ez + br + ba + bz:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    annotate_x(er[0].axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(filename)

    for diag in er + ea + ez + br + ba + bz:
        diag.clear()


if __name__ == "__main__":
    ncols=3
    nrows=4

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.2, 1, 1.2, 1], figure=fig)

    er = []
    ea = []
    ez = []

    br = []
    ba = []
    bz = []

    planes = planes[1:]

    ea.append(electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\phi}$"))
    er.append(electric_field("Y", subplot(fig, gs, 1, 0), "$E_r$"))
    ez.append(electric_field("Y", subplot(fig, gs, 2, 0), "$E_z$"))

    ea.append(electric_field("Z", subplot(fig, gs, 0, 1), "$E_{\\phi}$"))
    er.append(electric_field("Z", subplot(fig, gs, 1, 1), "$E_r$"))
    ez.append(electric_field("Z", subplot(fig, gs, 2, 1), "$E_z$"))

    ba.append(electric_field("Y", subplot(fig, gs, 0, 2), "$B_{\\phi}$"))
    br.append(electric_field("Y", subplot(fig, gs, 1, 2), "$B_r$"))
    bz.append(magnetic_field("Y", subplot(fig, gs, 2, 2), "$|B|$"))

    ba.append(electric_field("Z", subplot(fig, gs, 0, 3), "$B_{\\phi}$"))
    br.append(electric_field("Z", subplot(fig, gs, 1, 3), "$B_r$"))
    bz.append(magnetic_field("Z", subplot(fig, gs, 2, 3), "$B_z$"))

    br[0].vmin_vmax = (-B0, +B0)

    bz[0].cmap = unsigned_cmap
    bz[0].vmin_vmax = (0, 2 * B0)

    res_dir = f"{params_path}/Fields"
    mkdir(res_dir)

    offset = 5
    t0 = 0 * offset # int(3000 / (dts * dt))
    t_range = create_t_range(t0, int(time / dts), offset)


    for t in t_range:
        plot_fields3(t)
