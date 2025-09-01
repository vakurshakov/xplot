#!/usr/bin/env python3

from plot import *

def plot_fields3(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    def get_parsed_fields_xy(fr, fa, name):
        return (
            get_parsed_field(fr[0], name, "Y", "x", t),
            get_parsed_field(fa[0], name, "Y", "y", t),
            *get_parsed_field(fr[1], name, "Z", "", t))

    er[0].data, ea[0].data, er[1].data, ea[1].data = get_parsed_fields_xy(er, er, "E")

    br[0].data = parse_file(get_parsed_file(t, br[0].path_to_file), 0) - BORY
    ba[0].data = parse_file(get_parsed_file(t, ba[0].path_to_file), 1) - B0AY
    br[1].data = parse_file(get_parsed_file(t, br[1].path_to_file), 1) - BORX
    ba[1].data = parse_file(get_parsed_file(t, ba[1].path_to_file), 0) - B0AX

    def get_parsed_fields_z(fz, name):
        return (
            get_parsed_field(fz[0], name, planes[0], "z", t),
            get_parsed_field(fz[1], name, planes[1], "z", t))

    ez[0].data, ez[1].data = get_parsed_fields_z(ez, "E")
    bz[0].data, bz[1].data = get_parsed_fields_z(bz, "B")

    bz[0].data = np.sqrt(np.square(bz[0].data) + np.square(br[0].data + BORY))

    for diag in er + ea + ez + br + ba + bz:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    br[0].title.set_bbox(bbox)
    ba[0].title.set_bbox(bbox)
    br[1].title.set_bbox(bbox)
    ba[1].title.set_bbox(bbox)
    bz[0].title.set_bbox(bbox)

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=0.99, bbox=bbox, fontsize=big)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(filename)

    for diag in er + ea + ez + br + ba + bz:
        diag.clear()

B0AX = magnetic_field("X")
B0AX = parse_file(get_parsed_file(0, B0AX.path_to_file), 0)

B0AY = magnetic_field("Y")
B0AY = parse_file(get_parsed_file(0, B0AY.path_to_file), 1)

BORX = magnetic_field("X")
BORX = parse_file(get_parsed_file(0, BORX.path_to_file), 1)

BORY = magnetic_field("Y")
BORY = parse_file(get_parsed_file(0, BORY.path_to_file), 0)

if __name__ == "__main__":
    ncols=4
    nrows=3

    fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    er = []
    ea = []
    ez = []

    br = []
    ba = []
    bz = []

    planes = planes[1:]

    er.append(electric_field("Y", subplot(fig, gs, 0, 0), "$E_r$"))
    ea.append(electric_field("Y", subplot(fig, gs, 0, 1), "$E_{\\phi}$"))
    ez.append(electric_field("Y", subplot(fig, gs, 0, 2), "$E_z$"))

    er.append(electric_field("Z", subplot(fig, gs, 1, 0), "$E_r$"))
    ea.append(electric_field("Z", subplot(fig, gs, 1, 1), "$E_{\\phi}$"))
    ez.append(electric_field("Z", subplot(fig, gs, 1, 2), "$E_z$"))

    vmap = (-2e-2, +2e-2)
    br.append(magnetic_field("Y", subplot(fig, gs, 2, 0), "$\\delta B_x$", vmap))
    ba.append(magnetic_field("Y", subplot(fig, gs, 2, 1), "$\\delta B_y$", vmap))
    bz.append(magnetic_field("Y", subplot(fig, gs, 2, 2), "$|B|$", (-B0, +B0)))

    br.append(magnetic_field("X", subplot(fig, gs, 3, 0), "$\\delta B_y$", vmap))
    ba.append(magnetic_field("X", subplot(fig, gs, 3, 1), "$\\delta B_x$", vmap))
    bz.append(magnetic_field("Z", subplot(fig, gs, 3, 2), "$B_z$", (0, B0), unsigned_cmap))

    bz[0].cmap = unsigned_cmap
    bz[0].vmin_vmax = (0, 2 * B0)

    res_dir = f"{params_path}/Fields"
    mkdir(res_dir)

    offset = 25
    t0 = 0
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range[::-1]:
        plot_fields3(t)
