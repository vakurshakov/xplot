#!/usr/bin/env python

from lib_common import *
from plot_fields import *
from plot_particles import *

if __name__ == "__main__":
    ncols=2
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1.2], height_ratios=[1] * nrows, figure=fig)

    t0 = int(time / dts)

    ni_y = particles_field("Ions", "Dens", "Y", subplot(fig, gs, 0, 0), "$n$", (0,2), unsigned_cmap)
    ne_y = particles_field("Electrons", "Dens", "Y", None, "", (0,0))
    ni_y.data = \
        parse_file(f"{get_prefix(t0)}/{ni_y.path_to_file}_{str(t0).zfill(4)}") - \
        parse_file(f"{get_prefix(t0)}/{ne_y.path_to_file}_{str(t0).zfill(4)}")

    # ni_z = particles_field("Ions", "Dens", "Z", subplot(fig, gs, 1, 0), "$n$", (0,2), unsigned_cmap)
    # ne_z = particles_field("Electrons", "Dens", "Z", None, "", (0,0))
    # ni_z.data = \
    #     parse_file(f"{get_prefix(t0)}/{ni_z.path_to_file}_{str(t0).zfill(4)}") - \
    #     parse_file(f"{get_prefix(t0)}/{ne_z.path_to_file}_{str(t0).zfill(4)}")

    bz = magnetic_field("Z", subplot(fig, gs, 1, 0), "$B_z$")
    bz.cmap = unsigned_cmap
    bz.vmin_vmax = (0, 0.2)
    bz.data = get_parsed_field(bz, "B", planes[1], "z", t0)
    bz.data = np.abs(bz.data)

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    for d in [ ni_y, bz ]:
        d.axes_position.set_aspect(1)
        d.draw(add_cbar=True)
        d.draw_info()

    # annotate_x(ni.axes_position, "$t / \\tau = {" f"{t0 * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()

    filename = f"{res_dir}/title.png"
    print(filename, t0, "[dts]", f"{t0 * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

