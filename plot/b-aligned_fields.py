#!/usr/bin/env python3

from lib_common import *

def plot_baligned_fields(t):
    for tt in range(t, t + offset):
        update_data(tt)

    for diag in [ep]: # [ep, etr, etz]:
        diag.data /= offset
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    annotate_x(ep.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(f"{res_dir}/baligned_field_{str(t // offset).zfill(4)}.png")

    for diag in [ep]: # [ep, etr, etz]:
        diag.clear()

def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep, "E", "Y", "x", t)
    _ez = get_parsed_field(ep, "E", "Y", "z", t)

    dot = (_er * br + _ez * bz)

    agg(ep, dot)

    # _er -= ep.data * br
    # _ez -= ep.data * bz
    # agg(etr, _er)
    # agg(etz, _ez)

if __name__ == "__main__":
    ncols=1
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    ep = electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\|}$")
    # etr = electric_field("Y", subplot(fig, gs, 1, 0), "$E_{\\perp, r}$")
    # etz = electric_field("Y", subplot(fig, gs, 2, 0), "$E_{\\perp, z}$")

    for diag in [ep]:
        diag.vmin_vmax = (-0.6e-3, +0.6e-3)

    b = magnetic_field("Y")
    br = get_parsed_field(b, "B", "Y", "x", 0)
    bz = get_parsed_field(b, "B", "Y", "z", 0)
    b = np.hypot(br, bz)
    br = np.divide(br, b, where=(b > 1e-3), out=np.zeros_like(b))
    bz = np.divide(bz, b, where=(b > 1e-3), out=np.zeros_like(b))

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    offset = 100
    t0 = 25 * offset
    plot_baligned_fields(t0)
