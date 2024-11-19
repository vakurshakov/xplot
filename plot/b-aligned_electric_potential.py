#!/usr/bin/env python3

from lib_common import *

def plot_fields_b_aligned(t):
    for tt in range(t, t + offset):
        update_data(tt)

    ep.data /= offset
    ep.axes_position.set_aspect(1)
    ep.draw(add_cbar=True)
    ep.draw_info()

    zs = np.arange(0, data_shape["Y"][1])[:10] * dx
    xc = data_shape["Y"][0] // 2

    ax = phi.axes_position
    e_c = ep.data[:10,xc] / T_i
    ax.plot(zs, e_c, label="$E_{\\|}(z) / T_i$")
    ax.plot(zs, -cumulative_trapezoid(e_c, zs, dx, initial=0), label="$\\varphi(z) / T_i$")
    ax.legend(fontsize=0.8 * ssmol, loc="lower right")
    phi.draw_info()

    annotate_x(ep.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(f"{res_dir}/{str(t // offset).zfill(4)}.png")

    for diag in [ep, phi]:
        diag.clear()

def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep, "E", "Y", "x", t)
    _ez = get_parsed_field(ep, "E", "Y", "z", t)
    agg(ep, (_er * br + _ez * bz))


if __name__ == "__main__":
    ncols=2
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 2], height_ratios=[1] * nrows, figure=fig)

    ep = electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\|}$")
    ep.vmin_vmax = (-0.6e-3, +0.6e-3)

    phi = electric_field("Y")
    phi.axes_position = subplot(fig, gs, 1, 0)
    phi.set_axes_args(title="$\\varphi(z) = -\\int_0^z E_{\\|}(\\zeta) d\\zeta$")

    b = magnetic_field("Y")
    br = get_parsed_field(b, "B", "Y", "x", 0)
    bz = get_parsed_field(b, "B", "Y", "z", 0)
    b = np.hypot(br, bz)
    br = np.divide(br, b, where=(b > 1e-3), out=np.zeros_like(b))
    bz = np.divide(bz, b, where=(b > 1e-3), out=np.zeros_like(b))

    res_dir = f"{params_path}/EPotential_BAligned"
    mkdir(res_dir)

    offset = 100
    t0 = 25 * offset

    for t in create_t_range(t0, int(time / dts), offset):
        plot_fields_b_aligned(t)
        break
