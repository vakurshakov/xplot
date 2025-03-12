#!/usr/bin/env python3

from plot import *

def plot_baligned_electric_potential(t):
    for tt in range(t, t + offset):
        update_data(tt)

    ep.data /= offset
    ep.axes_position.set_aspect(1)
    ep.draw(add_cbar=True)
    ep.draw_info()

    xc = data_shape["Y"][0] // 2
    xl = xc + 25
    w1 = -24
    w2 = +2
    xs1 = select_magnetic_line(xl + w1)
    xs2 = select_magnetic_line(xl + w2)
    zs = np.arange(0, data_shape["Y"][1])

    e_l = np.zeros(data_shape["Y"][1])
    for z in zs:
        rs = (np.arange(xs1[z], xs2[z] + 1) - xc) * dx
        area = cumulative_trapezoid(2 * np.pi * rs, rs, initial=0)[-1]
        e_l[z] += cumulative_trapezoid(ep.data[z, xs1[z]:(xs2[z] + 1)] * (2 * np.pi * rs), rs, initial=0)[-1] / (T_i * area)
        e_l[z] += cumulative_trapezoid(ep.data[z, (2 * xc - (xs2[z] + 1)):(2 * xc - xs1[z])] * (2 * np.pi * rs), rs, initial=0)[-1] / (T_i * area)
    e_l /= 2

    # e_l = ep.data[:,xc] / T_i

    ax = ep.axes_position
    ax.plot(+(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(+(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(-(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(-(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)

    ax = phi.axes_position
    ax.plot(zs * dz, e_l, label="$E_{\\|}(z) / T_i$", linewidth=3)
    ax.plot(zs * dz, -cumulative_trapezoid(e_l, zs * dx, initial=0), label="$\\varphi(z) / T_i$", linewidth=3)
    ax.legend(fontsize=0.8 * ticksize, loc="upper left")
    phi.draw_info()

    annotate_x(ep.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(f"{res_dir}/baligned_epotential_{str(t // offset).zfill(4)}.png")

    for diag in [ep, phi]:
        diag.clear()


def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep, "E", "Y", "x", t)
    _ez = get_parsed_field(ep, "E", "Y", "z", t)

    dot = (_er * br + _ez * bz)
    ep.data = agg(ep.data, np.divide(dot, b, where=(b > 1e-3), out=np.zeros_like(b)))


def select_magnetic_line(xl):
    xs = data_shape["Y"][0]
    zs = data_shape["Y"][1]
    xc = xs // 2
    zc = zs // 2

    # 2 * np.pi can be removed since comparison is relative
    b_f0 = np.sum(bz[zc,xc:xl] * np.arange(0, xl - xc))

    xmap = np.zeros(zs, dtype=np.int32)
    for z in np.arange(0, zs):
        b_fz = 0
        for x in np.arange(xc, xs):
            b_fz += bz[z, x] * (x - xc)
            xmap[z] = x

            if (b_fz >= b_f0):
                break
    return xmap


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

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    offset = 100
    t0 = 25 * offset
    plot_baligned_electric_potential(t0)
