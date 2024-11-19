#!/usr/bin/env python3

from lib_common import *

def plot_baligned_info(t):
    for tt in range(t, t + offset):
        update_data(tt)

    for diag in [ep, n, vp, vt, prr, pzz]:
        diag.data /= offset

    xc = data_shape["Y"][0] // 2
    xl = xc + 25
    w1 = -2 # 4
    w2 = +2
    xs1 = select_magnetic_line(xl + w1)
    xs2 = select_magnetic_line(xl + w2)
    zs = np.arange(0, data_shape["Y"][1])

    for diag in [ep, b]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

        ax = diag.axes_position
        ax.plot(+(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(+(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)

    zs = np.arange(0, data_shape["Y"][1])
    e_l = average_over_magnetic_tube(ep.data, zs, xs1, xs2, xc) / T_i
    n_l = average_over_magnetic_tube(n.data, zs, xs1, xs2, xc)
    vp_l = average_over_magnetic_tube(vp.data, zs, xs1, xs2, xc) / np.sqrt(T_i / mi_me)
    vt_l = average_over_magnetic_tube(vt.data, zs, xs1, xs2, xc) / np.sqrt(T_i / mi_me)
    prr_l = average_over_magnetic_tube(prr.data, zs, xs1, xs2, xc) / T_i
    pzz_l = average_over_magnetic_tube(pzz.data, zs, xs1, xs2, xc) / T_i

    n_zero_tolerance = 1e-3
    vp_l = np.divide(vp_l, n_l, where=(n_l > n_zero_tolerance), out=np.zeros_like(vp_l))
    vt_l = np.divide(vt_l, n_l, where=(n_l > n_zero_tolerance), out=np.zeros_like(vt_l))

    prr_l = np.divide(prr_l, n_l, where=(n_l > n_zero_tolerance), out=np.zeros_like(prr_l))
    pzz_l = np.divide(pzz_l, n_l, where=(n_l > n_zero_tolerance), out=np.zeros_like(pzz_l))

    ax = phi.axes_position
    ax.plot(zs * dz, e_l, label="$E_{\\|}(z) / T_i$", linewidth=3)
    ax.plot(zs * dz, -cumulative_trapezoid(e_l, zs * dx, initial=0), label="$\\varphi(z) / T_i$", linewidth=3)
    ax.legend(fontsize=0.8 * ssmol, loc="upper left")

    n.axes_position.plot(zs * dz, n_l, linewidth=2)
    vp.axes_position.plot(zs * dz, vp_l, linewidth=2)
    vt.axes_position.plot(zs * dz, vt_l, linewidth=2)
    prr.axes_position.plot(zs * dz, prr_l, linewidth=2)
    pzz.axes_position.plot(zs * dz, pzz_l, linewidth=2)

    for diag in [phi, n, vp, vt, prr, pzz]:
        diag.draw_info()

    # phi.axes_position.set_xlim(0, 50)
    print(np.max(-cumulative_trapezoid(e_l, zs * dx, initial=0)))

    annotate_x(ep.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()
    fig.savefig(f"{res_dir}/baligned_info_{str(t // offset).zfill(4)}.png")

    for diag in [ep, phi]:
        diag.clear()


def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    def align(fr, fz):
        dot = (fr * br + fz * bz)
        return np.divide(dot, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))

    _er = get_parsed_field(ep, "E", "Y", "x", t)
    _ez = get_parsed_field(ep, "E", "Y", "z", t)
    ep.data = agg(ep.data, align(_er, _ez))

    n.data = agg(n.data, get_parsed_scalar(n, t))
    prr.data = agg(prr.data, get_parsed_scalar(prr, t))
    pzz.data = agg(pzz.data, get_parsed_scalar(pzz, t))

    _vr = get_parsed_field(vp, "E", "Y", "x", t)
    _vz = get_parsed_field(vt, "E", "Y", "z", t)
    _vp = align(_vr, _vz)
    vp.data = agg(vp.data, _vp)

    _vtr = _vr - _vp * np.divide(br, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    _vtz = _vz - _vp * np.divide(bz, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    vt.data = agg(vt.data, np.hypot(_vtr, _vtz))


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

def average_over_magnetic_tube(data, zs, xs1, xs2, xc):
    f_l = np.zeros(data_shape["Y"][1])
    for z in zs:
        rs = (np.arange(xs1[z], xs2[z] + 1) - xc) * dx
        area = cumulative_trapezoid(2 * np.pi * rs, rs, initial=0)[-1]
        f_l[z] += cumulative_trapezoid(data[z, xs1[z]:(xs2[z] + 1)] * (2 * np.pi * rs), rs, initial=0)[-1] / area
        f_l[z] += cumulative_trapezoid(data[z, (2 * xc - (xs2[z] + 1)):(2 * xc - xs1[z])] * (2 * np.pi * rs), rs, initial=0)[-1] / area
    return f_l / 2


if __name__ == "__main__":
    ncols=3
    nrows=3

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1.5, 1.5], height_ratios=[1] * nrows, figure=fig)

    ep = electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\|}$")
    ep.vmin_vmax = (-0.6e-3, +0.6e-3)

    phi = electric_field("Y", subplot(fig, gs, 1, 0))
    phi.set_axes_args(title="$\\varphi(z) = -\\int_0^z E_{\\|}(\\zeta) d\\zeta$")

    n = particles_field("Ions", "Dens", "Y", subplot(fig, gs, 2, 0))
    vp = particles_field("Ions", "Current", "Y", subplot(fig, gs, 1, 1))
    vt = particles_field("Ions", "Current", "Y", subplot(fig, gs, 2, 1))
    prr = particles_field("Ions", "Pxx", "Y", subplot(fig, gs, 2, 2))
    pzz = particles_field("Ions", "Pzz", "Y", subplot(fig, gs, 1, 2))

    n.set_axes_args(title="$n_i / n_0$") #, (0, 1.5))
    vp.set_axes_args(title="$v_{\\|}^i / v_T$") #, (-0.02, +0.02))
    vt.set_axes_args(title="$v_{\\perp}^i / v_T$") #, (-0.02, +0.02))
    prr.set_axes_args(title="$\\Pi_{rr}^i / (n_i T_i)$") #, (0, +0.04))
    pzz.set_axes_args(title="$\\Pi_{zz}^i / (n_i T_i)$") #, (0, +0.04))

    b = magnetic_field("Y", subplot(fig, gs, 0, 1), "$|B|$")
    br = get_parsed_field(b, "B", "Y", "x", 0)
    bz = get_parsed_field(b, "B", "Y", "z", 0)
    b.data = np.hypot(br, bz)

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    offset = 100
    t0 = int(25 * offset) # * (2 / 3))
    plot_baligned_info(t0)
