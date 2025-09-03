#!/usr/bin/env python3

from plot import *

def plot_baligned_electric_potential(t):
    for tt in range(t, t + offset):
        update_data(tt)

    ep.data /= offset

    xc = data_shape["Y"][0] // 2
    zs = np.arange(0, data_shape["Y"][1])

    def calc_phi_r(r0, w1, w2):
        xl = xc + r0
        xs1 = select_magnetic_line(bz, xl + w1)
        xs2 = select_magnetic_line(bz, xl + w2)

        e_l1 = np.zeros(data_shape["Y"][1])
        e_l2 = np.zeros(data_shape["Y"][1])
        for z in zs:
            rs = (np.arange(xs1[z], xs2[z] + 1) - xc) * dx
            area = cumulative_trapezoid(2 * np.pi * rs, rs, initial=0)[-1]
            e_l1[z] += cumulative_trapezoid(ep.data[z, xs1[z]:(xs2[z] + 1)] * (2 * np.pi * rs), rs, initial=0)[-1] / (T_i * area)
            e_l2[z] += cumulative_trapezoid(ep.data[z, (2 * xc - (xs2[z] + 1)):(2 * xc - xs1[z])] * (2 * np.pi * rs), rs, initial=0)[-1] / (T_i * area)
        return e_l1, e_l2, -cumulative_trapezoid(e_l1, zs * dz, initial=0), -cumulative_trapezoid(e_l2, zs * dz, initial=0)

    w1 = -2
    w2 = +2
    xl_min = w2
    xl_max = int(ex / dx - 2)

    phi_r.data = np.zeros((len(zs), 2 * xl_max))
    _, _, d1, d2 = calc_phi_r(xl_min, w1, w2)
    
    for i in range(w1, w2):
        phi_r.data[:, xl_max+i] = (d1 + d2) / 2

    for i in np.arange(w2, xl_max):
        _, _, d1, d2 = calc_phi_r(i, w1, w2)
        phi_r.data[:, xl_max-i] = d1
        phi_r.data[:, xl_max+i] = d2
    phi_r.data[:, 0] = d1

    def draw_phi_r(phi_r, data):
        ax = phi_r.axes_position
        ax.plot(zs * dz, (data[0] + data[1]) / 2, linewidth=3)
        ax.plot(zs * dz, data[2], linewidth=3)
        ax.plot(zs * dz, data[3], linewidth=3)
        # ax.legend(fontsize=0.8 * ssmol, loc="lower left")
        ax.grid(alpha=0.6)
        phi_r.draw_info()
    draw_phi_r(phi_1, calc_phi_r(xl_min, w1, w2))
    draw_phi_r(phi_2, calc_phi_r(xl_max, w1, w2))

    ep.axes_position.set_aspect(1)
    ep.draw(add_cbar=True)
    ep.draw_info()

    phi_r.axes_position.set_aspect(1)
    phi_r.draw(add_cbar=True)
    phi_r.draw_info()

    def draw_phi_info(xl):
        xs1 = select_magnetic_line(bz, xc + xl + w1)
        xs2 = select_magnetic_line(bz, xc + xl + w2)
        zmin = (zs * dz)[0]
        zmax = (zs * dz)[-1]
        zl = [zmin, zmax]
        ax = ep.axes_position
        ax.plot(-(xs1 - xc) * dx - 0.4, zs * dz, color="black", linewidth=2)
        ax.plot(+(xs1 - xc) * dx + 0.4, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs2 - xc) * dx - 0.2, zs * dz, color="black", linewidth=2)
        ax.plot(+(xs2 - xc) * dx + 0.2, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs1 - xc) * dx - 0.2, zs * dz, color="C1", linewidth=2)
        ax.plot(+(xs1 - xc) * dx + 0.2, zs * dz, color="C2", linewidth=2)
        ax.plot(-(xs2 - xc) * dx      , zs * dz, color="C1", linewidth=2)
        ax.plot(+(xs2 - xc) * dx      , zs * dz, color="C2", linewidth=2)
        ax = phi_r.axes_position
        ax.plot([-(xl + w1) * dx + 0.15]*2, zl, color="black", linewidth=2)
        ax.plot([+(xl + w1) * dx - 0.15]*2, zl, color="black", linewidth=2)
        ax.plot([-(xl + w2) * dx + 0.15]*2, zl, color="black", linewidth=2)
        ax.plot([+(xl + w2) * dx - 0.16]*2, zl, color="black", linewidth=2)
        ax.plot([-(xl + w1) * dx + 0.10]*2, zl, color="C1", linewidth=2)
        ax.plot([+(xl + w1) * dx - 0.10]*2, zl, color="C2", linewidth=2)
        ax.plot([-(xl + w2) * dx + 0.10]*2, zl, color="C1", linewidth=2)
        ax.plot([+(xl + w2) * dx - 0.10]*2, zl, color="C2", linewidth=2)
    draw_phi_info(xl_max)
    draw_phi_info(xl_min)

    fig.suptitle(f"$t / \\tau = {{{t * dts / tau:.2f}}}$", size=big, bbox=bbox, y=0.985)

    fig.tight_layout()

    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    fig.savefig(filename)
    print("---------- Processed", filename)

    for diag in [ep, phi_r, phi_1, phi_2]:
        diag.clear()


def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep, "E", "Y", "x", t)
    _ez = get_parsed_field(ep, "E", "Y", "z", t)

    dot = (_er * br + _ez * bz)
    ep.data = agg(ep.data, np.divide(dot, b, where=(b > 1e-3), out=np.zeros_like(b)))


if __name__ == "__main__":
    ncols=4
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1, 1.5, 1.5], height_ratios=[1] * nrows, figure=fig)

    ep = electric_field("Y", subplot(fig, gs, 0, 0), "")
    ep.vmin_vmax = (-2e-3, +2e-3)
    ep.axes_args["title"] = "$E_{\\|}(x, z)$"

    phi_r = Field(None, subplot(fig, gs, 1, 0), cmap=signed_cmap)
    phi_r.vmin_vmax=(-3, 3)
    bx = -28 * dx
    ex = +28 * dx
    by = boundaries["Y"][2]
    ey = boundaries["Y"][3]
    phi_r.boundaries = (bx, ex, by, ey)
    phi_r.axes_args["title"] = "$\\varphi(\\xi, z)$"
    phi_r.axes_args["xlim"] = (bx, ex)
    phi_r.axes_args["ylim"] = (by, ey)
    phi_r.axes_args["xlabel"] = "$\\xi,~c/\\omega_{pe}$"
    phi_r.axes_args["ylabel"] = "$z,~c/\\omega_{pe}$"
    phi_r.axes_args["xticks"] = np.linspace(bx, ex, 5)
    phi_r.axes_args["yticks"] = np.linspace(by, ey, 5)

    phi_1 = electric_field("Y")
    phi_1.axes_position = subplot(fig, gs, 2, 0)
    phi_1.axes_args["title"] = "$\\varphi(z) = -\\int_0^z E_{\\|}(\\zeta) d\\zeta,~\\xi = \\rm 2$"
    phi_1.axes_args["xlabel"] = "$z,~c/\\omega_{pe}$"
    phi_1.axes_args["ylim"] = (-1, 5)

    phi_2 = electric_field("Y")
    phi_2.axes_position = subplot(fig, gs, 3, 0)
    phi_2.axes_args["title"] = f"$\\varphi(z),~\\xi = \\pm {ex-2:.0f}$"
    phi_2.axes_args["xlabel"] = "$z,~c/\\omega_{pe}$"
    phi_2.axes_args["ylim"] = (-1, 5)

    b = magnetic_field("Y")
    br = get_parsed_field(b, "B", "Y", "x", 0)
    bz = get_parsed_field(b, "B", "Y", "z", 0)
    b = np.hypot(br, bz)

    res_dir = f"{params_path}/B-Aligned_electric_potential"
    mkdir(res_dir)

    offset = 10
    t_range = create_t_range(offset, int(time / dts) - offset, offset)

    for t in t_range[::-1]:
        plot_baligned_electric_potential(t)

