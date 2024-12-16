#!/usr/bin/env python3

from final import *

def plot_baligned_info(t):
    for tt in range(t, t + offset):
        update_data(tt)

    for diag in [ep, jai, jae]:
        diag.data /= offset

    xc = data_shape["Y"][0] // 2
    zc = data_shape["Y"][1] // 2
    xl = xc + 1
    w1 = -1 # 4
    w2 = +1
    xs1 = select_magnetic_line(bz, xl + w1)
    xs2 = select_magnetic_line(bz, xl + w2)
    zs = np.arange(0, data_shape["Y"][1])

    for diag in [ep, jai, jae]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True, cbar_exponential=True)
        diag.draw_info()

    ax = ep.axes_position
    ax.plot(+(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(+(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(-(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
    ax.plot(-(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)

    annotate_x(jai.axes_position, "$t / \\tau = {" f"{t * dts / tau:.1f}" "}$", y=1.2)
    print(f"{int(t * dts / dt)} [dt] {int(t)} [dts] {t * dts / tau:.3f} [tau]")

    fig.tight_layout()

    filename = f"{res_dir}/fields_{int(t * dts / dt)}dt_{t * dts / tau:.1f}tau"
    fig.savefig(filename + ".pdf")

    for diag in [ep, jai, jae]:
        diag.clear()


def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep.path_to_file, "E", "Y", "x", t)
    _ez = get_parsed_field(ep.path_to_file, "E", "Y", "z", t)
    ep.data = agg(ep.data, align(_er, _ez, br, bz, b.data))

    jai.data = agg(jai.data, get_parsed_field(jai.path_to_file, "E", "Z", "", t)[1])
    jae.data = agg(jae.data, get_parsed_field(jae.path_to_file, "E", "Z", "", t)[1])


if __name__ == "__main__":
    ncols=3
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1.2, 1.2], height_ratios=[1] * nrows, figure=fig)

    v = 6e-4
    ep = electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\|}$")
    ep.vmin_vmax = (-v, v)

    v = 6e-3
    jai = particles_field("Ions", "Current", "Z", subplot(fig, gs, 1, 0), "$J_{\\phi}^i$", (-v, +v))
    jae = particles_field("Electrons", "Current", "Z", subplot(fig, gs, 2, 0), "$J_{\\phi}^e$", (-v, +v))

    b = magnetic_field("Y")
    br = get_parsed_field(b.path_to_file, "B", "Y", "x", 0)
    bz = get_parsed_field(b.path_to_file, "B", "Y", "z", 0)
    b.data = np.hypot(br, bz)

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    offset = 100
    ts = [
        int(2.0 * tau / dts),
        int(3.0 * tau / dts),
        int(3.5 * tau / dts),
        int(4.0 * tau / dts),
    ]

    for t in ts:
        plot_baligned_info(t)
