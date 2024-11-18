#!/usr/bin/env python3

from lib_common import *
from plot_fields import *

from scipy.integrate import cumulative_trapezoid

def plot_electric_potential(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    def time_average(func, t0, dt):
        result = np.zeros_like(func(t0))
        for t in range(t0 - dt, t0 + 1):
            result += func(t)
        return result / dt

    er.data = phi_averaged(time_average(lambda t0: get_parsed_field(er, "E", "Z", "", t0), t, 10)[0], R_MAP)

    def avg1(d): return sliding_average(d, 3)
    def avg2(d): return sliding_average(d, 7)

    rs = (np.arange(0, data_shape["Z"][0] // 2) + 0.5) * dx
    avg1_rs = avg1(rs)
    avg2_rs = avg2(rs)

    ax = er.axes_position
    ax.plot(rs, er.data / T_i, label="$E_r(r) / T_i$")
    ax.plot(rs, -cumulative_trapezoid(er.data, rs, dx, initial=0) / T_i, label="$\\varphi(r) / T_i$")

    annotate_x(ax, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2, size=smol)

    v = (-2, 10)
    er.draw_info(ylim=v, yticks=np.linspace(*v, 7))
    ax.set_xlabel("$r,~c/\\omega_{pe}$", fontsize=ssmol)
    ax.legend(fontsize=ssmol * 0.8, loc="upper left")
    ax.grid()

    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

    for diag in [ er ]:
        diag.clear()

if __name__ == "__main__":
    ncols=1
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.5, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    er = electric_field("Z")
    er.axes_position = subplot(fig, gs, 0, 0)
    er.set_axes_args(title="Electric potential, $\\varphi(r) = -\\int_0^r E_r(r) dr$")

    res_dir = f"{params_path}/Electric_potential"
    mkdir(res_dir)

    offset = 5

    t0 = 10 # 800 * 5
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range:
        plot_electric_potential(t)
