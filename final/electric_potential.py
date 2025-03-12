#!/usr/bin/env python3

from lib_common import *
from plot_fields import *

from scipy.integrate import cumulative_trapezoid

def plot_electric_potential():
    rs = (np.arange(0, data_shape["Z"][0] // 2) + 0.5) * dx

    def time_average(func, t0, dt):
        result = np.zeros_like(func(t0))
        for t in range(t0 - dt, t0 + 1):
            result += func(t)
        return result / dt

    ts = [ 1, 3, 4.47 ]
    for i, t in enumerate(ts):
        er.axes_position = subplot(fig, gs, i, 0)

        tt = int(t * tau / dts)
        er.data = get_parsed_field(er, "E", "X", "y", tt)
        er.draw(add_cbar=True)
        er.cbar = None
        er.draw_info()
        ax = er.axes_position

        annotate_x(ax, f"$t = {t:.1f}\,\\tau$", y=0.9, x=0.19, size=labelsize)

        e_avg = phi_averaged(time_average(lambda t0: get_parsed_field(electric_field("Z"), "E", "Z", "", t0), tt, 2)[0], R_MAP)
        phi.axes_position.plot(rs, -cumulative_trapezoid(e_avg, rs, dx, initial=0) / T_i, label=f"$t = {t:.1f}\,\\tau$", linewidth=2)

    v = (0, 4)
    ax = phi.axes_position
    phi.draw_info(xlim=(0, 30), ylim=v, yticks=np.linspace(*v, 5))
    ax.set_xlabel("$r,~c/\\omega_{pe}$", fontsize=labelsize)
    ax.legend(fontsize=ticksize, loc="upper left")
    ax.grid()

    fig.tight_layout()

    filename = f"{res_dir}/figure001"
    fig.savefig(filename + ".png")
    fig.savefig(filename + ".pdf")

if __name__ == "__main__":
    ncols=4
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    er = electric_field("X", None, "")
    er.set_axes_args(title="$E_r(y,\,z)$")

    phi = electric_field("X")

    phi.axes_position = subplot(fig, gs, 3, 0)
    phi.set_axes_args(title="$e \\varphi(r) / T_i$")

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    plot_electric_potential()
