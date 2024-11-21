#!/opt/anaconda3/bin/python3
##/usr/bin/python3

from lib_common import *
from plot_fields import *
from plot_particles import *

B0 = 0.191

def plot_spatial_beta():
    ts = [ 1, 3, 4.47 ]
    for i, t in enumerate(ts):
        tt = int(t * tau / dts)

        bz.axes_position = subplot(fig, gs, i, 0)
        beta.axes_position = subplot(fig, gs, i, 1)

        bz.data = get_parsed_field(bz, "B", "Z", "z", tt)
        prr_i.data = get_parsed_scalar(prr_i, tt)
        prr_e.data = get_parsed_scalar(prr_e, tt)

        beta.data = 2 * (prr_i.data + prr_e.data) / (B0 * B0)

        bz.draw(add_cbar=True)
        bz.cbar = None
        bz.draw_info()
        bz.axes_position.set_aspect(1)
        annotate_x(bz.axes_position, f"$t = {t:.1f}\,\\tau$", y=0.9, x=0.18, size=smol)

        beta.draw(add_cbar=True)
        beta.cbar = None
        beta.draw_info()
        beta.axes_position.set_aspect(1)
        annotate_x(beta.axes_position, f"$t = {t:.1f}\,\\tau$", y=0.9, x=0.18, size=smol)

    ts = np.arange(0, int(time / dts) + 1, 1) * dts / tau
    b_min = np.load(f"{res_dir}/b_min_t.npy", allow_pickle=True)
    beta_max = np.load(f"{res_dir}/beta_max_t.npy", allow_pickle=True)

    b_t = Field(None, subplot(fig, gs, 3, 0))
    b_t.set_axes_args(
        title="\\rm Magnetic field",
        ylabel="$B_{z,\,\\min}$",
        xlabel="$t / \\tau$",
    )
    ax = b_t.axes_position
    ax.plot(ts, b_min, linewidth=2)
    ax.grid()
    b_t.draw_info(xlim=(0, 5), ylim=(0, 0.2))

    beta_t = Field(None, subplot(fig, gs, 3, 1))
    beta_t.set_axes_args(
        title="\\rm Plasma beta",
        ylabel="$\\beta_{\\max}$",
        xlabel="$t / \\tau$",
    )
    ax = beta_t.axes_position
    ax.plot(sliding_average(ts), sliding_average(beta_max), linewidth=2)
    ax.grid()
    beta_t.draw_info(xlim=(0, 5), ylim=(0, 1.5))

    fig.tight_layout()
    # plt.show()

    filename = f"{res_dir}/figure001"
    fig.savefig(filename + ".png")
    fig.savefig(filename + ".pdf")


if __name__ == "__main__":
    ncols=4
    nrows=2

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    bz = magnetic_field("Z", None, "")
    bz.set_axes_args(title="$B_z(x,\,y)$")

    beta = Field(prefix, None, vmin_vmax=(0, 1), cmap=unsigned_cmap)
    generate_info(beta, "Z", "")
    beta.set_axes_args(title="$2 \\Pi_{rr}(x,\,y) / B_0^2$")

    prr_i = particles_field("Ions", pressures["Prr"], "Z")
    prr_e = particles_field("Electrons", pressures["Prr"], "Z")

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    plot_spatial_beta()
