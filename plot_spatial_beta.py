#!/usr/bin/python3

from lib_common import *
from plot_fields import *
from plot_particles import *

B0 = 0.191

def plot_spatial_beta(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    bz.data = get_parsed_field(bz, "B", "Z", "z", t)
    prr_i.data = get_parsed_scalar(prr_i, t)
    prr_e.data = get_parsed_scalar(prr_e, t)

    beta.data = 2 * (prr_i.data + prr_e.data) / (B0 * B0)

    fig.suptitle(f"$t / \\tau = {t * dts / tau:.3f}$", x=0.51, bbox=bbox, fontsize=smol)

    bz.draw(add_cbar=True)
    bz.draw_info()
    bz.axes_position.set_aspect(1)

    beta.draw(add_cbar=True)
    beta.draw_info()
    beta.axes_position.set_aspect(1)

    fig.tight_layout()
    plt.show()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

    for diag in [bz, beta]:
        diag.clear()

if __name__ == "__main__":
    ncols=2
    nrows=1

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    bz = magnetic_field("Z", subplot(fig, gs, 0, 0), "")
    bz.set_axes_args(title="$B_z(x,\,y)$")

    beta = Field(None, subplot(fig, gs, 1, 0), cmap=unsigned_cmap)
    generate_info(beta, "Z", "")
    beta.set_axes_args(title="$2 \\Pi_{rr}(x,\,y) / B_0^2$")

    prr_i = particles_field("Ions", pressures["Prr"], "Z")
    prr_e = particles_field("Electrons", pressures["Prr"], "Z")

    res_dir = f"{params_path}/Spatial_beta"
    mkdir(res_dir)

    offset = 5
    t0 = 0
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range:
        plot_spatial_beta(t)
