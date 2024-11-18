#!/usr/bin/env python3

from lib_common import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"./{params_path}/Other"
mkdir(res_dir)

Omega_i = B0 / mi_me
v_e = Field(f"./{params_path}/Data1D/v_e_r_t", subplot(fig, gs, 0, 0), (0, data_shape[0] // 2 * dx, 0, time / tau), signed_cmap, (-1, +1))
v_i = Field(f"./{params_path}/Data1D/v_i_r_t", subplot(fig, gs, 1, 0), (0, data_shape[0] // 2 * dx, 0, time / tau), signed_cmap, (-1, +1))

axes_args = {
    "xlabel": "$r,~c / \\omega_{pi}$",
    "ylabel": "$t,~\\tau$",
}

v_e.set_axes_args(title="$\\langle \\omega_e \\rangle (r, t), \\Omega_i$", **axes_args)
v_i.set_axes_args(title="$\\langle \\omega_i \\rangle (r, t), \\Omega_i$", **axes_args)

# tau
tmin = 0.30
tmax = 2.00

# c / omega_pi
wx = 0.05
x0 = 3.00
xmin = x0 - wx
xmax = x0 + wx

fig.text(0.066, 0.85, "$\\langle \\omega_s \\rangle = \\langle J_s / n_s \\rangle_{\\phi} / r$", size=15)

for d in [v_e, v_i]:
    d.data = np.load(f"{d.path_to_file}.npy", allow_pickle=True)

    rs = (np.arange(d.data.shape[0]) + 0.1) * dx
    d.data /= rs[:,np.newaxis]
    d.data /= Omega_i

    d.draw(add_cbar=True)
    d.draw_info()

    ax = d.axes_position

    def plot_rectangle(xmin, xmax, ymin, ymax):
        ax.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin], c="black", linewidth=0.8)

    plot_rectangle(xmin, xmax, tmin, tmax)

    selected = d.data[int(tmin * tau / dts):int(tmax * tau / dts), int(xmin / dx):int(xmax / dx)]
    v_max = np.max(selected)
    v_min = np.min(selected)
    v_avg = np.mean(selected)

    ax.text(xmin + 0.1, tmax + 0.25, "${\\rm min} = " f"{v_min:.3f}$")
    ax.text(xmin + 0.1, tmax + 0.50, "${\\rm max} = " f"{v_max:.3f}$")
    ax.text(xmin + 0.1, tmax + 0.75, "${\\rm avg} = " f"{v_avg:.3f}$")


fig.tight_layout()

fig.savefig(f"{res_dir}/velocities.png")
fig.clear()
