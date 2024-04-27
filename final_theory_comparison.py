#!/usr/bin/python3

from final_common import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.18, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)


arg_tc = {
    "xlabel": "$t, \\tau$",
    "xlim": (0, 9),
    "xticks": np.linspace(0, 9, 7),
}

b_t = Field(None, subplot(fig, gs, 0, 0))
b_t.set_axes_args(
    title="${\\rm Magnetic~field},~b$",
    **arg_tc,
)

n_t = Field(None, subplot(fig, gs, 1, 0))
n_t.set_axes_args(
    title="${\\rm Ion~density},~n_i$",
    ylim=(0, 2),
    **arg_tc,
)

T_t = Field(None, subplot(fig, gs, 2, 0))
T_t.set_axes_args(
    title="${\\rm Ion~temperature},~T_i$",
    ylabel="{\\rm keV}",
    **arg_tc,
)

b = np.load(f"{params_path}/Data1D/b_t.npy")
pr = np.load(f"{params_path}/Data1D/pri_t.npy")
pa = np.load(f"{params_path}/Data1D/pai_t.npy")
ni = np.load(f"{params_path}/Data1D/ni_t.npy")

Ti = (pr + pa) / 2
Ti = np.divide(Ti, ni, out=np.zeros_like(Ti), where=(np.abs(ni) > 1e-3)) * ((B0 * B0 / 2) * 511)  # to keV

b_t.data = b
T_t.data = Ti
n_t.data = ni

ts = np.arange(0, int(9 * tau / dts) + 1) * dts / tau
_ts = sliding_average(ts, 8)

lw = 3

for diag in [b_t, n_t, T_t]:
    diag.axes_position.plot(_ts, sliding_average(diag.data, 8), label="${\\rm PIC}$", linewidth=lw)

xi = B0 * B0 / (2 * (T_i + T_e))

b_th = np.exp(- ts / (2 * xi))
b_t.axes_position.plot(ts, b_th, label="${\\rm theory}$", linewidth=lw)

n_th = 2 * xi * (1 - b_th)
n_t.axes_position.plot(ts, n_th, label="${\\rm theory}$", linewidth=lw)

T_th = (T_i * 511) * (1 + b_th) / 2
T_t.axes_position.plot(ts, T_th, label="${\\rm theory}$", linewidth=lw)

for diag in [b_t, n_t, T_t]:
    ax = diag.axes_position
    ax.grid(alpha=0.3)
    diag.draw_info()

for ax, letter in zip(fig.axes, "abc"):
    annotate_x(ax, "${\\rm " + letter + "}$", 0.08, 0.88)

b_t.axes_position.legend(loc="upper right", fontsize=ssmol * 1.1)
T_t.axes_position.legend(loc="upper right", fontsize=ssmol * 1.1)
n_t.axes_position.legend(loc="lower right", fontsize=ssmol * 1.1)

fig.tight_layout()
fig.tight_layout()

fig.savefig(f"{res_dir}/figure6.pdf")
