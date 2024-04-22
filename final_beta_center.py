#!/usr/bin/python3

from final_common import *

ncols=1
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.18, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

set_big(42)
set_smol(40)
set_ssmol(36)

plot = Field(None, subplot(fig, gs, 0, 0), boundaries, unsigned_cmap, (0, +0.02))

plot.set_axes_args(
    title="${\\rm Plasma~beta}$",
    xlabel="$t, \\tau$",
    xlim=(0, 9),
    xticks=np.linspace(0, 9, 7),

    ylim=(0, 1.3),
    # yticks=np.linspace(0, 9, 5),
)

b = np.load(f"{params_path}/Data1D/b_t.npy")
pr = np.load(f"{params_path}/Data1D/pr_t.npy")
pd = np.load(f"{params_path}/Data1D/pd_t.npy")

pb = 1 - np.square(b)

ts = np.arange(0, int(9 * tau / dts) + 1) * dts / tau

w = 10
_ts = sliding_average(ts, w)

ax = plot.axes_position
ax.plot(_ts, sliding_average(pr, w), label="$\\beta_r$")
ax.plot(_ts, sliding_average(pr - pd, w), label="$\\beta_r - \\Delta \\beta$")
ax.plot(ts, pb, label="$\\beta_b$")

tmin = int(3.5 * tau / dts)
tmax = int(6.5 * tau / dts)
ax.fill_between(
    sliding_average(ts[tmin:tmax], w),
    sliding_average(pr[tmin:tmax], w),
    np.zeros_like(sliding_average(pr[tmin:tmax], w)),
    color="grey", 
    alpha=0.2,
    zorder=-1
)

ax.annotate(
    "${\\rm linear~stage~of}~m=3~{\\rm growth}$",
    (4.8, 0.9),
    (3.1, 0.5),
    **arg_anno,
)

ax.legend(fontsize=ssmol, framealpha=1)
ax.grid(alpha=0.3)
plot.draw_info()

fig.tight_layout()
fig.tight_layout()

fig.savefig(f"{res_dir}/figure5.pdf")
