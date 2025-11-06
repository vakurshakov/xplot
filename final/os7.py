#!/usr/bin/env python3

from final import *
from xplot.fourier.fourier import *

ncols=1
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1]*nrows, figure=fig)
ax = subplot(fig, gs, 0, 0)

compare_models = [
    ("Damping", "../T11_MergeV"),
    ("Conductive, circle", "../ConductiveWall_FixRotor/Circle"),
    ("Conductive, square", "../ConductiveWall_FixRotor/Square"),
]

for (label, path) in compare_models:
    F_at = prepare_field_phit(None, "ni", "", 10, 1.0, path)
    F_mw, w, m = prepare_field_mw(ax, "", 10, 1200, F_at.data)
    F_mw.data = np.abs(F_mw.data)

    def mean(d):
        wc=len(w)//2
        ww=3
        return d[-1, :] # np.mean(d[wc-ww:wc+ww,:], axis=0)

    ax.plot(m, mean(F_mw.data), label=label, linewidth=2)

set_big(28)
set_smol(26)
set_ssmol(24)

mmax = 8
ymax = 1500
F_mw.axes_args["title"] = "$\\delta n_i(m, \\, \\omega = 0, \\, r = 10)$"
F_mw.axes_args["xlabel"] = "$m,~{\\rm units}$"
F_mw.axes_args["xlim"] = (-mmax, +mmax)
F_mw.axes_args["xticks"] = np.linspace(-mmax, +mmax, 9)
F_mw.axes_args.pop("ylabel")
F_mw.axes_args.pop("ylim") # (0, ymax)
F_mw.axes_args.pop("yticks") # np.linspace(0, ymax, 6)
F_mw.draw_info()

ax.legend(fontsize=ssmol * 0.64, loc="lower left", framealpha=1.0)
ax.grid(alpha=0.6)

fig.tight_layout()
fig.savefig(f"{res_dir}/os7.pdf")