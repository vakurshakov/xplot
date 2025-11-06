#!/usr/bin/env python3

from plot import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"{params_path}/Other"
mkdir(res_dir)

flux = Field(None, subplot(fig, gs, 0, 0))

flux.axes_args["title"] = "$\\delta \\Phi = \\int \\! B_z(t) \\, dS - \\int \\! B_z(0) \\, dS$"
flux.axes_args["xlabel"] = "time, $t$"
flux.axes_args["xlim"] = (0, 3)
flux.axes_args["xticks"] = np.linspace(0, 3, 5)

flux.draw_info()
ax = flux.axes_position

def draw(path, label):
    data = np.load(f"{path}/Collection/flux_t.npy", allow_pickle=True)
    # data -= data[0]
    ax.plot(np.arange(len(data)) * dt / tau, data, label=label)

draw("../ConductiveWall_FixRotor/Square", "square")
draw("../ConductiveWall_FixRotor/Circle", "circle")
draw("../t11_np_1000", "damping")

ax.legend(fontsize=ssmol)

fig.tight_layout()
fig.savefig(f"{res_dir}/flux_comparison.png")