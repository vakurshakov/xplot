#!/usr/bin/env python3

from final import *

ncols=4
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

y_max = data_shape["Z"][1]
xc    = data_shape["Z"][0] // 2
xw    = 2

ys = (np.arange(0, y_max) - y_max // 2) * dy

nmap = (0, 4)
bzmap = (0, 0.2)
bpmap = (0, 1.5)

bz = [
    magnetic_field("Z", subplot(fig, gs, 0, 0), "", bzmap, unsigned_cmap),
    magnetic_field("Z", subplot(fig, gs, 1, 0), "", bzmap, unsigned_cmap),
    magnetic_field("Z", subplot(fig, gs, 2, 0), "", bzmap, unsigned_cmap),
]

s = "Ions"
ms = (mi_me if s == "Ions" else 1)
es = ((+1)  if s == "Ions" else (-1))

ns = [
    particles_field(s, "Density", "Z", subplot(fig, gs, 0, 1), f"$n_{s[0].lower()}$", nmap, unsigned_cmap),
    particles_field(s, "Density", "Z", subplot(fig, gs, 1, 1), f"$n_{s[0].lower()}$", nmap, unsigned_cmap),
    particles_field(s, "Density", "Z", subplot(fig, gs, 2, 1), f"$n_{s[0].lower()}$", nmap, unsigned_cmap),
]

bzl = magnetic_field(         "Z", subplot(fig, gs, 3, 0), "")
bpl = particles_field("", "", "Z", subplot(fig, gs, 3, 1), "")
bzl.data = [None]*len(bz)
bpl.data = [None]*len(bz)

def update_args(d, title, map):
    d.axes_args.update(
        title=title,
        xlabel="$y,~c/\\omega_{pe}$",
        ylabel=None,
        xlim=(-30, 30),
        ylim=map,
        xticks=np.linspace(-30, 30, 5),
        yticks=np.linspace(*map, 5),
    )

update_args(bzl, "Magnetic field, $B_z(y)$", bzmap)
update_args(bpl, "Plasma beta, $\\beta(y)$", bpmap)

pr = particles_field(s, "Prr", "Z")
pa = particles_field(s, "Ppp", "Z")
pr.data = [None]*len(bz)
pa.data = [None]*len(bz)

TAU = [4, 7, int(time / tau)]

def line_average(d):
    return np.mean(d[:,xc-xw:xc+xw], axis=1)

for i, tt in enumerate(TAU):
    t = int(tt * tau / dts)

    bz[i].axes_args["title"] = f"$B_z(x,\\,y,\\,t/\\tau = {tt})$"
    ns[i].axes_args["title"] = f"$n_{s[0].lower()}(x,\\,y,\\,t/\\tau = {tt})$"

    bz[i].data = get_parsed_field(bz[i], "B", "Z", "z", t) 
    ns[i].data = get_parsed_scalar(ns[i], t) 
    pr.data[i] = line_average(get_parsed_scalar(pr, t))
    pa.data[i] = line_average(get_parsed_scalar(pa, t))

    bzl.data[i] = line_average(bz[i].data)
    bpl.data[i] = (pr.data[i] + pa.data[i]) / np.square(bzl.data[i])

for d in [ *bz, *ns ]:
    d.draw(add_cbar=True)
    d.draw_info()
    d.axes_position.set_aspect(1)

for i, tt in enumerate(TAU):
    bzl.axes_position.plot(ys, bzl.data[i], label=f"$t/\\tau = {tt}$", linewidth=3)
    bpl.axes_position.plot(ys, bpl.data[i], label=f"$t/\\tau = {tt}$", linewidth=3)

for d in [ bzl, bpl ]:
    d.draw_info()
    d.axes_position.grid(alpha=0.6)

bzl.axes_position.legend(fontsize=ssmol, framealpha=1, loc="lower right")
bpl.axes_position.legend(fontsize=ssmol, framealpha=1, bbox_to_anchor=(0.76,0.985))

fig.tight_layout()
fig.savefig(f"{res_dir}/os4.pdf")
