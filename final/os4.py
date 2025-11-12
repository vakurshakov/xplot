#!/usr/bin/env python3

from final import *

set_big(big*1.1)
set_smol(smol*1.1)
set_ssmol(ssmol*1.1)

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
bpmap = (0, 1)

bz = [
    magnetic_field("Z", subplot(fig, gs, 0, 0), "", bzmap, unsigned_cmap),
    magnetic_field("Z", subplot(fig, gs, 1, 0), "", bzmap, unsigned_cmap),
    magnetic_field("Z", subplot(fig, gs, 2, 0), "", bzmap, unsigned_cmap),
]

ns = [
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 0, 1), "$n_i$", nmap, unsigned_cmap),
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 1, 1), "$n_i$", nmap, unsigned_cmap),
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 2, 1), "$n_i$", nmap, unsigned_cmap),
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

update_args(bzl, "\\rm Magnetic field, $B_z(y)$", bzmap)
update_args(bpl, "\\rm Plasma beta, $\\beta(y)$", bpmap)

pr_i = particles_field("Ions", "Prr", "Z")
pa_i = particles_field("Ions", "Ppp", "Z")
pr_e = particles_field("Electrons", "Prr", "Z")
pa_e = particles_field("Electrons", "Ppp", "Z")
pr_i.data = [None]*len(bz)
pa_i.data = [None]*len(bz)
pr_e.data = [None]*len(bz)
pa_e.data = [None]*len(bz)

TAU = [4, 7, int(time / tau)]

def line_average(d):
    return np.mean(d[:,xc-xw:xc+xw], axis=1)

bz0 = line_average(get_parsed_field(bz[0], "B", "Z", "z", 0))
nmax = [None]*len(bz)

phi = np.linspace(0, 2 * np.pi)
xr = rr * np.cos(phi)
yr = rr * np.sin(phi)

for i, tt in enumerate(TAU):
    t = int(tt * tau / dts)

    bz[i].axes_args["title"] = f"$B_z(x,\\,y,\\,t/\\tau = {tt})$"
    ns[i].axes_args["title"] = f"$n_i(x,\\,y,\\,t/\\tau = {tt})$"

    bz[i].data = get_parsed_field(bz[i], "B", "Z", "z", t) 
    ns[i].data = get_parsed_scalar(ns[i], t) 
    nmax[i] = np.argmax(ns[i].data)

    pd =  line_average(get_parsed_scalar(pr_i, t)) + line_average(get_parsed_scalar(pa_i, t))
    pd += line_average(get_parsed_scalar(pr_e, t)) + line_average(get_parsed_scalar(pa_e, t))
    pd /= 2

    bzl.data[i] = line_average(bz[i].data)
    bpl.data[i] = pd / (np.square(bz0) / 2)

    def draw(d):
        d.draw(add_cbar=True)
        d.draw_info()
        ax = d.axes_position
        ax.set_aspect(1)
        ax.plot(xr, yr, linewidth=3, linestyle="--", c="black")
        
        y0, x0 = np.unravel_index(nmax[i], ns[i].data.shape)
        x0 = (x0 - xc) * dx
        y0 = (y0 - xc) * dx
        ax.scatter(x0, y0, c="black", s=150*1.6, zorder=10)
        ax.scatter(x0, y0, c="red",   s=150*1.0, zorder=10)

    draw(bz[i])
    draw(ns[i])

for i, tt in enumerate(TAU):
    bzl.axes_position.plot(ys, bzl.data[i], label=f"$t/\\tau = {tt}$", linewidth=3)
    bpl.axes_position.plot(ys, bpl.data[i], label=f"$t/\\tau = {tt}$", linewidth=3)

for d in [ bzl, bpl ]:
    d.draw_info()
    d.axes_position.grid(alpha=0.6)

bzl.axes_position.legend(fontsize=Fonts.ssmol, framealpha=1, loc="lower right")
bpl.axes_position.legend(fontsize=Fonts.ssmol, framealpha=1, bbox_to_anchor=(0.76,0.985))

dlmap = [
    (bz[0], "a"),
    (bz[1], "b"),
    (bz[2], "c"),

    (ns[0], "d"),
    (ns[1], "e"),
    (ns[2], "f"),

    (bzl, "g"),
    (bpl, "h"),
]

for (d, l) in dlmap:
    annotate_x(d.axes_position, f"\\rm {l}", -0.1, 1.15, size=Fonts.big)

fig.tight_layout()
fig.savefig(f"{res_dir}/os4.pdf")
