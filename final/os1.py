#!/usr/bin/env python3

from final import *

set_big(big*1.1)
set_smol(smol*1.1)
set_ssmol(ssmol*1.1)

ncols=4
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1.3, 1], figure=fig)

nmap = (0, 3.0)
bmap = (0, 0.2)
jmap = (-0.01, +0.01)

niy = particles_field("Ions",       "Density", "Y", subplot(fig, gs, 0, 0), "$n_i$", nmap, unsigned_cmap)
niz = particles_field("Ions",       "Density", "Z", subplot(fig, gs, 0, 1), "$n_i$", nmap, unsigned_cmap)

jepy = particles_field("Electrons", "Current", "Y", subplot(fig, gs, 1, 0), "$J_{\\phi}^e$", jmap)
jepz = particles_field("Electrons", "Current", "Z", subplot(fig, gs, 1, 1), "$J_{\\phi}^e$", jmap)
jizy = particles_field("Ions",      "Current", "Y", subplot(fig, gs, 2, 0), "$J_z^i$", jmap)
jipz = particles_field("Ions",      "Current", "Z")

jpzl = particles_field("", "", "Z", subplot(fig, gs, 2, 1), "")
nizl = particles_field("", "", "Z", subplot(fig, gs, 3, 0), "")
bzzl = magnetic_field(         "Z", subplot(fig, gs, 3, 1), "")

jpzl.axes_args["title"] = "\\rm Current density"
nizl.axes_args["title"] = "\\rm Ion density"
bzzl.axes_args["title"] = "\\rm Magnetic field, $B_z$"

TAU = 4
TIME = int(TAU * tau / dts)

niy.data = get_parsed_scalar(niy, TIME)
niz.data = get_parsed_scalar(niz, TIME)

jepy.data = get_parsed_field(jepy, "E", "Y", "y", TIME)
jepz.data = get_parsed_field(jepz, "E", "Z", "y", TIME)[1]
jizy.data = get_parsed_field(jizy, "E", "Y", "z", TIME)
jipz.data = get_parsed_field(jipz, "E", "Z", "y", TIME)[1]

for diag in [niy, niz, jepy, jepz, jizy]:
    if (diag == niz or diag == jepz):
        diag.axes_position.set_aspect(1)
    diag.draw(add_cbar=True)
    diag.draw_info()

rs = np.arange(0, data_shape["Z"][0] // 2) * dx

ax = jpzl.axes_position
jepl = phi_averaged(jepz.data, R_MAP)
jipl = phi_averaged(jipz.data, R_MAP)
ax.plot(rs, jepl,        label="$J_{\\phi}^e$",               linewidth=3)
ax.plot(rs, jipl,        label="$J_{\\phi}^i$",               linewidth=3)
ax.plot(rs, jepl + jipl, label="$J_{\\phi}^e + J_{\\phi}^i$", linewidth=3)
ax.legend(loc="upper right", fontsize=smol * 0.9)

for n in np.arange(0, TAU+1):
    t = int(n * tau / dts) 
    args = dict(linewidth=(3 if t != TIME else 5), linestyle=(None if t != TIME else "--"), zorder=10)
    nizl.axes_position.plot(rs, phi_averaged(get_parsed_scalar(niz, t), R_MAP),                **args)
    bzzl.axes_position.plot(rs, phi_averaged(get_parsed_field(bzzl, "B", "Z", "z", t), R_MAP), **args)

def info(map):
    return dict(
        xlabel="$r,~c/\\omega_{pe}$",
        ylabel=None,
        xlim=(0, 30),
        ylim=map,
        xticks=np.linspace(0, 30, 7),
        yticks=np.linspace(*map, 5),
    )

jpzl.draw_info(**info((jmap)))
nizl.draw_info(**info(nmap))
bzzl.draw_info(**info(bmap))

dlmap = [
    (niy, "a"),
    (niz, "b"),
    (jepy, "c"),
    (jepz, "d"),
    (jizy, "e"),
    (jpzl, "f"),
    (nizl, "g"),
    (bzzl, "h"),
]

for (d, l) in dlmap:
    y = 1.1
    if l in "bdfh": y=1.15
    annotate_x(d.axes_position, f"\\rm {l}", -0.1, y)

fig.tight_layout()
fig.savefig(f"{res_dir}/os1.pdf")
