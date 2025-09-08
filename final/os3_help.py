#!/usr/bin/env python3

from final import *

ncols=3
nrows=2

fig = plt.figure(figsize=(8 * ncols, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.1] * nrows, figure=fig)

pmap = (0,     +0.005)
vmap = (-0.01, +0.01)

s = sorts[1]

pr = particles_field(s, "Prr",     "Y", subplot(fig, gs, 0, 0), "$\\Pi_{rr}$", pmap, unsigned_cmap)
pz = particles_field(s, "Pzz",     "Y", subplot(fig, gs, 0, 1), "$\\Pi_{zz}$", pmap, unsigned_cmap)
jr = particles_field(s, "Current", "Y", subplot(fig, gs, 1, 0), "$v_r$",       vmap)
jz = particles_field(s, "Current", "Y", subplot(fig, gs, 1, 1), "$v_z$",       vmap)
n  = particles_field(s, "Density", "Y")

pcr = Field("", subplot(fig, gs, 2, 0), cmap=unsigned_cmap, vmin_vmax=pmap)
generate_info(pcr, "Y", "$P_{\\perp}$")

pcz = Field("", subplot(fig, gs, 2, 1), cmap=unsigned_cmap, vmin_vmax=pmap)
generate_info(pcz, "Y", "$P_{\\|}$")

TAU = 4
TIME = int(TAU * tau / dts)

pr.data = get_parsed_scalar(pr, TIME)
pz.data = get_parsed_scalar(pz, TIME)
jr.data = get_parsed_field(jr, "E", "Y", "x", TIME)
jz.data = get_parsed_field(jz, "E", "Y", "z", TIME)
n.data  = get_parsed_scalar(n, TIME)

jr.data = np.divide(jr.data, n.data, where=(np.abs(n.data) > 1e-3), out=np.zeros_like(n.data))
jz.data = np.divide(jz.data, n.data, where=(np.abs(n.data) > 1e-3), out=np.zeros_like(n.data))

ms = (mi_me if s == "Ions" else 1)
es = ((+1)  if s == "Ions" else (-1))
pcr.data = pr.data - ms * (n.data / es) * np.square(jr.data)
pcz.data = pz.data - ms * (n.data / es) * np.square(jz.data)

for d in [ pr, pz, jr, jz, pcr, pcz ]:
    d.draw(add_cbar=True)
    d.draw_info()

fig.suptitle(f"{s}, $t / \\tau = {TAU:.3f}$", y=0.99, bbox=bbox, fontsize=big)
fig.tight_layout(rect=(0, 0, 1, 0.99))

fig.savefig(f"{res_dir}/os3_help_{s.lower()}.png")