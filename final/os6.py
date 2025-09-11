#!/usr/bin/env python3

from final import *

ncols=2
nrows=2

fig = plt.figure(figsize=(8 * ncols, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1.1]*nrows, figure=fig)

nmap = (0, 3.0)
jmap = (-0.01, +0.01)


TAU = [4, int(time / tau)]

for i, tt in enumerate(TAU):
  t = int(tt * tau / dts)

  ne = particles_field("Electrons",  "Density", "X", subplot(fig, gs, i, 0), "$n_e$", nmap, unsigned_cmap)
  ne.data = get_parsed_scalar(ne, t) * (-1)
  ne.draw(add_cbar=True)
  ne.draw_info()

  jez = particles_field("Electrons", "Current", "X", subplot(fig, gs, i, 1), "$J_z^e$", jmap)
  jez.data = get_parsed_field(jez, "E", "X", "z", t)
  jez.draw(add_cbar=True)
  jez.draw_info()

fig.tight_layout()
fig.savefig(f"{res_dir}/os6.pdf")
