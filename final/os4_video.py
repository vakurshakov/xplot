#!/usr/bin/env python3

from final import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.0, 8 * nrows * 1.02))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

bzmap = (0, 0.2)
nmap = (0, 4)

bz = magnetic_field("Z", subplot(fig, gs, 0, 0), "$B_z$", bzmap, unsigned_cmap)
ni = particles_field("Ions", "Density", "Z", subplot(fig, gs, 1, 0), "$n_i$", nmap, unsigned_cmap)

bz.axes_position.set_aspect(1)
ni.axes_position.set_aspect(1)

res_dir = f"{params_path}/Final/os4_video"
mkdir(res_dir)

offset = 10
t_range = create_t_range(0, int(time / dts), offset)

for t in t_range:
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        continue

    bz.data = get_parsed_field(bz, "B", "Z", "z", t) 
    ni.data = get_parsed_scalar(ni, t) 

    for d in [bz, ni]:
        d.draw(add_cbar=True)
        d.draw_info()

    fig.suptitle(f"$t / \\tau = {t * dts / tau:.3f}$", x=0.53, y=0.98, bbox=bbox, fontsize=big)
    
    if t == t_range[0]:
        fig.tight_layout(w_pad=3, rect=(0, 0, 1, 0.99))

    fig.savefig(filename)

    for d in [bz, ni]:
        d.clear()