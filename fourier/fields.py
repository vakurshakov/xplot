#!/usr/bin/env python3

from fourier import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"{params_path}/Spectra"
mkdir(res_dir)

for name, (title, max_at, max_mw) in named_props:
    for r in r_range:
        F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"${title}(\\phi, t, r = {r:.2f})$", r, max_at)
        F_at.draw(add_cbar=True)
        F_at.draw_info()

        F_mw, w, m = prepare_field_mw(subplot(fig, gs, 1, 0), f"${title}(m, \\omega, r = {r:.2f})$", r, max_mw, F_at.data)
        F_mw.data = np.abs(F_mw.data)
        F_mw.draw(add_cbar=True)
        F_mw.draw_info()

        draw_common_lines_mw(F_mw.axes_position, r)

        filename = f"{res_dir}/{name}_r={r:.2f}"
        print(filename)

        fig.tight_layout()
        fig.savefig(f"{filename}.png")
        fig.clear()
