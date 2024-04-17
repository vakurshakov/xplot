#!/usr/bin/python3

from plot_fourier_common import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"{params_path}/Spectras"
mkdir(res_dir)

for name, (title, max_at, max_mw) in name.items():
    for r in r_range:
        F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"${title}(\\phi,\,t,\,r = {r:.2f})$", r, max_at)
        F_at.draw(add_cbar=True)
        F_at.draw_info()

        # Creating Fourier
        F_mw, w, m = prepare_field_mw(subplot(fig, gs, 1, 0), f"${title}(m,\,\\omega,\,r = {r:.2f})$", r, max_mw, F_at.data)
        F_mw.data = np.abs(F_mw.data)

        F_mw.draw(add_cbar=True)
        F_mw.draw_info()

        draw_common_lines_mw(F_mw.axes_position)
        
        fig.tight_layout()

        filename = f"{res_dir}/{name}_r={r:.2f}"
        print(filename)

        fig.savefig(f"{filename}.png")
        fig.clear()
