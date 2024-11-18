#!/usr/bin/env python3

from plot.fourier_common import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"{params_path}/Spectras_filtered"
mkdir(res_dir)

Omega_i = B0 / mi_me
m0 = 3
w0 = 0.68
print(f"w: {w0:4f} [1/tau], {w0 / tau:4f} [omega_pe], {w0 / (Omega_i * tau):4f} [Omega_i]")

def at_lines(ax):
    phi = np.linspace(0, 2 * np.pi)
    const = -10
    ax.plot(phi, (- m0 * phi - const) / w0, linestyle='--', c='r', linewidth=0.8, alpha=1)

def kw_lines(ax):
    draw_common_lines_mw(ax)
    ax.plot([-mmax, +mmax], [w0 / tau, w0 / tau], linestyle="-", color="r", linewidth=0.5)


for name, (title, max_at, max_mw) in names.items():
    for r in r_range:
        F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"${title}(\\phi,\,t,\,r = {r:.2f})$", r, max_at)
        F_at.draw(add_cbar=True)
        F_at.draw_info()

        # Creating Fourier
        F_mw, w, m = prepare_field_mw(subplot(fig, gs, 2, 0), f"${title}(m,\,\\omega,\,r = {r:.2f})$", r, max_mw, F_at.data)

        F_mw_complex = F_mw.data
        F_mw_complex[:, np.where(np.abs(m) > m0)] = 0.0
        F_mw.data = np.abs(F_mw_complex)

        F_mw.draw(add_cbar=True)
        F_mw.draw_info()

        # Clearing and making filtered map
        F_at_filtered = prepare_field_phit(subplot(fig, gs, 1, 0), name, f"${title}(\\phi,\,t,\,r = {r:.2f}),~m = {m0:0d}$", r, max_at * 0.5)
        F_at_filtered.data = inverse_fourier_transform(F_mw_complex)[0]
        F_at_filtered.draw(add_cbar=True)
        F_at_filtered.draw_info()

        kw_lines(F_mw.axes_position)
        at_lines(F_at_filtered.axes_position)

        fig.tight_layout()

        filename = f"{res_dir}/{name}_r={r:.2f}"
        print(filename)

        fig.savefig(f"{filename}.png")
        fig.clear()
