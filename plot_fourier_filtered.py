#!/usr/bin/python3

from lib_common import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"./{params_path}/Spectras_filtered"
mkdir(res_dir)

name = {
    "Bz":   ("\\delta B_z", 0.02,  100),
    "Er":   ("E_r",         0.02,  100),
    "Ea":   ("E_{\\phi}",   0.02,  100),
    "Jr_e": ("J_r^e",       0.02,  15),
    "Ja_e": ("J_{\\phi}^e", 0.02,  15),
    "Jr_i": ("J_r^i",       0.002, 5),
    "Ja_i": ("J_{\\phi}^i", 0.002, 5),
}

r_range = reduce_array(np.arange(1.00, 3.76, 0.25), rank, proc)

tmin = int(4 * tau / dts)
tmax = int(8 * tau / dts)

f_tmin = tmin * dts / tau
f_tmax = tmax * dts / tau

# the minimum of all min(m) is 56
mmax = 50
wmax = 0.03

Omega_i = B0 / mi_me
m0 = 3
w0 = 0.68
print(f"w: {w0:4f} [1/tau], {w0 / tau:4f} [omega_pe], {w0 / (Omega_i * tau):4f} [Omega_i]")

def at_lines(ax):
    phi = np.linspace(0, 2 * np.pi)
    const = -10

    ax.plot(phi, (- m0 * phi - const) / w0, linestyle='--', c='r', linewidth=0.8, alpha=1)

def kw_lines(ax):
    ax.plot([-m0, -m0], [-wmax, +wmax], linestyle="-", color="r", linewidth=0.5)
    ax.plot([+m0, +m0], [-wmax, +wmax], linestyle="-", color="r", linewidth=0.5)
    ax.plot([-mmax, +mmax], [w0 / tau, w0 / tau], linestyle="-", color="r", linewidth=0.5)


for name, (title, at_max, mw_max) in name.items():
    for r in r_range:
        F_at = Field(f"./{params_path}/Data1D/{name}_phi_t_r={r:.2f}", subplot(fig, gs, 0, 0), None, signed_cmap, (-at_max, +at_max))
        F_at.set_axes_args(
            title=f"${title}(\\phi,\,t,\,r = {r:.2f})$",
            ylabel="$t,~\\tau$",
            xlabel="$\\phi,~{\\rm rad}$",
        )

        F_at.data = np.load(f"{F_at.path_to_file}.npy", allow_pickle=True)
        F_at.data = F_at.data[tmin:tmax,:]
        F_at.boundaries = (0.0, 2 * np.pi, f_tmin, f_tmax)
        F_at.set_axes_args(
            xlim=(0.0, 2 * np.pi),
            ylim=(f_tmin, f_tmax),
            xticks=np.linspace(0.0, 2 * np.pi, 5),
            yticks=np.linspace(f_tmin, f_tmax, 5),
            xticklabels=["$0$", "$\\pi / 2$", "$\\pi$", "$3 \\pi / 2$", "$2 \\pi$"],
        )

        if ("Bz" in F_at.path_to_file):
            F_at.data -= F_at.data[0,:]

        F_at.draw(add_cbar=True)
        F_at.draw_info()

        # Creating Fourier
        F_mw = Field(None, subplot(fig, gs, 2, 0), None, signed_cmap, (0.0, +mw_max))
        F_mw.set_axes_args(
            title=f"${title}(m,\,\\omega,\,r = {r:.2f})$",
            ylabel="$\\omega,~\\omega_{pe}$",
            xlabel="$m,~{\\rm units}$",
        )

        F_mw_complex, w, k = fourier_transform(F_at.data)
        m = np.round(k)

        F_mw.data = np.abs(F_mw_complex)
        F_mw.data[:, np.where(np.abs(m) > m0)] = 0.0

        F_mw.boundaries = (min(m), max(m), min(w), max(w))
        F_mw.set_axes_args(
            xlim=(-mmax, +mmax),
            ylim=(-wmax, +wmax),
            xticks=np.linspace(-mmax, +mmax, 5),
            yticks=np.linspace(-wmax, +wmax, 5),
        )

        F_mw.draw(add_cbar=True)
        F_mw.draw_info()


        # Clearing and making filtered map
        F_at_filtered = Field(None, subplot(fig, gs, 1, 0), None, signed_cmap, (-at_max * 0.8, +at_max * 0.8))
        F_at_filtered.set_axes_args(**F_at.axes_args)
        F_at_filtered.set_axes_args(title=f"${title}(\\phi,\,t,\,r = {r:.2f}),~m = {m0:0d}$")

        F_mw_complex[:, np.where(np.abs(m) > m0)] = 0.0
        F_at_filtered.data = inverse_fourier_transform(F_mw_complex)[0]
        F_at_filtered.boundaries = (0.0, 2 * np.pi, f_tmin, f_tmax)
        F_at_filtered.draw(add_cbar=True)
        F_at_filtered.draw_info()
        
        kw_lines(F_mw.axes_position)
        at_lines(F_at_filtered.axes_position)

        fig.tight_layout()

        filename = f"{res_dir}/{name}_r={r:.2f}"
        print(filename)

        fig.savefig(f"{filename}.png")
        fig.clear()