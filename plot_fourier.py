#!/usr/bin/python3

from lib_common import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"./{params_path}/Spectras"
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

def fourier_transform_w(data):
    f_data = np.fft.ifftshift(data)
    f_data = np.fft.fft2(f_data)
    f_data = np.fft.fftshift(f_data)

    # we seek for decomposition by exp(-iwt + ikx)
    f_data = np.abs(f_data)

    # np.fft uses linear frequenses, so we multiply the arrays by (2 * np.pi)
    shape = data.shape
    w = np.fft.fftfreq(shape[0], d=dts) * (2 * np.pi)
    k = np.fft.fftfreq(shape[1], d=(2 * np.pi /  shape[1])) * (2 * np.pi)

    w = np.fft.fftshift(w)
    k = np.fft.fftshift(k)

    return f_data, w, k
        
def draw_helping_lines(ax):
    m0 = 3
    ax.plot([-m0, -m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=0.5)
    ax.plot([+m0, +m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=0.5)
    
    Omega_i = B0 / mi_me
    Omega_e = B0

    for m in range(1, 5):
        ax.plot([-mmax, +mmax], [+m * Omega_i, +m * Omega_i], linestyle="--", color="r", linewidth=0.5)

    ax.plot([-mmax, +mmax], [+Omega_e, +Omega_e], linestyle="--", color="r", linewidth=0.5)
    ax.plot([-mmax, +mmax], [+0.0, +0.0], linestyle="--", color="r", linewidth=0.5)

    # drift 
    omega0_pi = 1 / np.sqrt(mi_me)
    ms = np.linspace(-mmax, +mmax, 100)

    # we're drawing spectrum in untits of (omega_pi / c, omega_pe)
    # ax.plot(ms, ms * (+0.419) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, ms * (+0.200) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, ms * (+1.505) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, (ms / r) * omega0_pi * np.sqrt(T_i / mi_me), linestyle='--', c='r', linewidth=0.8, alpha=1)

for name, (title, at_max, mw_max) in name.items():
    for r in r_range:
        F_at = Field(f"./{params_path}/Data1D/{name}_phi_t_r={r:.2f}", subplot(fig, gs, 0, 0), None, signed_cmap, (-at_max, +at_max))
        F_at.set_axes_args(
            title=f"${title}(\\phi,\,t,\,r = {r:.2f})$",
            ylabel="$t,~\\tau$",
            xlabel="$\\phi,~{\\rm rad}$",
        )

        F_at.data = np.load(f"{F_at.path_to_file}.npy", allow_pickle=True)

        tmin = int(4 * tau / dts)
        tmax = int(9 * tau / dts)
        F_at.data = F_at.data[tmin:tmax,:]

        tmin *= dts / tau
        tmax *= dts / tau
        F_at.boundaries = (0.0, 2 * np.pi, tmin, tmax)
        F_at.set_axes_args(
            xlim=(0.0, 2 * np.pi),
            ylim=(tmin, tmax),
            xticks=np.linspace(0.0, 2 * np.pi, 5),
            yticks=np.linspace(tmin, tmax, 5),
            xticklabels=["$0$", "$\\pi / 2$", "$\\pi$", "$3 \\pi / 2$", "$2 \\pi$"],
        )

        if ("Bz" in F_at.path_to_file):
            F_at.data -= F_at.data[0,:]

        F_at.draw(add_cbar=True)
        F_at.draw_info()

        # Creating Fourier
        F_mw = Field(None, subplot(fig, gs, 1, 0), None, signed_cmap, (0.0, +mw_max))
        F_mw.set_axes_args(
            title=f"${title}(m,\,\\omega,\,r = {r:.2f})$",
            ylabel="$\\omega,~\\omega_{pe}$",
            xlabel="$m,~{\\rm units}$",
        )

        F_mw.data, w, k = fourier_transform_w(F_at.data)
        m = np.round(k)

        F_mw.boundaries = (min(m), max(m), min(w), max(w))

        # the minimum of all min(m) is 56
        mmax = 50
        wmax = 0.03
        F_mw.set_axes_args(
            xlim=(-mmax, +mmax),
            ylim=(-wmax, +wmax),
            xticks=np.linspace(-mmax, +mmax, 5),
            yticks=np.linspace(-wmax, +wmax, 5),
        )

        F_mw.draw(add_cbar=True)
        F_mw.draw_info()

        draw_helping_lines(F_mw.axes_position)
        
        fig.tight_layout()

        filename = f"{res_dir}/{name}_r={r:.2f}"
        print(filename)

        fig.savefig(f"{filename}.png")
        fig.clear()
