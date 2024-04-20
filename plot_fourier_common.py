#!/usr/bin/python3

from lib_common import *

names = {
    "Bz":   ("\\delta B_z", 0.02,  150),
    "Er":   ("E_r",         0.02,  80),
    "Ea":   ("E_{\\phi}",   0.02,  80),
    "Jr_e": ("J_r^e",       0.02,  100),
    "Ja_e": ("J_{\\phi}^e", 0.02,  100),
    "Jr_i": ("J_r^i",       0.002, 5),
    "Ja_i": ("J_{\\phi}^i", 0.002, 5),
}

r_range = reduce_array(np.arange(1.00, 3.76, 0.25), rank, proc)

tmin_data = int(3.245 * tau / dts)
tmin = int(3.5 * tau / dts) 
tmax = int(6.0 * tau / dts)

f_tmin = tmin * dts / tau
f_tmax = tmax * dts / tau
args_phit = {
    "ylabel": "$t,~\\tau$",
    "ylim": (f_tmin, f_tmax),
    "yticks": np.linspace(f_tmin, f_tmax, 6),

    "xlabel": "$\\phi,~{\\rm rad}$",
    "xlim": (0.0, 2 * np.pi),
    "xticks": np.linspace(0.0, 2 * np.pi, 5),
    "xticklabels": ["$0$", "$\\pi / 2$", "$\\pi$", "$3 \\pi / 2$", "$2 \\pi$"],
}

mmax = 50
wmax = 0.03
args_mw = {
    "ylabel": "$\\omega,~\\omega_{pe}$",
    "ylim": (-wmax, +wmax),
    "yticks": np.linspace(-wmax, +wmax, 5),

    "xlabel": "$m,~{\\rm units}$",
    "xlim": (-mmax, +mmax),
    "xticks": np.linspace(-mmax, +mmax, 5),
}

m0 = 3
w0 = 0.68
Omega_i = B0 / mi_me
Omega_e = B0
omega0_pi = 1 / np.sqrt(mi_me)

print(f"w: {w0:4f} [1/tau], {w0 / tau:4f} [omega_pe], {w0 / (Omega_i * tau):4f} [Omega_i]")

def draw_common_lines_mw(ax):
    ax.plot([-m0, -m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=0.5)
    ax.plot([+m0, +m0], [-wmax, +wmax], linestyle="--", color="r", linewidth=0.5)

    for m in range(1, 5):
        ax.plot([-mmax, +mmax], [+m * Omega_i, +m * Omega_i], linestyle="--", color="r", linewidth=0.5)

    ax.plot([-mmax, +mmax], [+Omega_e, +Omega_e], linestyle="--", color="r", linewidth=0.5)
    ax.plot([-mmax, +mmax], [+0.0, +0.0], linestyle="--", color="r", linewidth=0.5)

    # drift
    ms = np.linspace(-mmax, +mmax, 100)

    # we're drawing spectrum in untits of (omega_pi / c, omega_pe)
    # ax.plot(ms, ms * (+0.419) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, ms * (+0.200) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, ms * (+1.505) * Omega_i, linestyle='-', c='r', linewidth=0.8, alpha=1)
    # ax.plot(ms, (ms / r) * omega0_pi * np.sqrt(T_i / mi_me), linestyle='--', c='r', linewidth=0.8, alpha=1)

def prepare_field_phit(ax, name, title, r, max_phit):
    F_at = Field(f"{params_path}/Data1D/{name}_phi_t_r={r:.2f}", ax, None, signed_cmap, (-max_phit, +max_phit))
    F_at.set_axes_args(title=title, **args_phit)

    F_at.data = np.load(f"{F_at.path_to_file}.npy", allow_pickle=True)
    F_at.data = F_at.data[(tmin - tmin_data):(tmax - tmin_data),:]
    F_at.boundaries = (0.0, 2 * np.pi, f_tmin, f_tmax)

    if ("Bz" in F_at.path_to_file):
        F_at.data -= F_at.data[0,:]

    return F_at

def prepare_field_mw(ax, title, r, max_mw, data):
    F_mw = Field(None, ax, None, signed_cmap, (0.0, +max_mw))
    F_mw.set_axes_args(title=title, **args_mw)

    F_mw.data, w, k = fourier_transform(data)
    m = np.round(k)

    # the minimum of all min(m) is 56
    F_mw.boundaries = (min(m), max(m), min(w), max(w))
    return F_mw, w, m
