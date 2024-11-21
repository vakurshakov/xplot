import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from lib_common import *

named_props = [
    ["b",   ("\\delta B_z", 0.02,  150)],
    ["er",  ("E_r",         0.02,  40 )],
    ["ea",  ("E_{\\phi}",   0.02,  40 )],
    ["jri", ("J_r^i",       0.002, 5  )],
    ["jai", ("J_{\\phi}^i", 0.002, 5  )],
    ["jre", ("J_r^e",       0.02,  100)],
    ["jae", ("J_{\\phi}^e", 0.02,  100)],
]

r0    = 10
rmax  = 15
rstep = 2.5
r_range = reduce_array(np.arange(r0, rmax + rstep, rstep))

data_tmin = int(0 * tau / dts)
d_tmin = int(0 * tau / dts)
d_tmax = int(3 * tau / dts)
f_tmin = d_tmin * dts / tau
f_tmax = d_tmax * dts / tau

mmax = 50
wmax = 0.02

args_phit = {
    "ylabel": "$t,~\\tau$",
    "ylim": (f_tmin, f_tmax),
    "yticks": np.linspace(f_tmin, f_tmax, 7),

    "xlabel": "$\\phi,~{\\rm rad}$",
    "xlim": (0.0, 2 * np.pi),
    "xticks": np.linspace(0.0, 2 * np.pi, 5),
    "xticklabels": ["$0$", "$\\pi / 2$", "$\\pi$", "$3 \\pi / 2$", "$2 \\pi$"],
}

args_mw = {
    "ylabel": "$\\omega,~\\omega_{pe}$",
    "ylim": (-wmax, +wmax),
    "yticks": np.linspace(-wmax, +wmax, 5),

    "xlabel": "$m,~{\\rm units}$",
    "xlim": (-mmax, +mmax),
    "xticks": np.linspace(-mmax, +mmax, 5),
}

Omega_i = B0 / mi_me
Omega_e = B0
omega0_pi = 1 / np.sqrt(mi_me)
print("Omega_i:", Omega_i, "[w_pe]", "Omega_e:", Omega_e, "[w_pe]", "omega0_pi:", omega0_pi, "[w_pe]")

def draw_common_lines_mw(ax, r=r0):
    # we're drawing spectrum in untits of (omega_pi / c, omega_pe)
    text_mult = 1.1
    text_size = ssmol * 0.6

    ms = range(0, 5)
    for m in ms:
        ax.plot([-mmax, +mmax], [m*Omega_i, m*Omega_i], linestyle="--", color="r", linewidth=0.5)

        if (m*Omega_i > wmax and m > ms[1]) or m == ms[-1]:
            ax.text(mmax / 2, m*Omega_i * text_mult, "$m \\Omega_i$", color="r", fontsize=text_size)
            break

    if (Omega_e < wmax):
        ax.plot([-mmax, +mmax], [Omega_e, Omega_e], linestyle="--", color="r", linewidth=0.5)
        ax.text(-mmax / 2, Omega_e * text_mult, "$Omega_e$", color="r", fontsize=text_size)

    # drift
    # ms = np.linspace(-mmax, +mmax, 100)
    # ax.plot(ms, (ms / r) * np.sqrt(T_i / mi_me), linestyle='--', c='r', linewidth=0.8, alpha=1)


def prepare_field_phit(ax, name, title, r, max_phit):
    F_at = Field(f"{params_path}/Collection/{name}_phit_r={r:.2f}", ax, (0.0, 2 * np.pi, f_tmin, f_tmax), signed_cmap, (-max_phit, +max_phit))
    F_at.set_axes_args(title=title, **args_phit)

    F_at.data = np.load(f"{F_at.path_to_file}.npy", allow_pickle=True)

    assert ((d_tmax - d_tmin) < F_at.data.shape[0]), f"Incorrect fourier window is chosen: d_tmin {d_tmin}, d_tmax {d_tmax}, F_at.data.shape {F_at.data.shape[0]}"
    F_at.data = F_at.data[(d_tmin - data_tmin):(d_tmax - data_tmin),:]

    if ("B_z" in F_at.axes_args["title"]):
        F_at.data -= F_at.data[0,:]

    return F_at


def prepare_field_mw(ax, title, r, max_mw, data):
    F_mw = Field(None, ax, None, signed_cmap, (0.0, +max_mw))
    F_mw.set_axes_args(title=title, **args_mw)

    F_mw.data, w, k = fourier_transform(data)
    m = np.round(k)

    F_mw.boundaries = (min(m), max(m), min(w), max(w))
    return F_mw, w, m
