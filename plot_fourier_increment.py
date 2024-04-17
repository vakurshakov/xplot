#!/usr/bin/python3

from lib_common import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"./{params_path}/Other"
mkdir(res_dir)

r_range = [2.25, 2.5, 2.75, 3.0] # reduce_array(np.arange(1.00, 3.76, 0.25), rank, proc)

tmin = int(4 * tau / dts)
tmax = int(9 * tau / dts)

f_tmin = tmin * dts / tau
f_tmax = tmax * dts / tau

# the minimum of all min(m) is 56
mmax = 50
wmax = 0.03

Omega_i = B0 / mi_me

m0 = 3
w0 = 0.68
print(f"w: {w0:4f} [1/tau], {w0 / tau:4f} [omega_pe], {w0 / (Omega_i * tau):4f} [Omega_i]")
        
names = {
    "Ea":   ("E_{\\phi}",   1e-5),
}

arg_2d = {
    "xlabel": "$\\phi,~{\\rm rad}$",
    "xlim": (0.0, 2 * np.pi),
    "xticks": np.linspace(0.0, 2 * np.pi, 5),
    "xticklabels": ["$0$", "$\\pi / 2$", "$\\pi$", "$3 \\pi / 2$", "$2 \\pi$"],

    "ylabel": "$t,~\\tau$",
    "ylim": (f_tmin, f_tmax),
    "yticks": np.linspace(f_tmin, f_tmax, 5),
}

arg_1d = {
    "xlim": (f_tmin, f_tmax),
    "xticks": np.linspace(f_tmin, f_tmax, 5),
}

for name, (title, at_max) in names.items():
    for r in r_range:
        # Original map
        F_at = Field(f"./{params_path}/Data1D/{name}_phi_t_r={r:.2f}", subplot(fig, gs, 0, 0), None, unsigned_cmap, (0, at_max))
        F_at.set_axes_args(title=f"$|{title}(\\phi,\,t,\,r = {r:.2f})|^2$", **arg_2d)
        F_at.boundaries = (0.0, 2 * np.pi, f_tmin, f_tmax)

        f_data = np.load(f"{F_at.path_to_file}.npy", allow_pickle=True)
        f_data = f_data[tmin:tmax,:]
        F_at.data = np.square(f_data)

        F_at.draw(add_cbar=True)
        F_at.draw_info()
        
        # Clearing and making filtered map
        F_at_filtered = Field(None, subplot(fig, gs, 1, 0), None, unsigned_cmap, (0, 1e-7))
        F_at_filtered.set_axes_args(title=f"$|{title}(\\phi,\,t,\,r = {r:.2f})|^2,~m = {m0:0d}$", **arg_2d)
        F_at_filtered.boundaries = (0.0, 2 * np.pi, f_tmin, f_tmax)

        F_mw_complex, w, k = fourier_transform(f_data)
        F_mw_complex[:, np.where(np.abs(k - m0) > 0.5)] = 0.0
        ff_data = inverse_fourier_transform(F_mw_complex)[0]
        F_at_filtered.data = np.square(ff_data)

        F_at_filtered.draw(add_cbar=True)
        F_at_filtered.draw_info()

        # Energy within field component
        Avg = Field(None, subplot(fig, gs, 2, 0), None, unsigned_cmap, (0, at_max))
        Avg.set_axes_args(
            title=f"$\\ln(\\langle |{title}(\\phi,\,t,\,r = {r:.2f})|^2 \\rangle)$",
            ylim=(-21, -11),
            yticks=np.linspace(-21, -11, 5),
            **arg_1d,
        )

        ts = np.arange(tmin, tmax) * dts / tau
        
        f_avg = np.log(np.mean(np.square(f_data), axis=1))
        Avg.axes_position.plot(ts, f_avg, label="$\\forall m$")
        
        ff_avg = np.log(np.mean(np.square(ff_data), axis=1))
        Avg.axes_position.plot(ts, ff_avg, label="$m = 3$")

        Avg.draw_info()
        Avg.axes_position.legend(loc="upper right", fontsize=25)

        # Util
        phi = np.linspace(0, 2 * np.pi)
        const = -11
        F_at_filtered.axes_position.plot(phi, (- m0 * phi - const) / w0, linestyle='--', c='r', linewidth=2, alpha=1)
        F_at_filtered.axes_position.text(np.pi, 4.5, f"$\\omega = {w0 / (Omega_i * tau):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)

        if (name == "Ea"):
            y0 = -20.2
            x0 = f_tmin
            gamma = 0.6
            Avg.axes_position.plot(ts, -20.2 + 2 * gamma * (ts - f_tmin), linestyle='--', c='r', linewidth=2)
            Avg.axes_position.text(6.0, -15.5, f"$\\Gamma = {(gamma / (Omega_i * tau)):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)
        if (name == "Er"):
            gamma = 0.35
            Avg.axes_position.plot(ts, -19.8 + 2 * gamma * (ts - f_tmin), linestyle='--', c='r', linewidth=2)
            Avg.axes_position.text(5.5, -16, f"$\\Gamma = {(gamma / (Omega_i * tau)):.2f}\,\\Omega_i$", fontsize=smol, bbox=bbox)

        fig.tight_layout()
        fig.tight_layout()

        filename = f"{res_dir}/Increment_{name}_r={r:.2f}"
        print(filename)

        fig.savefig(f"{filename}.png")
        fig.clear()
