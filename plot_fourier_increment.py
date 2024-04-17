#!/usr/bin/python3

from plot_fourier_common import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"{params_path}/Other"
mkdir(res_dir)

r_range = [2.25, 2.5, 2.75, 3.0] # reduce_array(np.arange(1.00, 3.76, 0.25), rank, proc)

m0 = 3
w0 = 0.68
print(f"w: {w0:4f} [1/tau], {w0 / tau:4f} [omega_pe], {w0 / (Omega_i * tau):4f} [Omega_i]")
        
names = {
    "Er": ("E_r", 1e-5),
    # "Ea": ("E_{\\phi}", 1e-5),
}

arg_1d = {
    "xlim": (f_tmin, f_tmax),
    "xticks": np.linspace(f_tmin, f_tmax, 5),
}

for name, (title, max_at) in names.items():
    for r in r_range:
        # Original map
        F_at = prepare_field_phit(subplot(fig, gs, 0, 0), name, f"$|{title}(\\phi,\,t,\,r = {r:.2f})|^2$", r, 0)
        f_data = F_at.data
        F_at.data = np.square(f_data)

        F_at.cmap = unsigned_cmap
        F_at.vmin_vmax = (0, max_at)
        F_at.draw(add_cbar=True)
        F_at.draw_info()
        

        # Clearing and making filtered map
        F_at_filtered = prepare_field_phit(subplot(fig, gs, 1, 0), name, f"$|{title}(\\phi,\,t,\,r = {r:.2f})|^2,~m = {m0:0d}$", r, 0)

        F_mw_complex, w, m = fourier_transform(f_data)
        F_mw_complex[:, np.where(np.abs(m) > m0)] = 0.0
        ff_data = inverse_fourier_transform(F_mw_complex)[0]
        F_at_filtered.data = np.square(ff_data)

        F_at_filtered.cmap = unsigned_cmap
        F_at_filtered.vmin_vmax = None # (0, max_at)
        F_at_filtered.draw(add_cbar=True)
        F_at_filtered.draw_info()


        # Energy within field component
        Avg = Field(None, subplot(fig, gs, 2, 0), None, unsigned_cmap, (0, max_at))
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
