#!/usr/bin/env python3

from plot import *

def plot_forces(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    def time_average(func, t0, dt):
        result = np.zeros_like(func(t0))
        for t in range(t0 - dt, t0):
            result += func(t)
        return result / dt

    n_i.data = phi_averaged(time_average(lambda t0: (+1) * get_parsed_scalar(n_i, t0), t, 10), R_MAP)
    prr_i.data = phi_averaged(time_average(lambda t0: get_parsed_scalar(prr_i, t0), t, 10), R_MAP)
    paa_i.data = phi_averaged(time_average(lambda t0: get_parsed_scalar(paa_i, t0), t, 10), R_MAP)

    n_e.data = phi_averaged(time_average(lambda t0: (-1) * get_parsed_scalar(n_e, t0), t, 10), R_MAP)
    prr_e.data = phi_averaged(time_average(lambda t0: get_parsed_scalar(prr_e, t0), t, 10), R_MAP)
    paa_e.data = phi_averaged(time_average(lambda t0: get_parsed_scalar(paa_e, t0), t, 10), R_MAP)

    ja_i.data = phi_averaged(time_average(lambda t0: get_parsed_field(ja_i, "E", "Z", "", t0), t, 10)[1], R_MAP)
    ja_e.data = phi_averaged(time_average(lambda t0: get_parsed_field(ja_e, "E", "Z", "", t0), t, 10)[1], R_MAP)

    er.data = phi_averaged(time_average(lambda t0: get_parsed_field(er, "E", "Z", "", t0), t, 10)[0], R_MAP)
    bz.data = phi_averaged(time_average(lambda t0: get_parsed_field(bz, "B", "Z", "z", t0), t, 10), R_MAP)

    def avg1(d): return sliding_average(d, 3)
    def avg2(d): return sliding_average(d, 7)

    rs = (np.arange(0, data_shape["Z"][0] // 2) + 0.5) * dx
    avg1_rs = avg1(rs)
    avg2_rs = avg2(rs)

    ax = forces_i.axes_position
    ax.plot(avg1_rs, np.gradient(avg1(prr_i.data), dx), label="$\\partial \\Pi^i_{rr} / \\partial r$", linewidth=2)
    ax.plot(avg1_rs, avg1(np.divide(prr_i.data - paa_i.data, rs)), label="$(\\Pi^i_{rr} - \\Pi^i_{\\phi \\phi}) / r$", linewidth=2)
    ax.plot(avg2_rs, avg2(n_i.data * er.data), label="$n_i E_r$", linewidth=2)
    ax.plot(rs, ja_i.data * bz.data, label="$J^i_{\\phi} B_z$", linewidth=2)
    ax.plot(avg2_rs, avg2(-(np.gradient(prr_i.data, dx) + np.divide(prr_i.data - paa_i.data, rs)) + n_i.data * er.data + ja_i.data * bz.data), label="res. force", linewidth=2, linestyle="--", color="red")

    ax = forces_e.axes_position
    ax.plot(avg1_rs, np.gradient(avg1(prr_e.data), dx), label="$\\partial \\Pi^e_{rr} / \\partial r$", linewidth=2)
    ax.plot(avg1_rs, avg1(np.divide(prr_e.data - paa_e.data, rs)), label="$(\\Pi^e_{rr} - \\Pi^e_{\\phi \\phi}) / r$", linewidth=2)
    ax.plot(avg2_rs, avg2(n_e.data * er.data), label="$n_e E_r$", linewidth=2)
    ax.plot(rs, ja_e.data * bz.data, label="$J^e_{\\phi} B_z$", linewidth=2)
    ax.plot(avg2_rs, avg2(-(np.gradient(prr_e.data, dx) + np.divide(prr_e.data - paa_e.data, rs)) - n_e.data * er.data + ja_e.data * bz.data), label="res. force", linewidth=2, linestyle="--", color="red")

    annotate_x(forces_i.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.12, size=smol)

    v = (-0.002, +0.001)
    forces_i.draw_info(xlim=(0, 30), ylim=v, yticks=np.linspace(*v, 7))
    forces_i.axes_position.legend(fontsize=ssmol * 0.8, loc="lower right")
    forces_i.axes_position.grid()
    forces_i.axes_position.set_title("Forces on ions", fontsize=smol)

    forces_e.draw_info(xlim=(0, 30), ylim=v, yticks=np.linspace(*v, 7))
    forces_e.axes_position.legend(fontsize=ssmol * 0.8, loc="lower right")
    forces_e.axes_position.grid()
    forces_e.axes_position.set_title("Forces on electrons", fontsize=smol)

    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

    for diag in [forces_i, forces_e]:
        diag.clear()

if __name__ == "__main__":
    ncols=1
    nrows=2

    fig = plt.figure(figsize=(8 * ncols * 1.5, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    forces_i = Field(prefix, subplot(fig, gs, 0, 0))
    forces_e = Field(prefix, subplot(fig, gs, 0, 1))

    n_i = particles_field("Ions", "Dens", "Z")
    ja_i = particles_field("Ions", "Current", "Z")
    prr_i = particles_field("Ions", pressures["Prr"], "Z")
    paa_i = particles_field("Ions", pressures["Paa"], "Z")

    n_e = particles_field("Electrons", "Dens", "Z")
    ja_e = particles_field("Electrons", "Current", "Z")
    prr_e = particles_field("Electrons", pressures["Prr"], "Z")
    paa_e = particles_field("Electrons", pressures["Paa"], "Z")

    er = electric_field("Z")
    bz = magnetic_field("Z")

    res_dir = f"{params_path}/Forces"
    mkdir(res_dir)

    offset = 5

    t0 = 10 # 800 * 5
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range:
        plot_forces(t)
