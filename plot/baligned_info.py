#!/usr/bin/env python3

from plot import *

def plot_baligned_info(t):
    for tt in range(t, t + offset):
        update_data(tt)

    for diag in [ep, ni, vpi, vti, prri, pzzi, ne, vpe, vte, prre, pzze]:
        diag.data /= offset

    xc = data_shape["Y"][0] // 2
    zc = data_shape["Y"][1] // 2
    xl = xc + 1
    w1 = -1 # 4
    w2 = +1
    xs1 = select_magnetic_line(xl + w1)
    xs2 = select_magnetic_line(xl + w2)
    zs = np.arange(0, data_shape["Y"][1])

    for diag in [ep, b]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

        ax = diag.axes_position
        ax.plot(+(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(+(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs1 - xc) * dx, zs * dz, color="black", linewidth=2)
        ax.plot(-(xs2 - xc) * dx, zs * dz, color="black", linewidth=2)

    zs = np.arange(0, data_shape["Y"][1])
    e_l = average_over_magnetic_tube(ep.data, zs, xs1, xs2, xc) # / T_e

    ni_l = average_over_magnetic_tube(ni.data, zs, xs1, xs2, xc)
    vpi_l = average_over_magnetic_tube(vpi.data, zs, xs1, xs2, xc) # / np.sqrt(T_e / mi_me)
    vti_l = average_over_magnetic_tube(vti.data, zs, xs1, xs2, xc) # / np.sqrt(T_e / mi_me)
    prri_l = average_over_magnetic_tube(prri.data, zs, xs1, xs2, xc) # / T_e
    pzzi_l = average_over_magnetic_tube(pzzi.data, zs, xs1, xs2, xc) # / T_e

    ne_l = -average_over_magnetic_tube(ne.data, zs, xs1, xs2, xc)
    vpe_l = -average_over_magnetic_tube(vpe.data, zs, xs1, xs2, xc) # / np.sqrt(T_e)
    vte_l = average_over_magnetic_tube(vte.data, zs, xs1, xs2, xc) # / np.sqrt(T_e)
    prre_l = average_over_magnetic_tube(prre.data, zs, xs1, xs2, xc) # / T_e
    pzze_l = average_over_magnetic_tube(pzze.data, zs, xs1, xs2, xc) # / T_e

    n_zero_tolerance = 1e-3
    vpi_l = np.divide(vpi_l, ni_l, where=(ni_l > n_zero_tolerance), out=np.zeros_like(vpi_l))
    vti_l = np.divide(vti_l, ni_l, where=(ni_l > n_zero_tolerance), out=np.zeros_like(vti_l))
    prri_l = np.divide(prri_l, ni_l, where=(ni_l > n_zero_tolerance), out=np.zeros_like(prri_l))
    pzzi_l = np.divide(pzzi_l, ni_l, where=(ni_l > n_zero_tolerance), out=np.zeros_like(pzzi_l))

    vpe_l = np.divide(vpe_l, ne_l, where=(ne_l > n_zero_tolerance), out=np.zeros_like(vpe_l))
    vte_l = np.divide(vte_l, ne_l, where=(ne_l > n_zero_tolerance), out=np.zeros_like(vte_l))
    prre_l = np.divide(prre_l, ne_l, where=(ne_l > n_zero_tolerance), out=np.zeros_like(prre_l))
    pzze_l = np.divide(pzze_l, ne_l, where=(ne_l > n_zero_tolerance), out=np.zeros_like(pzze_l))

    ax = phi.axes_position
    ax.plot(zs * dz, e_l, label="$E_{\\|}(z) / T_e$", linewidth=3)

    cos_map = np.ones_like(b.data[:,xl])
    for z, x in enumerate(0.5 * (xs1 + xs2)):
        cos_map[z] = np.divide(bz[z,int(x)], b.data[z,int(x)], where=(b.data[z,int(x)] > 1e-3))

    off = 5
    phi_l = -cumulative_trapezoid(e_l, zs * dz / cos_map, initial=0)

    _t = np.mean(pzze_l[zc-10:zc+10])
    print("prri:", np.mean(prri_l[zc-10:zc+10]),
          "pzzi:", np.mean(pzzi_l[zc-10:zc+10]),
          "prre:", np.mean(prre_l[zc-10:zc+10]),
          "pzze:", np.mean(pzze_l[zc-10:zc+10]))

    phi_l_th = (np.log(ne_l) * _t)[off:-off] / (T_e * 511)
    phi_l_th -= phi_l_th[0]

    ax.plot(zs * dz, phi_l, label="$\\varphi(z) / T_e$", linewidth=3)
    ax.plot((zs * dz)[off:-off], phi_l_th + np.mean(phi_l[:10]), label="$\\ln(n_e(z))$", linewidth=3)

    ni.axes_position.plot(zs * dz, ne_l, linewidth=2, label="electrons")
    ni.axes_position.plot(zs * dz, ni_l, linewidth=2, label="ions")
    ni.axes_position.plot(zs * dz, (ni_l - ne_l), linewidth=2, label="diff", color="magenta")

    vpi.axes_position.plot(zs * dz, vpe_l, zs * dz, vpi_l, linewidth=2)
    vti.axes_position.plot(zs * dz, vte_l, zs * dz, vti_l, linewidth=2)
    prri.axes_position.plot(zs * dz, prre_l, zs * dz, prri_l, linewidth=2)
    pzzi.axes_position.plot(zs * dz, pzze_l, zs * dz, pzzi_l, linewidth=2)

    phi.axes_position.legend(fontsize=0.8 * ssmol, loc="upper left")
    ni.axes_position.legend(fontsize=0.8 * ssmol, loc="upper left")

    for diag in [phi, ni, vpi, vti, prri, pzzi]:
        diag.draw_info()

    # phi.axes_position.set_xlim(0, 50)
    # n.axes_position.set_xlim(0, 50)
    print(np.max(phi_l))

    # dump("E_p", "me c wpe / e", int(t * dts / dt), e_l)
    # dump("Phi", "me c^2 / e", int(t * dts / dt), phi_l)
    # dump("IonsDensity", "n0", int(t * dts / dt), ni_l)
    # dump("IonsVelocity_pe", "c", int(t * dts / dt), vpi_l)
    # dump("IonsVelocity_te", "c", int(t * dts / dt), vti_l)
    # dump("IonsTemperature_rr", "n0 me c^2", int(t * dts / dt), prri_l)
    # dump("IonsTemperature_zz", "n0 me c^2", int(t * dts / dt), pzzi_l)
    # dump("ElectronsDensity", "n0", int(t * dts / dt), ne_l)
    # dump("ElectronsVelocity_pe", "c", int(t * dts / dt), vpe_l)
    # dump("ElectronsVelocity_te", "c", int(t * dts / dt), vte_l)
    # dump("ElectronsTemperature_rr", "n0 me c^2", int(t * dts / dt), prre_l)
    # dump("ElectronsTemperature_zz", "n0 me c^2", int(t * dts / dt), pzze_l)

    annotate_x(ep.axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)
    print(f"{int(t * dts / dt)} [dt] {int(t)} [dts] {t * dts / tau:.3f} [tau]")

    fig.tight_layout()
    fig.savefig(f"{res_dir}/baligned_info_t{str(t).zfill(4)}dts_xl{xl-xc}_w{w2-w1}.png")

def update_data(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    _er = get_parsed_field(ep.path_to_file, "E", "Y", "x", t)
    _ez = get_parsed_field(ep.path_to_file, "E", "Y", "z", t)
    ep.data = agg(ep.data, align(_er, _ez))

    # ions
    ni.data = agg(ni.data, get_parsed_scalar(ni, t))
    prri.data = agg(prri.data, get_parsed_scalar(prri, t))
    pzzi.data = agg(pzzi.data, get_parsed_scalar(pzzi, t))

    _vri = get_parsed_field(vpi.path_to_file, "E", "Y", "x", t)
    _vzi = get_parsed_field(vti.path_to_file, "E", "Y", "z", t)
    _vpi = align(_vri, _vzi)
    vpi.data = agg(vpi.data, _vpi)

    _vtri = _vri - _vpi * np.divide(br, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    _vtzi = _vzi - _vpi * np.divide(bz, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    vti.data = agg(vti.data, np.hypot(_vtri, _vtzi))

    # electrons
    ne.data = agg(ne.data, get_parsed_scalar(ne, t))
    prre.data = agg(prre.data, get_parsed_scalar(prre, t))
    pzze.data = agg(pzze.data, get_parsed_scalar(pzze, t))

    _vre = get_parsed_field(vpe.path_to_file, "E", "Y", "x", t)
    _vze = get_parsed_field(vte.path_to_file, "E", "Y", "z", t)
    _vpe = align(_vre, _vze)
    vpe.data = agg(vpe.data, _vpe)

    _vtre = _vre - _vpe * np.divide(br, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    _vtze = _vze - _vpe * np.divide(bz, b.data, where=(b.data > 1e-3), out=np.zeros_like(b.data))
    vte.data = agg(vte.data, np.hypot(_vtre, _vtze))


if __name__ == "__main__":
    ncols=3
    nrows=3

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1, 1.5, 1.5], height_ratios=[1] * nrows, figure=fig)

    ep = electric_field("Y", subplot(fig, gs, 0, 0), "$E_{\\|}$")
    ep.vmin_vmax = (-0.6e-3, +0.6e-3)

    phi = electric_field("Y", subplot(fig, gs, 1, 0))
    phi.set_axes_args(title="$\\varphi(z) = -\\int_0^z E_{\\|}(\\zeta) d\\zeta$")

    ni = particles_field("Ions", "Dens", "Y", subplot(fig, gs, 2, 0))
    vpi = particles_field("Ions", "Current", "Y", subplot(fig, gs, 1, 1))
    vti = particles_field("Ions", "Current", "Y", subplot(fig, gs, 2, 1))
    prri = particles_field("Ions", "Pxx", "Y", subplot(fig, gs, 2, 2))
    pzzi = particles_field("Ions", "Pzz", "Y", subplot(fig, gs, 1, 2))

    ne = particles_field("Electrons", "Dens", "Y")
    vpe = particles_field("Electrons", "Current", "Y")
    vte = particles_field("Electrons", "Current", "Y")
    prre = particles_field("Electrons", "Pxx", "Y")
    pzze = particles_field("Electrons", "Pzz", "Y")

    ni.set_axes_args(title="$n / n_0$")
    vpi.set_axes_args(title="$v_{\\|} / v_{T}$")
    vti.set_axes_args(title="$v_{\\perp} / v_{T}$")
    prri.set_axes_args(title="$\\Pi_{rr} / (n T)$")
    pzzi.set_axes_args(title="$\\Pi_{zz} / (n T)$")

    b = magnetic_field("Y", subplot(fig, gs, 0, 1), "$|B|$")
    br = get_parsed_field(b.path_to_file, "B", "Y", "x", 0)
    bz = get_parsed_field(b.path_to_file, "B", "Y", "z", 0)
    b.data = np.hypot(br, bz)

    res_dir = f"{params_path}/Other"
    mkdir(res_dir)

    offset = 100
    ts = [
        int(2.0 * tau / dts),
        int(3.0 * tau / dts),
        int(3.5 * tau / dts),
        int(4.0 * tau / dts),
    ]

    for t in ts:
        plot_baligned_info(t)
