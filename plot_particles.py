#!/usr/bin/python3

from lib_common import *

def particles_field(sort, diag_name, plane, subplot=None, title=None, v=None, cmap=signed_cmap):
    prefix = f"Particles/{sort}/Diag2D/{diag_name}Plane{plane}_{slices[plane][-1]}"
    field = Field(prefix, subplot, None, cmap, v)
    if title != None:
        generate_info(field, plane, title)
    return field

def plot_particles3(s, t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename):
        return

    jrs[0].data = get_parsed_field(jrs[0], "E", "Y", "x", t)
    jas[0].data = get_parsed_field(jas[0], "E", "Y", "y", t)
    jrs[1].data, jas[1].data = get_parsed_field(jrs[1], "E", "Z", "", t)

    for i, plane in enumerate(planes):
        e = (-1) if s == "Electrons" else (+1)
        ns[i].data = e * parse_file(f"{get_prefix(t)}/{ns[i].path_to_file}_{str(t).zfill(4)}")
        jzs[i].data = get_parsed_field(jzs[i], "E", plane, "z", t)

    for i, diag in enumerate(jrs + jas + jzs + ns):
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    annotate_x(jrs[0].axes_position, "$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=1.2)

    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(filename)

    for diag in jrs + jas + ns:
        diag.clear()


if __name__ == "__main__":
    planes = planes[1:]

    for s in sorts:
        ncols=4
        nrows=2

        fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.2))
        gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1.2, 1], figure=fig)

        v = 1
        ns = []
        ns.append(particles_field(s, "Dens", "Y", subplot(fig, gs, 0, 0), f"$n_{s[0].lower()}$", (0, v), unsigned_cmap))
        ns.append(particles_field(s, "Dens", "Z", subplot(fig, gs, 0, 1), f"$n_{s[0].lower()}$", (0, v), unsigned_cmap))

        v = 0.02
        jrs = []
        jrs.append(particles_field(s, "Current", "Y", subplot(fig, gs, 1, 0), "$J_r^"f"{s[0].lower()}""$", (-v, +v)))
        jrs.append(particles_field(s, "Current", "Z", subplot(fig, gs, 1, 1), "$J_r^"f"{s[0].lower()}""$", (-v, +v)))

        jas = []
        jas.append(particles_field(s, "Current", "Y", subplot(fig, gs, 2, 0), "$J_{\\phi}^"f"{s[0].lower()}""$", (-v, +v)))
        jas.append(particles_field(s, "Current", "Z", subplot(fig, gs, 2, 1), "$J_{\\phi}^"f"{s[0].lower()}""$", (-v, +v)))

        jzs = []
        jzs.append(particles_field(s, "Current", "Y", subplot(fig, gs, 3, 0), "$J_z^"f"{s[0].lower()}""$", (-v, +v)))
        jzs.append(particles_field(s, "Current", "Z", subplot(fig, gs, 3, 1), "$J_z^"f"{s[0].lower()}""$", (-v, +v)))

        res_dir = f"{params_path}/Info_{s}"
        mkdir(res_dir)

        offset = 5
        t0 = 0
        t_range = create_t_range(t0, int(time / dts), offset)

        for t in t_range:
            plot_particles3(s, t)
