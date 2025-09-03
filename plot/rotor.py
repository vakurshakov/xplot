#!/usr/bin/env python3

from plot import *

def plot_currents(t):
    filename = f"{res_dir}/{str(t // offset).zfill(4)}.png"
    if not timestep_should_be_processed(t, filename, False):
        return

    for j in jz:
        j.data = parse_file(get_parsed_file(t, j.path_to_file), 2)
    
    for i, j in enumerate(jz[:3]):
        jz[i].data += jz[i+3].data

    bxx = parse_file(get_parsed_file(t, get_fields_path('B', 'X')), 0) - B0X
    byy = parse_file(get_parsed_file(t, get_fields_path('B', 'Y')), 1) - B0Y
    bxz = parse_file(get_parsed_file(t, get_fields_path('B', 'Z')), 0)
    byz = parse_file(get_parsed_file(t, get_fields_path('B', 'Z')), 1)


    b[0].data = +np.gradient(bxx, dy, axis=1)
    b[1].data = -np.gradient(byy, dx, axis=1)
    b[2].data = +np.gradient(bxz, dy, axis=1) - np.gradient(byz, dx, axis=0)

    for diag in [ *jz[:3], *b ]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", y=0.99, bbox=bbox, fontsize=big)
    fig.tight_layout(rect=(0, 0, 1, 0.99))
    fig.savefig(filename)

    for diag in [ *jz[:3] ]:
        diag.clear()

B0X = magnetic_field("X")
B0X = parse_file(get_parsed_file(0, B0X.path_to_file), 0)

B0Y = magnetic_field("Y")
B0Y = parse_file(get_parsed_file(0, B0Y.path_to_file), 1)

if __name__ == "__main__":
    ncols=3
    nrows=2

    fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1.1))
    gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

    v = 0.02
    jz = [None]*6
    jz[0] = particles_field("Ions", "Current", "X", subplot(fig, gs, 0, 0), "\\rm sum, $J_z$", (-v, +v))
    jz[1] = particles_field("Ions", "Current", "Y", subplot(fig, gs, 1, 0), "\\rm sum, $J_z$", (-v, +v))
    jz[2] = particles_field("Ions", "Current", "Z", subplot(fig, gs, 2, 0), "\\rm sum, $J_z$", (-v, +v))

    jz[3] = particles_field("Electrons", "Current", "X")
    jz[4] = particles_field("Electrons", "Current", "Y")
    jz[5] = particles_field("Electrons", "Current", "Z")

    vmap = (-2e-2, +2e-2)
    b = [None]*3
    b[0] = magnetic_field("X", subplot(fig, gs, 0, 1), "$(\\partial_y B_x)$", vmap)
    b[1] = magnetic_field("Y", subplot(fig, gs, 1, 1), "$(\\partial_x B_y)$", vmap)
    b[2] = magnetic_field("Z", subplot(fig, gs, 2, 1), "$($\\rm rot$\\,B)_z$", vmap)

    res_dir = f"{params_path}/Rotor"
    mkdir(res_dir)

    offset = 25
    t0 = 0
    t_range = create_t_range(t0, int(time / dts), offset)

    for t in t_range[::-1]:
        plot_currents(t)
