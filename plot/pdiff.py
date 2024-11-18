#!/usr/bin/env python3

from lib_common import *

ncols=2
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

dPrr_i = Field(None, subplot(fig, gs, 0, 0), boundaries, signed_cmap, None)
dPrr_e = Field(None, subplot(fig, gs, 1, 0), boundaries, signed_cmap, None)

diag_2d = [dPrr_i, dPrr_e]

bx = boundaries[0] + buff * dx
ex = boundaries[1] - buff * dx
by = boundaries[2] + buff * dy
ey = boundaries[3] - buff * dy

arg_2d = {
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xlabel": "$x,~c/\\omega_{pi}$",
    "ylabel": "$y,~c/\\omega_{pi}$",
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

dPrr_i.set_axes_args(title="$\\Pi_{rr}^i - (B_v^2 - B_z^2) / 2$", **arg_2d)
dPrr_e.set_axes_args(title="$\\Pi_{rr}^e - (B_v^2 - B_z^2) / 2$", **arg_2d)

res_dir = f"./{params_path}/Other"
mkdir(res_dir)

t = int(8 * tau / dts)

prr_i = parse_file(get_particles_file("Ions", f"{pressures['Prr']}PlaneAvgZ", t))
prr_e = parse_file(get_particles_file("Electrons", f"{pressures['Prr']}PlaneAvgZ", t))
Bz = parse_file(get_fields_file(t), fields.index("Bz"))

pB = (B0 * B0 - Bz * Bz) / 2
dPrr_i.data = prr_i - pB
dPrr_e.data = prr_e - pB

# Pmax = np.max(dPrr_i.data + dPrr_e.data)
# dPrr_i.data /= Pmax
# dPrr_e.data /= Pmax

for diag in diag_2d:
    diag.axes_position.set_aspect(1)
    diag.draw(add_cbar=True)
    diag.draw_info()

fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", x=0.515, size=big, bbox=bbox)

fig.tight_layout()
fig.tight_layout()

filename = f"{res_dir}/pressures_diff"
print(filename)

fig.savefig(f"{filename}.png")
