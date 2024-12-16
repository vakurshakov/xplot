#!/usr/bin/env python3

from plot import *

### Figure parameters.

fig, gs = figure()

plot = magnetic_field("Y", subplot(fig, gs), "${\\rm Computational~domain,}~\\mathbf{B}_0$")
plot.vmin_vmax = (0, 0.5)

_bx = get_parsed_field(plot.path_to_file, "B", "Y", "x", 0)
_bz = get_parsed_field(plot.path_to_file, "B", "Y", "z", 0)

plot.data = np.sqrt(_bx * _bx + _bz * _bz)
plot.draw(add_cbar=True, cbar_ticks_num=6)
plot.draw_info()

ax = plot.axes_position

xhw = rr
z0 = 0.5 * data_shape["Y"][1] * dz
zhw = rz


ax.plot([-xhw, +xhw, +xhw, -xhw, -xhw], [z0-zhw, z0-zhw, z0+zhw, z0+zhw, z0-zhw], color="orange", alpha=1)
ax.fill_betweenx([z0-zhw, z0+zhw], [-xhw, -xhw], [+xhw, +xhw], color="orange", alpha=0.3)

ndots = 15
pad = 0.5
diff = 0.5
xs = np.random.uniform(-xhw + pad, +xhw - pad, ndots)
ys = np.random.uniform(z0-zhw + pad, z0+zhw - pad, ndots)
ax.scatter(xs - diff, ys - diff, s=100)
ax.scatter(xs, ys, s=100)

# def plot_coil(z, w):
#     rcoil = 60 * dx
#     rs = np.linspace(-rcoil, +rcoil, 50)
#     zp = (z + w) * np.ones_like(rs)
#     zm = (z - w) * np.ones_like(rs)
#     ax.plot(rs, zm, color="grey", alpha=0.4)
#     ax.plot(rs, zp, color="grey", alpha=0.4)
#     ax.fill_between(rs, zm, zp, color="grey", alpha=0.4)
# plot_coil(0, 0.1)
# plot_coil(5, 0.1)
# plot_coil(10, 0.1)


s = data_shape["Y"]
xs = (np.arange(0, s[0]) - 0.5 * s[0]) * dx
zhw = np.arange(0, s[1]) * dz
X, Y = np.meshgrid(xs, zhw)

ax.streamplot(
    X, Y, _bx, _bz,
    density=0.5,
    arrowsize=2,
    linewidth=2,
)

fig.tight_layout()
fig.tight_layout()
fig.subplots_adjust(wspace=0.3)

filename = f"{params_path}/Other/magnetic_field_scheme.png"

print(filename)
fig.savefig(f"{filename}")
