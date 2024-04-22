#!/usr/bin/python3

from final_common import *
from matplotlib.ticker import MultipleLocator

ncols=1
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.5, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"../Final"
mkdir(res_dir)

plot = Field(None, subplot(fig, gs, 0, 0), boundaries)
plot.set_axes_args(title="{\\rm Scheme}", **arg_2d)

ax = plot.axes_position

def plot_circle(ax, r0, xc, yc, **kwargs):
  phi = np.linspace(0, 2 * np.pi)
  ax.plot(xc + r0 * np.sin(phi), yc + r0 * np.cos(phi), **kwargs)


#incjection region
r0 = 15
plot_circle(ax, r0, 0, 0, linestyle='--', c='black', zorder=0)
ax.xaxis.set_minor_locator(MultipleLocator(15))
ax.yaxis.set_minor_locator(MultipleLocator(15))


# particles number
pp = 25

def gen(): return np.random.rand(pp)

# coords
pr = r0 * np.sqrt(gen())
pphi = 2 * np.pi * gen()

xs_e = pr * np.sin(pphi)
ys_e = pr * np.cos(pphi)

xs_i = xs_e + 1
ys_i = ys_e + 1

ax.scatter(xs_e, ys_e, s=90, zorder=1)
ax.scatter(xs_i, ys_i, s=90, zorder=2)

# velocities
def th(s): return np.sqrt(-2.0 * (s / 511) * np.log(gen()))

vx_e = th(T_e) * np.sin(2 * np.pi * gen())
vy_e = th(T_e) * np.sin(2 * np.pi * gen())
vx_i = th(T_i * mi_me) * np.sin(2 * np.pi * gen())
vy_i = th(T_i * mi_me) * np.sin(2 * np.pi * gen())

colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

arg_quiver = {
  "width": 0.003,
  "headwidth": 3,
  "headlength": 5,
}

pp_v = 4
ax.quiver(xs_e[:pp_v], ys_e[:pp_v], vx_e[:pp_v], vy_e[pp_v], color=colors[0], **arg_quiver, zorder=1)
ax.quiver(xs_i[:pp_v], ys_i[:pp_v], vx_i[:pp_v], vy_i[pp_v], color=colors[1], **arg_quiver, zorder=2)

# magnetic field
xc = +45
yc = -45
ax.scatter(xc, yc, s=80, c='r')
plot_circle(ax, 4, xc+0.2, yc-0.1, linestyle='-', c='r')

# parameters
x0 = 1.3
annotate_x(
  ax,
  "$\\tau = 6000~\\omega_{pe}^{-1}$\n"
  "$T_i = 10~{\\rm keV}$\n"
  "$T_e = 1~{\\rm keV}$\n"
  "$m_i / m_e = 100$",
  x=x0,
  y=0.7,
  size=ssmol
)
annotate_x(
  ax,
  "$\\Delta t = 1.5~\\omega_{pe}^{-1} \\approx 0.3~\\Omega_e^{-1}$\n"
  "$\\Delta x = 0.5~c/\\omega_{pe} \\approx 2~\\rho_e$",
  x=x0,
  y=0.4,
  size=ssmol,
)
annotate_x(
  ax,
  "$\\mathbf{B} = 0.2\,\\mathbf{e}_z$",
  x=x0,
  y=0.1,
  size=ssmol
)
ax.annotate("$R = 15$", (r0 * np.cos(3 * np.pi / 4) - 0.2, r0 * np.sin(3 * np.pi / 4) + 0.2), (-50, 30), **arg_anno)

ax.set_aspect(1)
ax.grid(alpha=0.5, zorder=0)

plot.draw_info()

fig.tight_layout()

fig.savefig(f"{res_dir}/scheme.pdf")
fig.clear()
