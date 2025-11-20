#!/usr/bin/env python3

from final import *
import matplotlib.patches as patches

set_big(10)
set_smol(10)
set_ssmol(9)

plt.rc('text', usetex=True)
plt.rc('axes', titlesize=10, labelsize=11)
plt.rc('xtick', labelsize=9.5)
plt.rc('ytick', labelsize=9.5)
plt.rc('legend', fontsize=9.5)
plt.rc('figure', titlesize=10)
plt.rc('lines', linewidth=1.3)

ncols=2
nrows=1

fig = plt.figure(figsize=(5.5, 2.5))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1.5, 1], height_ratios=[1], figure=fig)

Pxz = electric_field("Y", subplot(fig, gs, 0, 0), "")
Pxy = electric_field("Z", subplot(fig, gs, 1, 0), "")

Pxz.set_axes_args(
    xlim=Pxz.axes_args["ylim"],
    xlabel=Pxz.axes_args["ylabel"],
    xticks=Pxz.axes_args["yticks"],

    ylim=(-Nx*dx/2, +Nx*dx/2),
    ylabel=Pxz.axes_args["xlabel"],
    yticks=Pxz.axes_args["xticks"],
)

Pxy.set_axes_args(
    xlim=(-Nx*dx/2, +Nx*dx/2),
    ylim=(-Ny*dy/2, +Ny*dy/2),
)

Pxz.axes_args.pop("title")
Pxy.axes_args.pop("title")


b = magnetic_field("Y")
br = get_parsed_field(b, "B", "Y", "x", 0)
bz = get_parsed_field(b, "B", "Y", "z", 0)
b = np.hypot(br, bz)

xc = data_shape["Y"][0] // 2
zc = Nz*dz/2
zs = np.arange(0, data_shape["Y"][1])

Pxz.draw_info()
Pxy.draw_info()
Pxz.axes_position.set_aspect(1.9)
Pxy.axes_position.set_aspect(1)

zmin = (zs * dz)[0]
zmax = (zs * dz)[-1]
zl = [zmin, zmax]

ex = +28 * dx
w1 = -2
w2 = +2
xl_min = w2 + 3
xl_max = int(ex / dx - 1)
sw = 15

###
ax = Pxz.axes_position
xl = xl_max - 3
c = "C0"
w = w2
xsl = sliding_average(select_magnetic_line(bz, xc + xl + w), sw)
zsa = sliding_average(zs * dz, sw)
xsa = (xsl - xc) * dx
ax.plot(zsa, -xsa - 0.4, color="black", linewidth=1)
ax.plot(zsa, +xsa + 0.4, color="black", linewidth=1)
ax.plot(zsa, -xsa - 0.2, color=c,       linewidth=1)
ax.plot(zsa, +xsa + 0.2, color=c,       linewidth=1)

zm = 40
xsa0 = xsa[int((zc-zm)/dz)]
ax.plot([zc-zm, zc-zm], [-xsa0-1, -xsa0-10], linestyle=":", linewidth=0.8, color="black")
ax.plot([zc+zm, zc+zm], [-xsa0-1, -xsa0-10], linestyle=":", linewidth=0.8, color="black")
ax.arrow(zc, -xsa0-8.5, -zm+3, 0, linewidth=0.8, head_width=0.4, head_length=0.8)
ax.arrow(zc, -xsa0-8.5, +zm-3, 0, linewidth=0.8, head_width=0.4, head_length=0.8)
annotate_x(ax, "$80~c/\\omega_{pe}$", y=0.15, size=Fonts.ssmol*0.7, bbox=None)

xl = xl_min
c = "C1"
w = w2
xsl = sliding_average(select_magnetic_line(bz, xc + xl + w), sw)
zsa = sliding_average(zs * dz, sw)
xsa = (xsl - xc) * dx
ax.plot(zsa, -xsa - 0.4, color="black", linewidth=1)
ax.plot(zsa, +xsa + 0.4, color="black", linewidth=1)
ax.plot(zsa, -xsa - 0.2, color=c,       linewidth=1)
ax.plot(zsa, +xsa + 0.2, color=c,       linewidth=1)

ax.fill_between([zc-rz, zc+rz], [-rr, -rr], [+rr, +rr], color="lightgrey", zorder=10)#, alpha=0.9)
ax.plot([zc-rz, zc+rz, zc+rz, zc-rz, zc-rz], [-rr, -rr, +rr, +rr, -rr], color="black", linewidth=0.9, zorder=10)

rb=30
ax.plot(zl, [+rb, +rb], linestyle="--", linewidth=0.8, color="black")
ax.plot(zl, [-rb, -rb], linestyle="--", linewidth=0.8, color="black")

annotate_x(ax, "$y = 0~c/\\omega_{pe}$", size=Fonts.smol, bbox=None, y=1.05)

###
ax = Pxy.axes_position
phi = np.linspace(0, 2*np.pi, 100)
ax.plot(rr*np.cos(phi), rr*np.sin(phi), color="black", linewidth=0.8, zorder=10)
ax.plot(rb*np.cos(phi), rb*np.sin(phi), color="black", linestyle="--", linewidth=0.8, zorder=10)
circle = patches.Circle((0,0), rr, color="lightgrey", fill=True)
ax.add_patch(circle)

ax.plot([-rr, -rr], [-2, -15], linestyle=":", linewidth=0.8, color="black")
ax.plot([+rr, +rr], [-2, -15], linestyle=":", linewidth=0.8, color="black")
ax.arrow(0, -13.6, -rr+2, 0, linewidth=0.8, head_width=0.4, head_length=0.6)
ax.arrow(0, -13.6, +rr-2, 0, linewidth=0.8, head_width=0.4, head_length=0.6)
annotate_x(ax, "$20~c/\\omega_{pe} \\approx 6.6 \\rho_i$", y=0.20, size=Fonts.ssmol*0.7, bbox=None)

annotate_x(ax, "$z = 100~c/\\omega_{pe}$", size=Fonts.smol, bbox=None, y=1.05)

fs=Fonts.big*0.6
fig.text(0.25, 0.813, "\\rm absorbing layer", fontsize=fs, ha="center")
fig.text(0.25, 0.773, "\\rm vacuum", fontsize=fs, ha="center")
fig.text(0.543, 0.55, "\\rm ideally conducting wall", fontsize=fs, ha="center", va="center", rotation=90)

Pxz.axes_position.annotate(
    "\\rm region of \n\\rm plasma injection",
    xy=(105, 10), 
    xytext=(100, 18),
    arrowprops=dict(arrowstyle="->", linewidth=0.7), #shrink=0.05),
    fontsize=fs,
    multialignment="center",
)
             
fig.text(0.810, 0.813, "\\rm absorbing layer", fontsize=fs, ha="center")
fig.text(0.810, 0.770, "\\rm vacuum", fontsize=fs, ha="center")

Pxy.axes_position.annotate(
    "\\rm region of \n\\rm plasma injection",
    xy=(7, 7), 
    xytext=(10, 15),
    fontsize=fs,
    multialignment="center",
    arrowprops=dict(arrowstyle="->", linewidth=0.7), #shrink=0.05),
    bbox=dict(facecolor='white', edgecolor='white', boxstyle='round,pad=0.25'),
    zorder=30
)

fig.tight_layout()# rect=(-0.02, -0.07, 1, 1))
fig.savefig(f"{res_dir}/os0.pdf")
