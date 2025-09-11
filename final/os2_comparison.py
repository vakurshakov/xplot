#!/usr/bin/env python3

from final import *

ncols=3
nrows=1

fig = plt.figure(figsize=(8*ncols, 8*nrows))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1]*nrows, figure=fig)

emap = (-2e-3, +2e-3)

ep = electric_field("Y", subplot(fig, gs, 0, 0), "")
ep.vmin_vmax = emap
ep.axes_args["title"] = "$E_{\\|}(x, z)$"

bx = -28 * dx
ex = +28 * dx
by = boundaries["Y"][2]
ey = boundaries["Y"][3]

phi = Field(None, subplot(fig, gs, 1, 0), cmap=signed_cmap, vmin_vmax=(-3, 3))
phi.boundaries = (bx, ex, by, ey)
phi.axes_args["title"] = "$\\varphi(\\xi, z)$"
phi.axes_args["xlim"] = (bx, ex)
phi.axes_args["ylim"] = (by, ey)
phi.axes_args["xlabel"] = "$\\xi,~c/\\omega_{pe}$"
phi.axes_args["ylabel"] = "$z,~c/\\omega_{pe}$"
phi.axes_args["xticks"] = np.linspace(bx, ex, 5)
phi.axes_args["yticks"] = np.linspace(by, ey, 5)

phil = electric_field("Y")
phil.axes_position = subplot(fig, gs, 2, 0)
phil.axes_args["title"] = "$\\varphi(z)$"
phil.axes_args["xlabel"] = "$z,~c/\\omega_{pe}$"
phil.axes_args["xlim"] = (by, ey)
phil.axes_args["xticks"] = np.linspace(by, ey, 5)
phil.axes_args["ylim"] = (0, 4)
phil.axes_args["yticks"] = np.linspace(0, 4, 5)

b = magnetic_field("Y")
br = get_parsed_field(b, "B", "Y", "x", 0)
bz = get_parsed_field(b, "B", "Y", "z", 0)
b = np.hypot(br, bz)

TIME = 12_487 # 1/wpe
# TIME = 22_477 # 1/wpe
TIME = int(TIME / dts)

OFF = 20

for t in range(TIME - OFF//2, TIME + OFF//2):
    er = get_parsed_field(ep, "E", "Y", "x", t)
    ez = get_parsed_field(ep, "E", "Y", "z", t)
    dot = (er * br + ez * bz)
    ep.data = agg(ep.data, np.divide(dot, b, where=(b > 1e-3), out=np.zeros_like(b)))
ep.data /= OFF

xc = data_shape["Y"][0] // 2
zs = np.arange(0, data_shape["Y"][1])

w1 = -2
w2 = +2
xl_min = w2
xl_max = int(ex / dx - 2)

def calc_avg(data, r0):
    xl = xc + r0
    xs1 = select_magnetic_line(bz, xl + w1)
    xs2 = select_magnetic_line(bz, xl + w2)

    e_l1 = np.zeros(data_shape["Y"][1])
    e_l2 = np.zeros(data_shape["Y"][1])
    for z in zs:
        rs = (np.arange(xs1[z], xs2[z] + 1) - xc) * dx
        area = cumulative_trapezoid(2 * np.pi * rs, rs, initial=0)[-1]
        e_l1[z] += cumulative_trapezoid(data[z, xs1[z]:(xs2[z] + 1)] * (2 * np.pi * rs), rs, initial=0)[-1] / area
        e_l2[z] += cumulative_trapezoid(data[z, (2 * xc - (xs2[z] + 1)):(2 * xc - xs1[z])] * (2 * np.pi * rs), rs, initial=0)[-1] / area
    return e_l1, e_l2

def calc_phi(r0):
    d1, d2 = calc_avg(ep.data, r0)
    d1 = -cumulative_trapezoid(d1, zs * dz, initial=0) / T_i
    d2 = -cumulative_trapezoid(d2, zs * dz, initial=0) / T_i
    return d1, d2 

def calc_map(calc):
    data = np.zeros((len(zs), 2 * xl_max))
    d1, d2 = calc(xl_min)

    for i in range(w1, w2):
        data[:, xl_max+i] = (d1 + d2) / 2

    for i in np.arange(w2, xl_max):
        d1, d2 = calc(i)
        data[:, xl_max-i] = d1
        data[:, xl_max+i] = d2

    data[:, 0] = d1
    return data

phi.data = calc_map(calc_phi)

for d in [ ep, phi ]:
    d.axes_position.set_aspect(1)
    d.draw(add_cbar=True)
    d.draw_info()

def draw_linear(p, df):
    dmin = df(xl_min)
    dmin = (dmin[0] + dmin[1]) / 2
    dump("phi", "T_i", TIME, dmin)

    ax = p.axes_position
    ax.plot(zs * dz, dmin, label="$|\\xi| = 0$",  linewidth=3)
    ax.legend(loc="upper left", fontsize=ssmol*0.8)
    ax.grid(alpha=0.6)
    p.draw_info()

draw_linear(phil, calc_phi)

zmin = (zs * dz)[0]
zmax = (zs * dz)[-1]
zl = [zmin, zmax]

def line(xl, c, w):
    sw = 10
    xsl = sliding_average(select_magnetic_line(bz, xc + xl + w), sw)
    zsa = sliding_average(zs * dz, sw)
    ax = ep.axes_position
    ax.plot(-(xsl - xc) * dx - 0.4, zsa, color="black", linewidth=3)
    ax.plot(+(xsl - xc) * dx + 0.4, zsa, color="black", linewidth=3)
    ax.plot(-(xsl - xc) * dx - 0.2, zsa, color=c,       linewidth=3)
    ax.plot(+(xsl - xc) * dx + 0.2, zsa, color=c,       linewidth=3)
    ax = phi.axes_position
    ax.plot([-(xl + w) * dx + 0.15]*2, zl, color="black", linewidth=3)
    ax.plot([+(xl + w) * dx - 0.15]*2, zl, color="black", linewidth=3)
    ax.plot([-(xl + w) * dx + 0.10]*2, zl, color=c,       linewidth=3)
    ax.plot([+(xl + w) * dx - 0.10]*2, zl, color=c,       linewidth=3)

line(xl_min, "C0", w2)

fig.suptitle(f"$\\omega_{{pe}} t = {TIME * dts}$", x=0.515, y=0.975, bbox=bbox, fontsize=big)
fig.tight_layout(rect=(0, 0, 1, 0.98))

fig.savefig(f"{res_dir}/os2_comparison_{TIME}.pdf")
