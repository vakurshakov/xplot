#!/usr/bin/env python3

from final import *

ncols=3
nrows=2

fig = plt.figure(figsize=(8*ncols, 8*nrows))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1]*nrows, figure=fig)

emap = (-2e-3, +2e-3)
vmap = (-3.0, +3.0)

ep = electric_field("Y", subplot(fig, gs, 0, 0), "")
ep.vmin_vmax = emap
ep.axes_args["title"] = "$E_{\\|}(x, z)$"

vp = electric_field("Y", subplot(fig, gs, 0, 1), "")
vp.vmin_vmax = vmap 
vp.axes_args["title"] = "$v_{\\|}^i(x, z) / c_s$"

bx = -28 * dx
ex = +28 * dx
by = boundaries["Y"][2]
ey = boundaries["Y"][3]

def curvilinear(title, i, v):
    f = Field(None, subplot(fig, gs, *i), cmap=signed_cmap)
    f.vmin_vmax=v
    f.boundaries = (bx, ex, by, ey)
    f.axes_args["title"] = title
    f.axes_args["xlim"] = (bx, ex)
    f.axes_args["ylim"] = (by, ey)
    f.axes_args["xlabel"] = "$\\xi,~c/\\omega_{pe}$"
    f.axes_args["ylabel"] = "$z,~c/\\omega_{pe}$"
    f.axes_args["xticks"] = np.linspace(bx, ex, 5)
    f.axes_args["yticks"] = np.linspace(by, ey, 5)
    return f

phi = curvilinear("$\\varphi(\\xi, z)$", (1, 0), (-3, 3))
vpc = curvilinear("$v_{\\|}^i(\\xi, z) / c_s$", (1, 1), vmap)

def linear(title, i, v, nv):
    f = electric_field("Y")
    f.axes_position = subplot(fig, gs, *i)
    f.axes_args["title"] = title
    f.axes_args["xlabel"] = "$z,~c/\\omega_{pe}$"
    f.axes_args["xlim"] = (by, ey)
    f.axes_args["xticks"] = np.linspace(by, ey, 5)
    f.axes_args["ylim"] = v
    f.axes_args["yticks"] = np.linspace(*v, nv)
    return f

phil = linear("$\\varphi(z)$",        (2, 0), (0, 4), 5)
vpcl = linear("$v_{\\|}^i(z) / c_s$", (2, 1), vmap,   7)

b = magnetic_field("Y")
br = get_parsed_field(b, "B", "Y", "x", 0)
bz = get_parsed_field(b, "B", "Y", "z", 0)
b = np.hypot(br, bz)

TAU = 4
OFF = 20
TIME = int(TAU * tau / dts)

for t in range(TIME - OFF//2, TIME + OFF//2):
    er = get_parsed_field(ep, "E", "Y", "x", t)
    ez = get_parsed_field(ep, "E", "Y", "z", t)
    dot = (er * br + ez * bz)
    ep.data = agg(ep.data, np.divide(dot, b, where=(b > 1e-3), out=np.zeros_like(b)))

    niy = get_parsed_scalar(particles_field("Ions", "Density", "Y"), t)
    jiry = get_parsed_field(particles_field("Ions", "Current", "Y"), "E", "Y", "x", t)
    jizy = get_parsed_field(particles_field("Ions", "Current", "Y"), "E", "Y", "z", t)

    d = b * niy
    dot = (jiry * br + jizy * bz)
    vp.data = agg(vp.data, np.divide(dot, d, where=(d > 1e-5), out=np.zeros_like(d)))

cs = np.sqrt(T_e / mi_me)

ep.data /= OFF
vp.data /= OFF * cs

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

def calc_vpc(r0):
    return calc_avg(vp.data, r0)

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
vpc.data = calc_map(calc_vpc)

for d in [ ep, vp, phi, vpc ]:
    d.axes_position.set_aspect(1)
    d.draw(add_cbar=True)
    d.draw_info()

z_th = []
phil_th = []
vpcl_th = []

with open(f"{params_path}/Final/V_and_Fi.dat") as f:
    for l in f.readlines():
        v = l.split("\t")
        z_th.append(float(v[0]) * dz)
        phil_th.append(float(v[2]) / 2)
        vpcl_th.append(float(v[1]))

def draw_linear(p, df, d_th):
    ax = p.axes_position
    dmin = df(xl_min)
    dmax = df(xl_max)
    ax.plot(zs * dz, (dmin[0] + dmin[1]) / 2, label="$|\\xi| = 0$",  linewidth=3)
    ax.plot(zs * dz, (dmax[0] + dmax[1]) / 2, label="$|\\xi| = 12$", linewidth=3)
    ax.plot(z_th,    d_th,                    label="theory",        linewidth=3)
    ax.grid(alpha=0.6)
    ax.legend(loc="upper left", fontsize=ssmol*0.8)
    p.draw_info()

draw_linear(phil, calc_phi, phil_th)
draw_linear(vpcl, calc_vpc, vpcl_th)

zmin = (zs * dz)[0]
zmax = (zs * dz)[-1]
zl = [zmin, zmax]

def line(xl, c, w):
    sw = 10
    xsl = sliding_average(select_magnetic_line(bz, xc + xl + w), sw)
    zsa = sliding_average(zs * dz, sw)
    for ax in [ ep.axes_position, vp.axes_position ]:
        ax.plot(-(xsl - xc) * dx - 0.4, zsa, color="black", linewidth=3)
        ax.plot(+(xsl - xc) * dx + 0.4, zsa, color="black", linewidth=3)
        ax.plot(-(xsl - xc) * dx - 0.2, zsa, color=c,       linewidth=3)
        ax.plot(+(xsl - xc) * dx + 0.2, zsa, color=c,       linewidth=3)
    for ax in [ phi.axes_position, vpc.axes_position ]:
        ax.plot([-(xl + w) * dx + 0.15]*2, zl, color="black", linewidth=3)
        ax.plot([+(xl + w) * dx - 0.15]*2, zl, color="black", linewidth=3)
        ax.plot([-(xl + w) * dx + 0.10]*2, zl, color=c,       linewidth=3)
        ax.plot([+(xl + w) * dx - 0.10]*2, zl, color=c,       linewidth=3)

line(xl_min, "C0", w2)
line(xl_max, "C1", w1)
line(xl_max, "C1", w2)

fig.tight_layout()
fig.savefig(f"{res_dir}/os2.pdf")
