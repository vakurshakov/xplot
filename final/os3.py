#!/usr/bin/env python3

from final import *

set_big(big*1.08)
set_smol(smol*1.08)
set_ssmol(ssmol*1.08)

ncols=3
nrows=1

fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1.2, 1, 1], height_ratios=[1] * nrows, figure=fig)

r_max = data_shape["Z"][0] // 2
z_max = data_shape["Y"][1]
rc = r_max
rw = 2

rs = (np.arange(0, r_max) + 0.5) * dx
zs = np.arange(0, z_max) * dz

fmap = (-0.0002, +0.0002)
pmap = (-0.0005, +0.01)

forces     = Field("", subplot(fig, gs, 0, 0))
pressure_i = Field("", subplot(fig, gs, 1, 0))
pressure_e = Field("", subplot(fig, gs, 2, 0))

forces.set_axes_args(title="\\rm Longitudinal forces",            xlim=(0, 200), xticks=np.linspace(0, 200, 6), ylim=fmap, yticks=np.linspace(*fmap, 5))
pressure_i.set_axes_args(title="\\rm Radial pressure, ions",      xlim=(0, 30),  xticks=np.linspace(0, 30, 6),  ylim=pmap, yticks=np.linspace(0, pmap[1], 5))
pressure_e.set_axes_args(title="\\rm Radial pressure, electrons", xlim=(0, 30),  xticks=np.linspace(0, 30, 6),  ylim=pmap, yticks=np.linspace(0, pmap[1], 5))

def create_particles_fields(s):
    return \
        particles_field(s, "Density", "Y"), \
        particles_field(s, "Prr",     "Y"), \
        particles_field(s, "Pzz",     "Y"), \
        particles_field(s, "Density", "Z"), \
        particles_field(s, "Current", "Z"), \
        particles_field(s, "Prr",     "Z"), \
        particles_field(s, "Ppp",     "Z"), \
        particles_field(s, "Pzz",     "Z")

ny_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i, pzz_i = create_particles_fields("Ions")
ny_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e, pzz_e = create_particles_fields("Electrons")

erz = electric_field("Z")
bzy = magnetic_field("Y")
bzz = magnetic_field("Z")

TAU = 4
OFF = 20
TIME = int(TAU * tau / dts) # int(40000 / dts)

def line_mean(d):
    return np.mean(d[:,rc-rw:rc+rw], axis=1)

def parse_particles_data(t, ny, pry, pzy, nz, jaz, prz, paz, pzz):
    ny.data =  agg(ny.data,  line_mean(get_parsed_scalar(ny, t)))
    pry.data = agg(pry.data, line_mean(get_parsed_scalar(pry, t)))
    pzy.data = agg(pzy.data, line_mean(get_parsed_scalar(pzy, t)))

    nz.data  = agg(nz.data,  phi_averaged(get_parsed_scalar(nz,  t), R_MAP))
    prz.data = agg(prz.data, phi_averaged(get_parsed_scalar(prz, t), R_MAP))
    paz.data = agg(paz.data, phi_averaged(get_parsed_scalar(paz, t), R_MAP))
    pzz.data = agg(pzz.data, phi_averaged(get_parsed_scalar(pzz, t), R_MAP))
    jaz.data = agg(jaz.data, phi_averaged(get_parsed_field(jaz, "E", "Z", "x", t)[1], R_MAP))

for t in range(TIME - OFF//2, TIME + OFF//2):
    parse_particles_data(t, ny_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i, pzz_i)
    parse_particles_data(t, ny_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e, pzz_e)
    erz.data = agg(erz.data, phi_averaged(get_parsed_field(erz, "E", "Z", "x", t)[0], R_MAP))

bzy.data = line_mean(get_parsed_field(bzy, "B", "Y", "z", TIME))
bzz.data = phi_averaged(get_parsed_field(bzz, "B", "Z", "z", TIME), R_MAP)
bzz0     = phi_averaged(get_parsed_field(bzz, "B", "Z", "z", 0),    R_MAP)

for d in [ ny_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i, pzz_i, \
           ny_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e, pzz_e, erz ]:
    d.data /= OFF

def avg1(d):
    return sliding_average(d, 25)

pz = avg1(pzy_i.data + pzy_e.data)
dp = avg1(pry_i.data + pry_e.data) - pz
bz = avg1(bzy.data)
dpzz = np.gradient(pz, dz)
dpbz = dp * np.gradient(np.log(bz, where=(np.abs(bz) > 1e-3), out=np.zeros_like(bz)), dz)

ax = forces.axes_position
ax.plot(avg1(zs), -dpzz, label="$-\\partial_z \\Pi_{zz}$",                 linewidth=3)
ax.plot(avg1(zs), +dpbz, label="$(\\Pi_{rr} - \\Pi_{zz})~\\partial_z ln(B)$", linewidth=3)

d1b = np.gradient(bzy.data, dz)
d2b = np.gradient(d1b, dz)

zc = z_max // 2
zw = 2
b = np.mean(bzy.data[zc-zw:zc+zw])
d1b = np.mean(d1b[zc-zw:zc+zw])
d2b = np.mean(d2b[zc-zw:zc+zw])
R_frac = (1 / (2 * b)) * (d2b - (3 / 2) * np.square(d1b) / b)

def draw_pressure(d, nz, jaz, prz, paz, pzz):
    es = ((+1) if d == pressure_i else (-1))

    def integrate(data):
        data = cumulative_trapezoid(data, dx=dx, initial=0)
        data = -(data[-1] - data)
        return data

    dpra = integrate((prz.data - paz.data) / (rs + 0.1))
    dprz = integrate((prz.data - pzz.data) * R_frac * rs)
    pe   = integrate(nz.data * erz.data)
    pb   = integrate(jaz.data * bzz.data)

    ax = d.axes_position
    ax.plot(rs, prz.data,            label="$\\Pi_{rr}$",                          linewidth=4, color="C0")
    ax.plot(rs, dpra,                label="$\\Delta \\Pi_{r}$",                   linewidth=3, color="C1")
    ax.plot(rs, dprz,                label="$\\Delta \\Pi_{R}$",                   linewidth=3, color="C2")
    ax.plot(rs, pe * es,             label="$\\int_{\\infty}^r n E_r dr$",         linewidth=3, color="C3")
    ax.plot(rs, pb,                  label="$\\int_{\\infty}^r J_{\\phi} B_z dr$", linewidth=3, color="C4")
    ax.plot(rs, -dpra -dprz +pb +pe, label="\\rm control",                         linewidth=4, linestyle="--", color="red")

draw_pressure(pressure_i, nz_i, jaz_i, prz_i, paz_i, pzz_i)
draw_pressure(pressure_e, nz_e, jaz_e, prz_e, paz_e, pzz_e)

def get_loc(d):
    if d == forces:
        return dict(bbox_to_anchor=(0.46, 0.26))
    return dict(loc="upper right")

for d in [ forces, pressure_i, pressure_e ]:
    d.draw_info()
    ax = d.axes_position
    ax.legend(fontsize=Fonts.smol*0.85, framealpha=1, **get_loc(d)) #, loc=get_loc(d))
    ax.grid(alpha=0.6)
   
dlmap = [
    (forces, "a"),
    (pressure_i, "b"),
    (pressure_e, "c"),
]

for (d, l) in dlmap:
    annotate_x(d.axes_position, f"\\rm {l}", -0.1, 1.1, size=Fonts.big)

fig.tight_layout(w_pad=-1)
fig.savefig(f"{res_dir}/os3.pdf")
