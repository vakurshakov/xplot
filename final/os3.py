#!/usr/bin/env python3

from final import *

ncols=2
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.5, 8 * nrows * 1.1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

r_max = data_shape["Z"][0] // 2
z_max = data_shape["Y"][1]
rc = r_max
rw = 2

rs = (np.arange(0, r_max) + 0.5) * dx
zs = np.arange(0, z_max) * dz

fmap = (-0.0002, +0.0002)
pmap = (-0.002, +0.01)

forces_i = Field("",   subplot(fig, gs, 0, 0))
forces_e = Field("",   subplot(fig, gs, 0, 1))
pressure_i = Field("", subplot(fig, gs, 1, 0))
pressure_e = Field("", subplot(fig, gs, 1, 1))

def args(map, x, ny):
    return dict(xlim=(0, x), xticks=np.linspace(0, x, 6), ylim=map, yticks=np.linspace(*map, ny))

forces_i.set_axes_args(**args(fmap, 200, 5),   title = "Longitudinal forces, ions"     )
forces_e.set_axes_args(**args(fmap, 200, 5),   title = "Longitudinal forces, electrons")
pressure_i.set_axes_args(**args(pmap, 30, 7), title = "Radial pressure, ions"         )
pressure_e.set_axes_args(**args(pmap, 30, 7), title = "Radial pressure, electrons"    )

def create_particles_fields(s):
    return \
        particles_field(s, "Density", "Y"), \
        particles_field(s, "Current", "Y"), \
        particles_field(s, "Current", "Y"), \
        particles_field(s, "Prr",     "Y"), \
        particles_field(s, "Pzz",     "Y"), \
        particles_field(s, "Density", "Z"), \
        particles_field(s, "Current", "Z"), \
        particles_field(s, "Prr",     "Z"), \
        particles_field(s, "Ppp",     "Z")

ny_i, jry_i, jzy_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i = create_particles_fields("Ions")
ny_e, jry_e, jzy_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e = create_particles_fields("Electrons")

ezy = electric_field("Y")
erz = electric_field("Z")
bzy = magnetic_field("Y")
bzz = magnetic_field("Z")

TAU = 4
OFF = 20
TIME = int(TAU * tau / dts) # int(40000 / dts)

def line_mean(d):
    return np.mean(d[:,rc-rw:rc+rw], axis=1)

def parse_particles_data(t, ny, jry, jzy, pry, pzy, nz, jaz, prz, paz):
    ny.data =  agg(ny.data,  line_mean(get_parsed_scalar(ny, t)))
    pry.data = agg(pry.data, line_mean(get_parsed_scalar(pry, t)))
    pzy.data = agg(pzy.data, line_mean(get_parsed_scalar(pzy, t)))
    jry.data = agg(jry.data, line_mean(get_parsed_field(jzy, "E", "Y", "x", t)))
    jzy.data = agg(jzy.data, line_mean(get_parsed_field(jzy, "E", "Y", "z", t)))

    nz.data  = agg(nz.data,  phi_averaged(get_parsed_scalar(nz,  t), R_MAP))
    prz.data = agg(prz.data, phi_averaged(get_parsed_scalar(prz, t), R_MAP))
    paz.data = agg(paz.data, phi_averaged(get_parsed_scalar(paz, t), R_MAP))
    jaz.data = agg(jaz.data, phi_averaged(get_parsed_field(jaz, "E", "Z", "x", t)[1], R_MAP))

for t in range(TIME - OFF//2, TIME + OFF//2):
    parse_particles_data(t, ny_i, jry_i, jzy_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i)
    parse_particles_data(t, ny_e, jry_e, jzy_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e)

    ezy.data = agg(ezy.data, line_mean(get_parsed_field(ezy, "E", "Y", "z", t)))
    erz.data = agg(erz.data, phi_averaged(get_parsed_field(erz, "E", "Z", "x", t)[0], R_MAP))

bzy.data = line_mean(get_parsed_field(bzy, "B", "Y", "z", TIME))

bzz.data = phi_averaged(get_parsed_field(bzz, "B", "Z", "z", TIME), R_MAP)
bzz0     = phi_averaged(get_parsed_field(bzz, "B", "Z", "z", 0),    R_MAP)

for d in [ ny_i, jry_i, jzy_i, pry_i, pzy_i, nz_i, jaz_i, prz_i, paz_i, \
           ny_e, jry_e, jzy_e, pry_e, pzy_e, nz_e, jaz_e, prz_e, paz_e, \
           ezy, erz ]:
    d.data /= OFF

def avg1(d):
    return sliding_average(d, 25)

def draw_forces(d, ny, jry, jzy, pry, pzy):
    ms = (mi_me if d == forces_i else 1)
    es = ((+1)  if d == forces_i else (-1))

    # jzy.data = np.divide(jzy.data, ny.data, where=(np.abs(ny.data) > 1e-3), out=np.zeros_like(ny.data))
    # pry.data = pry.data - ms * (ny.data / es) * np.square(jry.data)
    # pzy.data = pzy.data - ms * (ny.data / es) * np.square(jzy.data)
    
    pz = avg1(pzy.data)
    bz = avg1(bzy.data)
    dpzz = np.gradient(pz, dz)
    dpbz = avg1(pry.data - pzy.data) * np.gradient(np.log(bz, where=(np.abs(bz) > 1e-3), out=np.zeros_like(bz)), dz)
    enEz = avg1(ny.data * ezy.data) * es
    
    ax = d.axes_position
    ax.plot(avg1(zs), dpzz, label="$\\partial_z P_{zz}$",                      linewidth=3)
    ax.plot(avg1(zs), dpbz, label="$(P_{rr} - P_{zz})~\\partial_z ln(B)$", linewidth=3)
    ax.plot(avg1(zs), enEz, label="$e n E_z$",                                  linewidth=3)
    ax.plot(avg1(zs), -dpzz -dpbz + enEz, linestyle="--", color="r",                                  linewidth=3)

draw_forces(forces_i, ny_i, jry_i, jzy_i, pry_i, pzy_i)
draw_forces(forces_e, ny_e, jry_e, jzy_e, pry_e, pzy_e)

def draw_pressure(d, nz, jaz, prz, paz):
    es = ((+1) if d == pressure_i else (-1))

    dpr = cumulative_trapezoid((prz.data - paz.data) / (rs + 0.1), dx=dx, initial=0)
    dpr = +(dpr - dpr[-1])

    pe = cumulative_trapezoid(nz.data * erz.data, dx=dx, initial=0)
    pe = +(pe - pe[-1]) * es

    pb = cumulative_trapezoid(jaz.data * bzz.data, dx=dx, initial=0)
    pb = +(pb - pb[-1])

    ax = d.axes_position
    ax.plot(rs, prz.data, label="$\\Pi_{rr}$",               linewidth=3)
    ax.plot(rs, dpr, label="$\\Delta \\Pi_{rr}$",            linewidth=3)
    ax.plot(rs, pe, label="$\\int_{\\infty}^r e n E_r dr$", linewidth=3)
    ax.plot(rs, pb, label="$\\int_{\\infty}^r J_{\\phi} B_z dr$",           linewidth=3)

draw_pressure(pressure_i, nz_i, jaz_i, prz_i, paz_i)
draw_pressure(pressure_e, nz_e, jaz_e, prz_e, paz_e)
    
def get_loc(d):
    if d in [ forces_i, forces_e ]:
        return "lower left"
    return "upper right"

for d in [ forces_i, forces_e, pressure_i, pressure_e ]:
    d.draw_info()
    ax = d.axes_position
    ax.legend(fontsize=ssmol, loc=get_loc(d))
    ax.grid(alpha=0.6)

fig.tight_layout()
fig.savefig(f"{res_dir}/os3.pdf")
