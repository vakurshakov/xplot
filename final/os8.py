#!/usr/bin/env python3

from final import *
from xplot.fourier.fourier import *

set_big(11)
set_smol(11)
set_ssmol(10)

plt.rc('text', usetex=True)
plt.rc('axes', titlesize=11, labelsize=12)
plt.rc('xtick', labelsize=10.5)
plt.rc('ytick', labelsize=10.5)
plt.rc('legend', fontsize=10.5)
plt.rc('figure', titlesize=11)
plt.rc('lines', linewidth=1.3)

import ConductiveWall_FixRotor.parameters_circle as ParamsCircle
import ConductiveWall_FixRotor.parameters_square as ParamsSquare

ncols=4
nrows=5

fig = plt.figure(figsize=(8, 10))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1, 1, 1, 0.06, 1.7], figure=fig)

nmap = (0, 5)
bzmap = (0, 0.2)
bpmap = (0, 1)

bz = [
    magnetic_field("Z", subplot(fig, gs, 0, 0), "", bzmap, unsigned_cmap), # t = 4, Circle
    magnetic_field("Z", subplot(fig, gs, 0, 1), "", bzmap, unsigned_cmap), # t = 7, Circle
    magnetic_field("Z", subplot(fig, gs, 0, 2), "", bzmap, unsigned_cmap), # t = 10, Circle

    magnetic_field("Z", subplot(fig, gs, 2, 0), "", bzmap, unsigned_cmap), # t = 4, Square
    magnetic_field("Z", subplot(fig, gs, 2, 1), "", bzmap, unsigned_cmap), # t = 7, Square
    magnetic_field("Z", subplot(fig, gs, 2, 2), "", bzmap, unsigned_cmap), # t = 10, Square
]

ns = [
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 1, 0), "", nmap, unsigned_cmap), # t = 4, Circle
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 1, 1), "", nmap, unsigned_cmap), # t = 7, Circle
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 1, 2), "", nmap, unsigned_cmap), # t = 10, Circle

    particles_field("Ions", "Density", "Z", subplot(fig, gs, 3, 0), "", nmap, unsigned_cmap), # t = 4, Square
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 3, 1), "", nmap, unsigned_cmap), # t = 7, Square
    particles_field("Ions", "Density", "Z", subplot(fig, gs, 3, 2), "", nmap, unsigned_cmap), # t = 10, Square
]

time = min(ParamsCircle.time, ParamsSquare.time)

TAU = [4, 7, int(time / tau)]

nmax = [None]*len(bz)

phi = np.linspace(0, 2 * np.pi)
xr = rr * np.cos(phi)
yr = rr * np.sin(phi)

def cbar(d, n, cax):
    fig = d.axes_position.figure
    cax.tick_params(labelsize=Fonts.ssmol, pad=6)
    fig.colorbar(
        d.im,
        orientation="horizontal",
        cax=cax,
        ticks=np.linspace(d.vmin_vmax[0], d.vmin_vmax[1], n),
    )

def draw(d, i, j):
    d.draw()

    if d == bz[2]:
        cbar(d, 5, subplot(fig, gs, 0, 3))
    if d == bz[5]:
        cbar(d, 5, subplot(fig, gs, 2, 3))
    if d == ns[2]:
        cbar(d, 6, subplot(fig, gs, 1, 3))
    if d == ns[5]:
        cbar(d, 6, subplot(fig, gs, 3, 3))
        
    d.draw_info()
    ax = d.axes_position
    ax.set_aspect(1)
    ax.plot(xr, yr, linewidth=1, linestyle="--", c="black")
    
    nmax[j] = np.argmax(ns[j].data)
    y0, x0 = np.unravel_index(nmax[j], ns[j].data.shape)
    xc = ns[j].data.shape[0] / 2
    x0 = (x0 - xc) * dx
    y0 = (y0 - xc) * dx
    ax.scatter(x0, y0, c="black", s=6*1.6, zorder=10)
    ax.scatter(x0, y0, c="red",   s=6*1.0, zorder=10)

bbox2 = dict(bbox)
bbox2["linewidth"] = 0.5

for i, tt in enumerate(TAU):
    t = int(tt * tau / dts)

    if i == 0:
        fig.text(0.29, 0.98, "\\rm Cylindrical conductive wall", size=Fonts.big, ha="center", bbox=bbox)
        fig.text(0.76, 0.98, "\\rm Square conductive wall", size=Fonts.big, ha="center", bbox=bbox)

    for j in (i, i + 3):
        bz[j].axes_args["title"] = "$B_z(x, \\, y)$"
        ns[j].axes_args["title"] = "$n_i(x, \\, y)$"

        if i < 2:
            bz[j].axes_args.pop("xlabel")
            ns[j].axes_args.pop("xlabel")
        
        ns[j].axes_args.pop("ylabel")

        if j == i + 3:
            bz[j].axes_args.pop("ylabel")


        if j == i:
            fig.text(0.043, 0.9605 - 0.1937 * i, f"$t/\\tau = {tt}$", size=Fonts.big, ha="center", bbox=bbox2)
            prefix = get_prefix(t, ParamsCircle.restart_timesteps, ParamsCircle.prefixes)        
        else:
            prefix = get_prefix(t, ParamsSquare.restart_timesteps, ParamsSquare.prefixes)        
        
        bz[j].data = get_parsed_field(bz[j], "B", "Z", "z", t, prefix) 
        ns[j].data = get_parsed_scalar(ns[j], t, prefix) 

        draw(bz[j], i, j)
        draw(ns[j], i, j)

############# FOURIER ############

TIME = int(10 * tau / dts)
TS = np.arange(0, TIME) * dts / tau
R = 10.0
MAX_M = 5

def prepare_field_m(ax, label, data):
    F_tm = Field(None, ax)

    F_tm.set_axes_args(
        xlabel="$t / \\tau$",
        xlim=(0, TIME * dts / tau),
        xticks=np.linspace(0, TIME * dts / tau, 6),
    )

    ## Merge it back into `fourier_transform()` ##
    f_data = np.fft.ifftshift(data, axes=1)
    f_data = np.fft.fft(f_data, axis=1)
    f_data = np.fft.fftshift(f_data, axes=1)
    F_tm.data = np.abs(f_data)

    shape = data.shape
    k = np.fft.fftfreq(shape[1], d=(2 * np.pi /  shape[1])) * (2 * np.pi)
    k = np.fft.fftshift(k)
    m = np.round(k)
    ##

    return F_tm, m

def avg(d):
    return sliding_average(d, 30)

def get_m(d, i):
    return d[:, np.argwhere(np.abs(m - i) < 0.5)[0][0]]


compare_models = [
    ("\\rm Conductive, circle", ParamsCircle.params_path),
    ("\\rm Conductive, square", ParamsSquare.params_path),
]

for j, (label, path) in enumerate(compare_models):
    ax = fig.add_subplot(gs[4, j*2:j*2+2])

    F_at = prepare_field_phit(None, "ni", "", R, 1.0, path)
    F_mw, m = prepare_field_m(ax, label, F_at.data)

    for i in np.arange(1, MAX_M + 1):
        ax.plot(avg(TS), avg(get_m(F_mw.data, i)), label=f"$m = {i}$", linewidth=1.5 if i == 1 else 1)
        ax.set_yscale("log")

    F_mw.draw_info()

    V = np.sqrt(T_i / mi_me)
    L = 80 / 2
    tt = TS[len(TS)//2 - 1000:]
    Ga = 0.5929 
    Gr = (Ga / tau) * (L / V)

    if j == 0:
        ax.plot(tt, 0.38 * np.exp(Ga * tt), linestyle="--", color="C0", linewidth=1.5, label=f"$\\Gamma \\approx {Gr:.2f} \\, \\Gamma_{{fl}}$")
        # ax.set_title(label, pad=0, y=1.05)
        ax.legend(loc="upper left", fontsize=8.5, framealpha=1.0) #bbox_to_anchor=(0.44, 0.55),
    elif j == 1:
        ax.plot(tt, 0.48 * np.exp(Ga * tt), linestyle="--", color="C0", linewidth=1.5)

    ax.set_ylim(1e-2, 1e+3)
    ax.grid(alpha=0.6)

fig.text(0.415, 0.290, f"$\\delta n_i(t, m, r = {R:.1f})$", size=Fonts.big, ha="center", bbox=bbox2)

fig.tight_layout(h_pad=0.1, w_pad=0.1, rect=(0, 0, 1, 0.98))
fig.savefig(f"{res_dir}/os8.pdf")
