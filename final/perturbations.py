#!/usr/bin/env python3

from final import *
from xplot.fourier.fourier import *

set_big(13)
set_smol(11)
set_ssmol(11)

plt.rc('text', usetex=True)
plt.rc('axes', titlesize=13, labelsize=13)
plt.rc('xtick', labelsize=11)
plt.rc('ytick', labelsize=11)
plt.rc('legend', fontsize=11.5)
plt.rc('figure', titlesize=13)
plt.rc('lines', linewidth=1.3)

ncols=1
nrows=1

fig = plt.figure(figsize=(6.0, 4.0))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1]*ncols, height_ratios=[1]*nrows, figure=fig)

compare_models = [
    ("Damping, circle", "../T11_MergeV"),
    # ("Conductive, circle", "../ConductiveWall_FixRotor/Circle"),
    # ("Conductive, square", "../ConductiveWall_FixRotor/Square"),
]

TIME = d_tmax
TS = np.arange(0, TIME) * dts / tau
R = 10.0
MAX_M = 5

def prepare_field_m(ax, label, data):
    F_tm = Field(None, ax)

    F_tm.set_axes_args(
        title=label,
        ylabel=f"$\\delta n_i(t, m, r = {R:.1f})$",
        xlabel="$t / \\tau$",
        xlim=(0, TIME * dts / tau),
        xticks=np.linspace(0, TIME * dts / tau, 7),
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

for j, (label, path) in enumerate(compare_models):
    ax = subplot(fig, gs, 0, j)

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
        ax.plot(tt, 0.60 * np.exp(Ga * tt), linestyle="--", color="C0", linewidth=1.5, label=f"$\\Gamma \\approx {Gr:.2f} \\, v_{{Ti}} / L$")
        ax.set_title(label, pad=0, y=1.05)
        ax.legend(bbox_to_anchor=(0.37, 0.55), framealpha=1.0) #

    ax.set_ylim(1e-2, 1e+3)
    ax.grid(alpha=0.6)

fig.tight_layout(pad=0.25)
fig.savefig(f"{res_dir}/perturbations_R{R:.1f}.pdf")