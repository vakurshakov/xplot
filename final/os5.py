#!/usr/bin/env python3

from final import *
from scipy.optimize import curve_fit

ncols=1
nrows=1

fig = plt.figure(figsize=(4.6 * ncols, 4.6 * nrows))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

# vs = np.load(f"{params_path}/Collection/ni_max_value_t.npy", allow_pickle=True)
xs = np.load(f"{params_path}/Collection/ni_max_x_t.npy", allow_pickle=True)
ys = np.load(f"{params_path}/Collection/ni_max_y_t.npy", allow_pickle=True)
rs = np.hypot(xs, ys)
ts = np.arange(0, len(xs)) * dts / tau

ax = subplot(fig, gs, 0, 0)
ax.set_title("\\rm Instability increment", fontsize=16)
ax.set_xlim(ts[0], int(ts[-1]))
ax.set_ylim(0, 10)
ax.set_yticks(np.linspace(*ax.get_ylim(), 6))
ax.tick_params(labelsize=14) #, pad=8)
ax.set_xlabel("\\rm time, $t / \\tau$", fontsize=15)
ax.set_ylabel("\\rm radius, $r$", fontsize=15)

def exponent(t, a, b, c):
  return a * np.exp(+ b * t) + c

popt, pcov = curve_fit(exponent, ts, rs)

# print("a =", popt[0])
# print("b =", popt[1])
# print("c =", popt[2])

V = np.sqrt(T_i / mi_me)
L = 80 / 2
Gamma = (popt[1] / tau) * (L / V)

ax.plot(ts, rs, label="$r^* = {\\rm argmax}(n_i)$")
ax.plot(ts, exponent(ts, *popt), linestyle="--", label=f"$\\Gamma \\approx {Gamma:.2f} \\, \\Gamma_{{fl}}$")

ax.legend(framealpha=1, loc="upper left", fontsize=15)
ax.grid(alpha=0.6)

# plt.show()

fig.tight_layout()
fig.savefig(f"{res_dir}/os5.pdf")
