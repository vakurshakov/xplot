#!/usr/bin/env python3

from plot import *

ncols=2
nrows=2

fig = plt.figure(figsize=(8 * ncols * 1.1, 8 * nrows * 1))
gs = GridSpec(ncols=ncols, nrows=nrows, width_ratios=[1] * ncols, height_ratios=[1] * nrows, figure=fig)

res_dir = f"./{params_path}/Fields_filtered"
mkdir(res_dir)

tmin = int(0 * tau / dts)
tmax = int(8 * tau / dts)

vE_max = +0.02
vFE_max = +0.01

m = 1 / 4
Omega_i = B0 / mi_me

def fourier_transform_w(data):
    f_data = np.fft.ifftshift(data, axes=0)
    f_data = np.fft.fft(f_data, axis=0)
    f_data = np.fft.fftshift(f_data, axes=0)

    # np.fft uses linear frequenses, so we multiply the arrays by (2 * np.pi)
    shape = data.shape
    w = np.fft.fftfreq(shape[0], d=dts) * (2 * np.pi)
    w = np.fft.fftshift(w)
    return f_data, w

def inverse_fourier_transform_w(f_data):
    data = np.fft.ifftshift(f_data, axes=0)
    data = np.fft.ifft(data, axis=0)
    data = np.fft.fftshift(data, axes=0)

    r_data = np.real(data)
    i_data = np.imag(data)
    return r_data, i_data

def fourier_filter(data):
    data_kw, w = fourier_transform_w(data)
    data_kw[np.where(np.abs(w) > m * Omega_i),:] = 0.0
    return inverse_fourier_transform_w(data_kw)

# TODO: MPI_Allgather
def recreate_field(path):
    y_range = range(data_shape[0])

    data = np.ndarray((tmax - tmin, len(y_range), data_shape[1]), dtype=np.float32)
    f_data = np.ndarray((tmax - tmin, len(y_range), data_shape[1]), dtype=np.float32)

    for y in y_range:
        print("recreating", path, f"{y:4d} [dy]")
        tmp = np.load(f"{path}_xt_y={y}.npy", allow_pickle=True)
        data[:,(y - y_range[0]),:] = tmp

        f_tmp = fourier_filter(tmp)[0]
        f_data[:,(y - y_range[0]),:] = f_tmp

    return data, f_data

Er = Field(f"./{params_path}/Data2D/Er", subplot(fig, gs, 0, 0), boundaries, signed_cmap, (-vE_max, vE_max))
Ea = Field(f"./{params_path}/Data2D/Ea", subplot(fig, gs, 1, 0), boundaries, signed_cmap, (-vE_max, vE_max))
Er_filtered = Field(None, subplot(fig, gs, 0, 1), boundaries, signed_cmap, (-vFE_max, vFE_max))
Ea_filtered = Field(None, subplot(fig, gs, 1, 1), boundaries, signed_cmap, (-vFE_max, vFE_max))

Er_full, Er_filtered_full  = recreate_field(Er.path_to_file)
Ea_full, Ea_filtered_full  = recreate_field(Ea.path_to_file)

bx = boundaries[0] + buff * dx
ex = boundaries[1] - buff * dx
by = boundaries[2] + buff * dy
ey = boundaries[3] - buff * dy

arg_2d = {
    "xlim": (bx, ex),
    "ylim": (by, ey),
    "xlabel": "$x,~c/\\omega_{pi}$",
    "ylabel": "$y,~c/\\omega_{pi}$",
    "xticks": np.linspace(bx, ex, 5),
    "yticks": np.linspace(by, ey, 5),
}

Er.set_axes_args(**arg_2d, title="$E_r(x, y)$")
Ea.set_axes_args(**arg_2d, title="$E_{\\phi}(x, y)$")

Er_filtered.set_axes_args(**arg_2d, title="$E_r(x, y),~\\omega < " f"{m:.2f}" "\,\\Omega_i$")
Ea_filtered.set_axes_args(**arg_2d, title="$E_{\\phi}(x, y),~\\omega < " f"{m:.2f}" "\,\\Omega_i$")

offset = 5
t_range = np.arange(tmin + rank * offset, tmax, proc * offset)

for t in t_range:
    filename = f"{res_dir}/{str(t // offset).zfill(4)}"
    if os.path.exists(filename):
        continue

    Er.data = Er_full[t - t_range[0]]
    Ea.data = Ea_full[t - t_range[0]]
    Er_filtered.data = Er_filtered_full[t - t_range[0]]  # real part
    Ea_filtered.data = Ea_filtered_full[t - t_range[0]]  # real part

    for diag in [Er, Ea, Er_filtered, Ea_filtered]:
        diag.axes_position.set_aspect(1)
        diag.draw(add_cbar=True)
        diag.draw_info()

    fig.suptitle("$t / \\tau = {" f"{t * dts / tau:.3f}" "}$", x=0.515, size=titlesize, bbox=bbox)

    fig.tight_layout()
    fig.tight_layout()

    print(filename, t, "[dts]", f"{t * dts / tau:.3f}", "[tau]")
    fig.savefig(f"{filename}.png")

    for diag in [Er, Ea]:
        diag.clear()
