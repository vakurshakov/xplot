import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

fancy = True
if (fancy):
    plt.rcParams.update({"text.usetex": True, "axes.formatter.use_mathtext": True})

cdict = {
    "red": (
        (0.0, 0.0, 0),
        (0.35, 0.5859375, 0.5859375),
        (0.4, 0, 0),    # blue
        (0.45, 0, 0),   # 0096ff
        (0.495, 1, 1),  # white
        (0.505, 1, 1),  # white
        (0.55, 1, 1),   # ff9600
        (0.65, 1, 1),
        (1.0, 0.5859375, 0.5859375),
    ),
    "green": (
        (0.0, 0.0, 0.0),
        (0.35, 0, 0),
        (0.4, 0, 0),                   # blue
        (0.45, 0.5859375, 0.5859375),  # 0096ff
        (0.495, 1, 1),                 # white
        (0.505, 1, 1),                 # white
        (0.55, 0.5859375, 0.5859375),  # ff9600
        (0.65, 0, 0),
        (1.0, 0, 0),
    ),
    "blue": (
        (0.0, 0.5859375, 0.5859375),
        (0.35, 1, 1),
        (0.4, 1, 1),    # blue
        (0.45, 1, 1),   # 0096ff
        (0.495, 1, 1),  # white
        (0.505, 1, 1),  # white
        (0.55, 0, 0),   # ff9600
        (0.65, 0.5859375, 0.5859375),
        (1.0, 0, 0),
    ),
}
signed_cmap = col.LinearSegmentedColormap("signed_cmap", cdict, N=256, gamma=1)

cdict = {
    "red": (
        (0.0, 1, 1),   # white
        (0.01, 1, 1),  # white
        (0.15, 0, 0),  # 0096ff
        (0.35, 0, 0),  # blue
        (0.55, 1, 1),  # ff9600
        (0.75, 1, 1),
        (1.0, 0.5859375, 0.5859375),
    ),
    "green": (
        (0.0, 1, 1),                   # white
        (0.01, 1, 1),                  # white
        (0.15, 0.5859375, 0.5859375),  # 0096ff
        (0.35, 0, 0),                  # blue
        (0.55, 0.5859375, 0.5859375),  # ff9600
        (0.75, 0, 0),
        (1.0, 0, 0),
    ),
    "blue": (
        (0.0, 1, 1),   # white
        (0.01, 1, 1),  # white
        (0.15, 1, 1),  # 0096ff
        (0.35, 1, 1),  # blue
        (0.55, 0, 0),  # ff9600
        (0.75, 0.5859375, 0.5859375),
        (1.0, 0, 0),
    ),
}
unsigned_cmap = col.LinearSegmentedColormap("unsigned_cmap", cdict, N=256, gamma=1)

if (fancy):
    big = 36
    smol = 34
    ssmol = 30
else:
    big = 30
    smol = 28
    ssmol = 24

def set_big(new_big):
    global big
    big = new_big

def set_smol(new_smol):
    global smol
    smol = new_smol

def set_ssmol(new_ssmol):
    global ssmol
    ssmol = new_ssmol


class Field:
    def __init__(
        self,
        path_to_file,
        axes_position=None,
        boundaries=(0, 1, 0, 1),
        cmap="plasma",
        vmin_vmax=(0, 1),
    ):
        self.path_to_file = path_to_file
        self.axes_position = axes_position
        self.data = None
        self.boundaries = boundaries
        self.cmap = cmap
        self.vmin_vmax = vmin_vmax
        self.cbar = None
        self.axes_args = {}
        self.im = None

    def read_data(self, t: int, **kwargs):
        data_shape = (
            int(self.boundaries[3] - self.boundaries[2]),
            int(self.boundaries[1] - self.boundaries[0]),
        )

        for name, arg in kwargs.items():
            if name == "data_shape":
                data_shape = arg
            else:
                raise RuntimeError("No argument named " + name)

        self.data = np.zeros(data_shape)

        with open(self.path_to_file + f"{t}.bin", "rb") as f:
            temp = np.fromfile(
                f, dtype=np.float32, count=int(data_shape[0] * data_shape[1]), offset=0
            )

            np.copyto(self.data, temp.reshape(data_shape))

    def draw(self, **kwargs):
        if self.boundaries == None:
            self.boundaries = (0, self.data.shape[1], 0, self.data.shape[0])

        if self.vmin_vmax == None:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            vmin = self.vmin_vmax[0]
            vmax = self.vmin_vmax[1]

        self.im = self.axes_position.imshow(
            self.data,
            cmap=self.cmap,
            interpolation="gaussian",
            origin="lower",
            aspect="auto",
            extent=(
                self.boundaries[0],
                self.boundaries[1],
                self.boundaries[2],
                self.boundaries[3],
            ),
            vmin=vmin,
            vmax=vmax,
        )


        add_cbar = None
        cbar_pad = 0.2
        cbar_orientation = "vertical"
        cbar_ticks_num = 5
        for name, arg in kwargs.items():
            if name == "add_cbar":
                add_cbar = arg
            if name == "cbar_orientation":
                cbar_orientation = arg
            if name == "cbar_pad":
                cbar_pad = arg
            if name == "cbar_ticks_num":
                cbar_ticks_num = arg

        if (self.cbar != None and self.vmin_vmax != None) or not add_cbar:
            return

        fig = self.axes_position.figure
        divider = make_axes_locatable(self.axes_position)

        cbar_side = "right" if (cbar_orientation == "vertical") else "bottom"
        cax = divider.append_axes(cbar_side, size="5%", pad=cbar_pad)
        cax.tick_params(labelsize=ssmol, pad=10)

        self.cbar = fig.colorbar(
            self.im,
            orientation=cbar_orientation,
            cax=cax,
            ticks=np.linspace(vmin, vmax, cbar_ticks_num),
        )

    def draw_box(self, **kwargs):
        if self.boundaries == None:
            self.boundaries = (0, self.data.shape[0], 0, self.data.shape[1])

        if self.vmin_vmax == None:
            vmin = self.data.min()
            vmax = self.data.max()
        else:
            vmin = self.vmin_vmax[0]
            vmax = self.vmin_vmax[1]

        for name, arg in kwargs.items():
            if name == "zdir":
                zdir = arg

        if zdir == 'x':
            y = np.arange(0, 153, 1) - 153 // 2
            z = np.arange(0, 103, 1)
            y, z = np.meshgrid(y, z)
            x = np.ones(self.data.shape) * (- 153 // 2)

        elif zdir == 'y':
            x = np.arange(0, 153, 1) - 153 // 2
            z = np.arange(0, 103, 1)
            x, z = np.meshgrid(x, z)
            y = np.ones(self.data.shape) * (+ 153 // 2)

        elif zdir == 'z':
            x = y = np.arange(0, 153, 1) - 153 // 2
            x, y = np.meshgrid(x, y)
            z = np.ones(self.data.shape) * 0

        ls = col.LightSource()

        # creating RGBA equivalent of our data
        rgba = ls.shade(
            self.data,
            cmap=self.cmap,
            vmin=vmin,
            vmax=vmax
        )

        # zeroing alpha-channel by hand
        rgba_epsilon = 1e-6
        rgba[np.abs(self.data) < rgba_epsilon, 3] = 0.1

        self.im = self.axes_position.plot_surface(
            x, y, z,
            facecolors=rgba,
            cmap=self.cmap,
            rstride=1,
            cstride=1,
            linewidth=1,
            antialiased=False
        )



    def clear(self):
        self.axes_position.cla()

        if self.cbar != None:
            self.cbar.remove()
            self.cbar = None

    def set_axes_args(self, **kwargs):
        supported_names = [
            "title",
            "xlim",   "ylim",   "zlim",
            "xticks", "yticks", "zticks",
            "xlabel", "ylabel", "zlabel",
            "xticklabels", "yticklabels", "zticklabels",
        ]

        for name, arg in kwargs.items():
            if name in supported_names:
                self.axes_args[name] = arg
            else:
                raise RuntimeError("No parameter named " + name)

    def draw_info(self, **kwargs):
        ax = self.axes_position
        ax.tick_params(labelsize=ssmol, pad=8)

        for name, arg in {**self.axes_args, **kwargs}.items():
            if name == "title":
                ax.set_title(arg, fontsize=big, pad=0, y=1.05)
            elif name == "xlim":
                ax.set_xlim(arg)
            elif name == "ylim":
                ax.set_ylim(arg)
            elif name == "zlim":
                ax.set_zlim(arg)
            elif name == "xlabel":
                ax.set_xlabel(arg, fontsize=smol, labelpad=12)
            elif name == "ylabel":
                ax.set_ylabel(arg, fontsize=smol, labelpad=10)
            elif name == "zlabel":
                ax.set_zlabel(arg, fontsize=smol, labelpad=16)
            elif name == "xticks":
                ax.set_xticks(arg)
            elif name == "yticks":
                ax.set_yticks(arg)
            elif name == "zticks":
                ax.set_zticks(arg)
            elif name == "xticklabels":
                ax.set_xticklabels(arg, fontsize=smol)
            elif name == "yticklabels":
                ax.set_yticklabels(arg, fontsize=smol)
            elif name == "zticklabels":
                ax.set_zticklabels(arg, fontsize=smol)
            else:
                raise RuntimeError("No parameter named " + name)


def return_parameters(parameters_file: str):
    with open(parameters_file, "r") as f:
        f.readline()
        TIME_dt_DTS = f.readline().split(" ")
        TIME = int(TIME_dt_DTS[0])
        dt = float(TIME_dt_DTS[1])
        DTS = int(TIME_dt_DTS[2])

        f.readline()
        bx_ex_dx = f.readline().split(" ")
        bx = int(bx_ex_dx[0])
        ex = int(bx_ex_dx[1])
        dx = float(bx_ex_dx[2])

        f.readline()
        by_ey_dy = f.readline().split(" ")
        by = int(by_ey_dy[0])
        ey = int(by_ey_dy[1])
        dy = float(by_ey_dy[2])

        f.readline()
        sof_char = f.readline().split(" ")
        sizeof_float = int(sof_char[0])

        return sizeof_float, TIME, dt, DTS, bx, ex, dx, by, ey, dy

### plot.py should be better parametrized

