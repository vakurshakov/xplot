import numpy as np

from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot_utils import *

# Main class for plotting
# @todo Should be better parametrized

plt.rcParams.update({"text.usetex": True, "axes.formatter.use_mathtext": True})

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
        cbar_exponential = False
        for name, arg in kwargs.items():
            if name == "add_cbar":
                add_cbar = arg
            if name == "cbar_orientation":
                cbar_orientation = arg
            if name == "cbar_pad":
                cbar_pad = arg
            if name == "cbar_ticks_num":
                cbar_ticks_num = arg
            if name == "cbar_exponential":
                cbar_exponential = arg

        if (self.cbar != None and self.vmin_vmax != None) or not add_cbar:
            return

        fig = self.axes_position.figure
        divider = make_axes_locatable(self.axes_position)

        cbar_side = "right" if (cbar_orientation == "vertical") else "bottom"
        cax = divider.append_axes(cbar_side, size="5%", pad=cbar_pad)
        cax.tick_params(labelsize=ticksize, pad=10)

        self.cbar = fig.colorbar(
            self.im,
            orientation=cbar_orientation,
            cax=cax,
            ticks=np.linspace(vmin, vmax, cbar_ticks_num),
        )

        if cbar_exponential:
            self.cbar.formatter.set_powerlimits((0, 0))
            self.cbar.ax.yaxis.get_offset_text().set_visible(False)
            annotate_x(self.cbar.ax, "$\\times\\,10^{" + str(find_exp(self.vmin_vmax[1])) + "}$", y=1.02, size=(ticksize * 0.8), bbox=None)


    def clear(self):
        self.axes_position.cla()


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
        ax.tick_params(labelsize=ticksize, pad=8)

        for name, arg in {**self.axes_args, **kwargs}.items():
            if name == "title":
                ax.set_title(arg, fontsize=titlesize, pad=0, y=1.05)
            elif name == "xlim":
                ax.set_xlim(arg)
            elif name == "ylim":
                ax.set_ylim(arg)
            elif name == "zlim":
                ax.set_zlim(arg)
            elif name == "xlabel":
                ax.set_xlabel(arg, fontsize=labelsize, labelpad=12)
            elif name == "ylabel":
                ax.set_ylabel(arg, fontsize=labelsize, labelpad=10)
            elif name == "zlabel":
                ax.set_zlabel(arg, fontsize=labelsize, labelpad=16)
            elif name == "xticks":
                ax.set_xticks(arg)
            elif name == "yticks":
                ax.set_yticks(arg)
            elif name == "zticks":
                ax.set_zticks(arg)
            elif name == "xticklabels":
                ax.set_xticklabels(arg, fontsize=labelsize)
            elif name == "yticklabels":
                ax.set_yticklabels(arg, fontsize=labelsize)
            elif name == "zticklabels":
                ax.set_zticklabels(arg, fontsize=labelsize)
            else:
                raise RuntimeError("No parameter named " + name)
