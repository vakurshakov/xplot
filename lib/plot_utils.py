from math import log10, floor

import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec

# Utilities to set font sizes externally
def set_big(new_big):
    global big
    big = new_big

def set_smol(new_smol):
    global smol
    smol = new_smol

def set_ssmol(new_ssmol):
    global ssmol
    ssmol = new_ssmol

def subplot(fig, gs, x=0, y=0):
    return fig.add_subplot(gs[x + y * gs.ncols])

bbox = dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.25')

def annotate_x(axis, annotation, x=0.5, y=1, size=big, ha="center", bbox=bbox):
    axis.annotate(
        annotation,
        xy=(x, y),
        xytext=(0, 1),
        xycoords="axes fraction",
        textcoords="offset points",
        ha=ha,
        va="baseline",
        size=size,
        bbox=bbox
    )

def annotate_y(axis, annotation):
    axis.annotate(
        annotation,
        xy=(0, 0.5),
        xytext=(-axis.yaxis.labelpad - 1, 0),
        xycoords=axis.yaxis.label,
        textcoords="offset points",
        ha="right",
        va="center",
        rotation=90,
        size=big,
    )

def figure(ncols=1, nrows=1, width_ratios=None, height_ratios=None):
    fig = plt.figure(figsize=(8 * ncols * 1.2, 8 * nrows * 1.2))

    if width_ratios == None:
        width_ratios = [1] * ncols
    if height_ratios == None:
        height_ratios = [1] * nrows

    gs = plt.GridSpec(
        ncols=ncols,
        nrows=nrows,
        width_ratios=width_ratios,
        height_ratios=height_ratios,
        figure=fig
    )

    return fig, gs

def find_exp(number):
    return int(floor(log10(abs(number))))

# Custom colormaps that can be used as `cmap` parameters anywhere
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
