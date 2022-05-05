from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt
import matplotlib as mpl
from importlib import resources
import pandas as pd
import seaborn as sns

# Color palettes
class colors:
    """
    Color palettes for plots
    """

    flatDesign = plt.cycler(
        color=["#e27a3d", "#344d5c", "#df5a49", "#43b29d", "#efc94d"]
    )

    firenze = plt.cycler(color=["#8E2800", "#468966", "#B64926", "#FFF0A5", "#FFB03B"])

    vitaminC = plt.cycler(color=["#FD7400", "#004358", "#FFE11A", "#1F8A70", "#BEDB39"])

    bella = plt.cycler(color=["#801637", "#047878", "#FFB733", "#F57336", "#C22121"])

    buddha = plt.cycler(color=["#192B33", "#FF8000", "#8FB359", "#FFD933", "#CCCC52"])

    elemental = plt.cycler(
        color=["#E64661", "#FFA644", "#998A2F", "#2C594F", "#002D40"]
    )

    carolina = plt.cycler(color=["#73839C", "#2E4569", "#AECCCF", "#D5957D", "#9C7873"])

    fourtyTwo = plt.cycler(
        color=["#2469A6", "#C4E1F2", "#F2E205", "#F2D22E", "#D9653B"]
    )

    terrazaverde = plt.cycler(
        color=["#DFE2F2", "#88ABF2", "#4384D9", "#56BFAC", "#D9B341"]
    )


def layout(
    fontSize=16,
    axTitleSize=16,
    axLabelSize=16,
    tickLabelSize=12,
    colors=colors.firenze,
):

    plt.rcParams["figure.constrained_layout.use"] = True
    plt.rcParams["savefig.dpi"] = 300

    plt.rc("figure", figsize=(8, 7), facecolor="white")

    # Text
    plt.rc("font", family="sans-serif", size=fontSize)

    # Legend
    plt.rc("legend", fontsize=fontSize / 1.5, fancybox=False, facecolor="white")

    # Axes
    plt.rc("xtick", direction="in", labelsize=tickLabelSize)
    plt.rc("ytick", direction="in", labelsize=tickLabelSize)
    plt.rc(
        "axes",
        grid=True,
        titlesize=axTitleSize,
        labelsize=axLabelSize,
        axisbelow=True,
        linewidth=1.5,
        prop_cycle=colors,
        facecolor="whitesmoke",
    )
    plt.rc("grid", color="snow")

    # Lines
    plt.rc("lines", linewidth=4, markersize=10, markeredgecolor="k")


def TAS(labels=False, fontsize="medium", **kwargs):
    """Returns a line plot element of classification of volcanic rocks
    in total-alkali vs silica plots
    """

    with resources.open_text("meltInc.static", "TAS.csv") as df:
        TAS = pd.read_csv(df)

    rock_labels = {
        "Picro-basalt": ["Picro\nbasalt", [41.7, 1.5]],
        "Basalt": ["Basalt", [47, 2.5]],
        "Basaltic andesite": ["Basaltic\nandesite", [53, 2.5]],
        "Andesite": ["Andesite", [58, 2.5]],
        "Dacite": ["Dacite", [65.5, 4]],
        "Trachy-basalt": ["Trachy-\nbasalt", [47.5, 5.5]],
        "Basaltic trachy-andesite": ["Basaltic\ntrachy-\nandesite", [51.6, 6.5]],
        "Trachy-andesite": ["Trachy-\nandesite", [56, 8]],
        "Trachyte": ["Trachyte", [64, 11]],
        "Tephrite": ["Tephrite", [43.5, 7]],
        "Phono-tephrite": ["Phono-\ntephrite", [47, 9.0]],
        "Tephri-phonolite": ["Tephri-\nphonolite", [51, 11]],
        "Phonolite": ["Phonolite", [55, 15]],
        "Foidite": ["Foidite", [45, 14]],
        "Rhyolite": ["Rhyolite", [72, 8.5]],
    }

    if labels:
        for _, rock in rock_labels.items():
            plt.text(
                *rock[1],
                rock[0],
                fontsize=fontsize,
                fontfamily="monospace",
                clip_on=True
            )

    for id in TAS.id.unique():
        plt.plot(
            TAS.loc[TAS.id == id, "x"],
            TAS.loc[TAS.id == id, "y"],
            "-",
            color="k",
            **kwargs
        )


def side_plots(fig, x_axis: bool, y_axis: bool, side=10, **kwargs):

    widths = [(100,), (100. - side, side)][y_axis]
    heights = [(100,), (side, 100. - side)][x_axis]

    left, right = 0.12, 0.97
    top, bottom = 0.97, 0.1
    spacing = 0.03

    grid = fig.add_gridspec(
        ncols=1 + y_axis,
        nrows=1 + x_axis,
        width_ratios=widths,
        height_ratios=heights,
        left=left,
        right=right,
        bottom=bottom,
        top=top,
        hspace=spacing,
        wspace=spacing,
    )

    ax = fig.add_subplot(grid[0 + x_axis, 0])

    axes = []

    if x_axis:        
        ax_kde_x = fig.add_subplot(grid[0, 0], sharex=ax)
        axes.append(ax_kde_x)

        ax_kde_x.spines["left"].set_visible(False)
        ax_kde_x.yaxis.set_visible(False)
        ax_kde_x.tick_params(axis="x", labelbottom=False, direction="in")


    if y_axis:
        ax_kde_y = fig.add_subplot(grid[0 + x_axis, 1], sharey=ax)  
        axes.append(ax_kde_y) 
        
        ax_kde_y.xaxis.set_visible(False)
        ax_kde_y.spines["bottom"].set_visible(False)
        ax_kde_y.tick_params(axis="y", labelleft=False, direction="in")

    
    for axis in axes:
        # axis.set_frame_on(False)
        axis.grid(False)
        axis.set_facecolor(color="white")

        for spine in ["right", "top"]:
            axis.spines[spine].set_visible(False)

    return ax, *axes

