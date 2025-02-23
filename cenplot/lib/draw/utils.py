import numpy as np
import matplotlib.pyplot as plt

from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

from .settings import PlotSettings
from ..utils import Unit
from ..track.types import LegendPosition, Track, TrackType, TrackPosition
from ..track.settings import DefaultTrackSettings


def create_subplots(
    tracks: list[Track],
    settings: PlotSettings,
    **kwargs: Any,
) -> tuple[Figure, np.ndarray, dict[int, int]]:
    """
    Generate a subplot figure from a list of `Track`s.

    # Args
    * `tracks`
        * Input `Track`s
    * `settings`
        * Plot settings.
    * `kwargs`
        * Additional arguments passed to `plt.subplots`

    # Returns
    * Figure, its axes, and a mapping of the idx of `tracks` to the row idx of axes.
    """
    track_props = []
    track_indices = {}
    requires_second_col = False
    legend_prop = settings.legend_prop

    track_idx = 0
    for i, track in enumerate(tracks):
        # Store index.
        # Only increment index if takes up a subplot axis.
        if track.pos == TrackPosition.Relative:
            track_indices[i] = track_idx
            track_idx += 1
            track_props.append(track.prop)
        # For each unique HOR monomer number, create a new track.
        # Divide the proportion of the image allocated between each mer track.
        elif track.opt == TrackType.HORSplit:
            if track.options.mode == "hor":
                n_subplots = track.data["name"].unique()
            else:
                n_subplots = track.data["mer"].unique()
            for j, _ in enumerate(n_subplots):
                track_indices[i + j] = track_idx
                track_props.append(track.prop)
                track_idx += 1
        else:
            track_indices[i] = track_idx - 1

        if not requires_second_col and track.options.legend:
            requires_second_col = True

    # Adjust columns and width ratio.
    num_cols = 2 if requires_second_col else 1
    if settings.legend_pos == LegendPosition.Left:
        width_ratios = (legend_prop, 1 - legend_prop) if requires_second_col else [1.0]
    else:
        width_ratios = (1 - legend_prop, legend_prop) if requires_second_col else [1.0]

    fig, axes = plt.subplots(
        # Count number of tracks
        len(track_props),
        num_cols,
        figsize=settings.dim,
        height_ratios=track_props,
        width_ratios=width_ratios,
        # Always return 2D ndarray
        squeeze=0,
        layout=settings.layout,
        **kwargs,
    )

    return fig, axes, track_indices


def merge_plots(
    figures: list[tuple[Figure, np.ndarray, list[str]]], outfile: str
) -> None:
    """
    Merge plots produced by `plot_one_cen`.

    # Args
    * `figures`
        * List of figures, their axes, and the name of the output files. Only pngs are concatentated.
    * `outfile`
        * Output merged file.
        * Either `png` or `pdf`

    # Returns
    * None
    """
    if outfile.endswith(".pdf"):
        with PdfPages(outfile) as pdf:
            for fig, _, _ in figures:
                pdf.savefig(fig)
    else:
        merged_images = np.concatenate(
            [
                plt.imread(file)
                for _, _, files in figures
                for file in files
                if file.endswith("png")
            ]
        )
        plt.imsave(outfile, merged_images)


def format_ax(
    ax: Axes,
    *,
    grid=False,
    xticks: bool = False,
    xticklabel_fontsize: float | str | None = None,
    yticks: bool = False,
    yticklabel_fontsize: float | str | None = None,
    spines: tuple[str, ...] | None = None,
) -> None:
    """
    Format an axis by removing elements and adjusting fontsize.
    """
    if grid:
        ax.grid(False)
    if xticks:
        ax.set_xticks([], [])
    if xticklabel_fontsize:
        for lbl in ax.get_xticklabels():
            lbl.set_fontsize(xticklabel_fontsize)
    if yticks:
        ax.set_yticks([], [])
    if yticklabel_fontsize:
        for lbl in ax.get_yticklabels():
            lbl.set_fontsize(yticklabel_fontsize)
    if spines:
        for spine in spines:
            ax.spines[spine].set_visible(False)


def set_both_labels(y_lbl: str, ax: Axes, track: Track):
    # Set y-label.
    if track.title:
        ax.set_ylabel(
            y_lbl,
            rotation="horizontal",
            ha="right",
            va="center",
            ma="center",
            fontsize=track.options.title_fontsize,
        )
    # Set x-label.
    if not track.options.hide_x:
        format_xaxis_ticklabels(ax, track)


def format_xaxis_ticklabels(ax: Axes, track: Track):
    """
    Format x-axis ticklabels with `Track.options.units_x`.
    """
    # Remove scientific notation.
    ax.xaxis.set_major_formatter("plain")
    new_xtick_labels = []
    units = Unit(track.options.units_x)
    xmin, xmax = ax.get_xlim()
    xticks, xticklabels = ax.get_xticks(), ax.get_xticklabels()
    for txt in xticklabels:
        x, _ = txt.get_position()
        # Convert units and round.
        new_x_txt = units.convert_value(x, 3)
        txt.set_text(new_x_txt)
        new_xtick_labels.append(txt)

    ax.set_xticks(xticks, new_xtick_labels)
    ax.set_xlabel(
        f"Position ({units.capitalize()})", fontsize=track.options.title_fontsize
    )
    ax.set_xlim(xmin, xmax)


def draw_uniq_entry_legend(
    ax: Axes,
    track: Track,
    ref_ax: Axes | None = None,
    ncols: int | None = DefaultTrackSettings.legend_ncols,
    **kwargs: Any,
) -> None:
    """
    Draw a legend with only unique entries.

    # Args
    * `ax`
        * Axis to plot on.
    * `track`
        * `Track` legend to plot.
    * `ref_ax`
        * Reference axis to use in getting legend handles and labels.
    * `ncols`
        * Number of columns for legend.
    * `kwargs`
        * Additional parameters passed to `ax.legend`

    # Returns
    * None
    """
    ref_ax = ref_ax if ref_ax else ax

    # Dedupe labels.
    handles, labels = ref_ax.get_legend_handles_labels()
    by_label: dict[str, Rectangle] = dict(zip(labels, handles))

    if not ncols:
        ncols = 4

    if not track.options.legend_title_only:
        legend = ax.legend(
            by_label.values(),
            by_label.keys(),
            ncols=ncols,
            # Set aspect ratio of handles so square.
            handlelength=1.0,
            handleheight=1.0,
            frameon=False,
            **kwargs,
        )

        # Set patches edge color manually.
        # Turns out get_legend_handles_labels will get all rect patches and setting linewidth will cause all patches to be black.
        for ptch in legend.get_patches():
            ptch.set_linewidth(1.0)
            ptch.set_edgecolor("black")
    else:
        legend = ax.legend([], [], frameon=False, loc="center left", alignment="left")

    # Set legend title.
    if track.options.legend_title:
        legend.set_title(track.options.legend_title)
        legend.get_title().set_fontsize(track.options.legend_title_fontsize)


def add_border(ax: Axes, height: float, zorder: float):
    xmin, xmax = ax.get_xlim()
    rect = Rectangle(
        (xmin, 0),
        xmax - xmin,
        height,
        fill=None,
        zorder=zorder,
    )
    ax.add_patch(rect)
