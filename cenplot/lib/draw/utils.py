from matplotlib.patches import Rectangle
import numpy as np
import matplotlib.pyplot as plt

from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

from ..utils import Unit
from ..track.types import LegendPosition, Track, TrackOption, TrackPosition
from ..track.settings import DefaultPlotSettings


def create_subplots(
    dfs_track: list[Track],
    width: float,
    height: float,
    legend_pos: LegendPosition,
    legend_prop: float,
    **kwargs: Any,
) -> tuple[Figure, np.ndarray, dict[int, int]]:
    track_props = []
    track_indices = {}
    requires_second_col = False

    track_idx = 0
    for i, track in enumerate(dfs_track):
        # Store index.
        # Only increment index if takes up a subplot axis.
        if track.pos == TrackPosition.Relative:
            track_indices[i] = track_idx
            track_idx += 1
            track_props.append(track.prop)
        # For each unique HOR monomer number, create a new track.
        # Divide the proportion of the image allocated between each mer track.
        elif track.opt == TrackOption.HORSplit:
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
    if legend_pos == LegendPosition.Left:
        width_ratios = (legend_prop, 1 - legend_prop) if requires_second_col else [1.0]
    else:
        width_ratios = (1 - legend_prop, legend_prop) if requires_second_col else [1.0]

    fig, axes = plt.subplots(
        # Count number of tracks
        len(track_props),
        num_cols,
        figsize=(width, height),
        height_ratios=track_props,
        width_ratios=width_ratios,
        # Always return 2D ndarray
        squeeze=0,
        **kwargs,
    )

    return fig, axes, track_indices


def merge_plots(figures: list[tuple[Figure, np.ndarray, str]], outfile: str) -> None:
    if outfile.endswith(".pdf"):
        with PdfPages(outfile) as pdf:
            for fig, _, _ in figures:
                pdf.savefig(fig)
    else:
        merged_images = np.concatenate([plt.imread(file) for _, _, file in figures])
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
    # Set correct units based on xlim.
    if not track.options.hide_x:
        xmin, xmax = ax.get_xlim()
        xlen = xmax - xmin
        if (xlen / 1_000_000) > 1:
            unit = Unit.Mbp
        else:
            unit = Unit.Bp
        ax.set_xlabel(
            f"Position ({unit.capitalize()})", fontsize=track.options.title_fontsize
        )


def draw_uniq_entry_legend(
    ax: Axes,
    track: Track,
    ref_ax: Axes | None = None,
    ncols: int = DefaultPlotSettings.legend_ncols,
    **kwargs: Any,
) -> None:
    ref_ax = ref_ax if ref_ax else ax
    # Dedupe labels.
    handles, labels = ref_ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))

    legend = ax.legend(
        by_label.values(), by_label.keys(), ncols=ncols, frameon=False, **kwargs
    )

    # Wrap text.
    for txt in legend.get_texts():
        txt.set_wrap(True)

    # Set legend title.
    if track.options.legend_title:
        legend.set_title(track.options.legend_title)
        legend.get_title().set_fontsize(track.options.legend_title_fontsize)
        legend.set_alignment("left")


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
