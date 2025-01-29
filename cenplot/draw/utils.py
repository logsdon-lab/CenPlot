import numpy as np
import matplotlib.pyplot as plt

from typing import Any

from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.backends.backend_pdf import PdfPages

from ..utils import Unit
from ..track.types import LegendPosition, Track, TrackOption, TrackPosition


def create_subplots(
    dfs_track: list[Track],
    width: float,
    height: float,
    legend_pos: LegendPosition,
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
            uniq_mers = track.data["mer"].unique()
            track_prop = track.prop / len(uniq_mers)
            for j, _ in enumerate(uniq_mers):
                track_indices[i + j] = track_idx
                track_props.append(track_prop)
                track_idx += 1
        else:
            track_indices[i] = track_idx - 1

        if not requires_second_col and track.options.legend:
            requires_second_col = True

    # Adjust columns and width ratio.
    num_cols = 2 if requires_second_col else 1
    if legend_pos == LegendPosition.Left:
        width_ratios = (0.2, 0.8) if requires_second_col else [1.0]
    else:
        width_ratios = (0.8, 0.2) if requires_second_col else [1.0]

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


def minimalize_ax(
    ax: Axes,
    *,
    grid=False,
    xticks: bool = False,
    yticks: bool = False,
    spines: tuple[str, ...] | None = None,
) -> None:
    if grid:
        ax.grid(False)
    if xticks:
        ax.set_xticks([], [])
    if yticks:
        ax.set_yticks([], [])
    if spines:
        for spine in spines:
            ax.spines[spine].set_visible(False)


def set_position_xlabel(ax: Axes):
    xmin, xmax = ax.get_xlim()
    xlen = xmax - xmin
    if (xlen / 1_000_000) > 1:
        unit = Unit.Mbp
    elif (xlen / 1_000) > 1:
        unit = Unit.Kbp
    else:
        unit = Unit.Bp
    ax.set_xlabel(f"Position ({unit.capitalize()})")


def draw_uniq_entry_legend(ax: Axes, ref_ax: Axes | None = None, **kwargs: Any) -> None:
    ref_ax = ref_ax if ref_ax else ax
    handles, labels = ref_ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), frameon=False, **kwargs)
