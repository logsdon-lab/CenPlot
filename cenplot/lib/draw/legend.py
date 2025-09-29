import sys
import numpy as np
from typing import Any
from matplotlib.axes import Axes
from matplotlib.artist import Artist

from ..track.types import Track
from ..draw.utils import format_ax


def draw_legend(
    ax: Axes,
    axes: np.ndarray,
    track: Track,
    ref_track_col: int,
) -> None:
    """
    Draw legend plot on axis for the given `Track`.

    # Args
    * `ax`
        * Axis to plot on.
    * `axes`
        * 2D `np.ndarray` of all axes to get reference axis.
    * `track`
        * Current `Track`.
    * `track_col`
        * Reference `Track` column.

    # Returns
    * None
    """
    if isinstance(track.options.index, int):
        ref_track_rows = [track.options.index]
    elif isinstance(track.options.index, list):
        ref_track_rows = track.options.index
    else:
        raise ValueError("Invalid type for reference legend indices.")

    all_label_handles: dict[Any, Artist] = {}
    for row in ref_track_rows:
        try:
            ref_track_ax: Axes = axes[row, ref_track_col]
        except IndexError:
            print(f"Reference axis index ({row}) doesn't exist.", sys.stderr)
            continue

        handles, labels = ref_track_ax.get_legend_handles_labels()
        labels_handles: dict[Any, Artist] = dict(zip(labels, handles))
        all_label_handles = all_label_handles | labels_handles

    # Provide custom order.
    # Keeps all elements as opposed to draw_uniq_entry_legend
    if track.options.legend_label_order:
        legend_labels = set(all_label_handles.keys())
        new_all_label_handles = {}
        for label in track.options.legend_label_order:
            if not all_label_handles.get(label):
                continue
            new_all_label_handles[label] = all_label_handles[label]

        remaining_labels = legend_labels.difference(new_all_label_handles.keys())
        for label in remaining_labels:
            new_all_label_handles[label] = all_label_handles[label]

        all_label_handles = new_all_label_handles

    # Some code dup.
    if not track.options.legend_title_only:
        legend = ax.legend(
            all_label_handles.values(),
            all_label_handles.keys(),
            ncols=track.options.legend_ncols if track.options.legend_ncols else 10,
            # Set aspect ratio of handles so square.
            handlelength=1.0,
            handleheight=1.0,
            frameon=False,
            fontsize=track.options.legend_fontsize,
            loc="center",
            alignment="center",
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

    format_ax(
        ax,
        grid=True,
        xticks=True,
        yticks=True,
        spines=("right", "left", "top", "bottom"),
    )
