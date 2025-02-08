import sys
import numpy as np

from matplotlib.axes import Axes

from ..track.types import Track
from ..draw.utils import draw_uniq_entry_legend, format_ax


def draw_legend(
    ax: Axes,
    axes: np.ndarray,
    track: Track,
    tracks: list[Track],
    track_row: int,
    track_col: int,
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
    * `tracks`
        * All tracks to get reference `Track`.
    * `track_row`
        * Reference track row.
    * `track_col`
        * Reference track col.

    # Returns
    * None
    """
    ref_track_row = (
        track.options.index if isinstance(track.options.index, int) else track_row - 1
    )
    try:
        ref_track_ax: Axes = axes[ref_track_row, track_col]
    except IndexError:
        print(f"Reference axis index ({ref_track_row}) doesn't exist.", sys.stderr)
        return None

    # TODO: Will not work with HOR split.
    legend_colname = (
        tracks[ref_track_row].options.mode
        if hasattr(tracks[ref_track_row].options, "mode")
        else "name"
    )
    try:
        srs_track = tracks[ref_track_row].data[legend_colname]
    except Exception:
        print(f"Legend column ({legend_colname}) doesn't exist in {track}.", sys.stderr)
        return None

    draw_uniq_entry_legend(
        ax,
        track,
        ref_track_ax,
        ncols=track.options.legend_ncols
        if track.options.legend_ncols
        else srs_track.n_unique(),
        loc="center",
        alignment="center",
    )
    format_ax(
        ax,
        grid=True,
        xticks=True,
        yticks=True,
        spines=("right", "left", "top", "bottom"),
    )
