import os
import sys
import numpy as np

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .hor import draw_hor, draw_hor_ort
from .label import draw_label
from .self_ident import draw_self_ident
from .bar import draw_bars
from .utils import create_subplots, minimalize_ax
from ..track.types import Track, TrackOption, TrackPosition, LegendPosition


def plot_one_cen(
    dfs_track: list[Track],
    outdir: str,
    outfmt: str,
    chrom: str,
    min_st_pos: int,
    max_end_pos: int,
    width: float,
    height: float,
    legend_pos: LegendPosition,
) -> tuple[Figure, np.ndarray, str]:
    # Show chrom trimmed of spaces for logs and filenames.
    print(f"Plotting {chrom}...", file=sys.stderr)

    # # Get min and max position of all tracks for this cen.
    #  = sys.maxsize
    # trk_maxtrk_min_st_end = 0
    # for trk in dfs_track:
    #     if trk.opt == TrackOption.SelfIdent:
    #         continue
    #     try:
    #         trk_min_st = min(trk.data["chrom_st"].min(), trk_min_st)
    #         trk_max_end = max(trk.data["chrom_end"].max(), trk_max_end)
    #     except TypeError:
    #         continue

    # # Scale height based on track length.
    # adj_height = height * (trk_max_end / max_end_pos)
    # height = height if adj_height == 0 else adj_height

    fig, axes, track_indices = create_subplots(
        dfs_track, width, height, legend_pos, constrained_layout=True
    )
    if legend_pos == LegendPosition.Left:
        track_col, legend_col = 1, 0
    else:
        track_col, legend_col = 0, 1

    track_labels: list[str] = []

    def get_track_label(track: Track, all_track_labels: list[str]) -> tuple[str, str]:
        esc_track_label = track.name.encode("unicode_escape").decode("utf-8")
        track_label = track.name.encode("ascii", "ignore").decode("unicode_escape")

        # Update track label for each overlap.
        if track.pos == TrackPosition.Overlap:
            try:
                track_label = f"{all_track_labels[-1]}\n{track_label}"
            except IndexError:
                pass

        return esc_track_label, track_label

    for zorder, track in enumerate(dfs_track):
        track_row = track_indices[track.name]
        esc_track_label, track_label = get_track_label(track, track_labels)

        try:
            track_ax: Axes = axes[track_row, track_col]
        except IndexError:
            print(
                f"Cannot get track ({track_row, track_col}) for {esc_track_label} with {track.pos} position."
            )
            continue
        try:
            legend_ax: Axes = axes[track_row, legend_col]
        except IndexError:
            legend_ax = None

        # Set xaxis limits
        track_ax.set_xlim(min_st_pos, max_end_pos)

        # Switch to line if different track option. {bar, label, ident, hor}
        if track.opt == TrackOption.HOR:
            draw_fn = draw_hor
        elif track.opt == TrackOption.HOROrt:
            draw_fn = draw_hor_ort
        elif track.opt == TrackOption.Label:
            draw_fn = draw_label
        elif track.opt == TrackOption.SelfIdent:
            draw_fn = draw_self_ident
        elif track.opt == TrackOption.Bar:
            draw_fn = draw_bars
        else:
            raise ValueError("Invlaid TrackOption. Unreachable.")

        draw_fn(
            ax=track_ax,
            legend_ax=legend_ax,
            track=track,
            zorder=zorder,
        )

        # Store label if more overlaps.
        track_labels.append(track_label)

        # Set label.
        # Allow chrom as title or name.
        track_label = chrom if track.options.chrom_as_title else track_label
        if track.options.title:
            track_ax.set_ylabel(
                track_label,
                rotation="horizontal",
                ha="right",
                va="center",
                ma="center",
            )

        if not legend_ax:
            continue

        if track.opt != TrackOption.SelfIdent or (
            track.opt == TrackOption.SelfIdent and not track.options.legend
        ):
            # Minimalize all legend cols
            minimalize_ax(
                legend_ax,
                grid=True,
                xticks=True,
                yticks=True,
                spines=("right", "left", "top", "bottom"),
            )
        else:
            minimalize_ax(
                legend_ax,
                grid=True,
                spines=("right", "top"),
            )

    # Add title
    # fig.suptitle(chrom)
    outfile = os.path.join(outdir, f"{chrom}.{outfmt}")

    # Pad between axes.
    fig.get_layout_engine().set(h_pad=0.5)
    fig.savefig(outfile, dpi=600, transparent=True)

    return fig, axes, outfile
