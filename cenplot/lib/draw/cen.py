import os
import logging
import numpy as np

from matplotlib.axes import Axes
from matplotlib.figure import Figure

from .settings import PlotSettings
from .hor import draw_hor, draw_hor_ort
from .label import draw_label
from .self_ident import draw_self_ident
from .bar import draw_bars
from .legend import draw_legend
from .utils import create_subplots, format_ax, set_both_labels
from ..io.utils import get_min_max_track
from ..track.types import Track, TrackType, TrackPosition, LegendPosition


def plot_one_cen(
    tracks: list[Track],
    outdir: str,
    chrom: str,
    settings: PlotSettings,
) -> tuple[Figure, np.ndarray, list[str]]:
    """
    Plot a single centromere figure from a list of `Track`s.

    # Args
    * `tracks`
        * List of tracks to plot. The order in the list determines placement on the figure.
    * `outdir`
        * Output directory.
    * `chrom`
        * Chromosome name to filter for in `Track.data`
    * `settings`
        * Settings for output plots.

    # Returns
    * Figure, its axes, and the output filename(s).

    # Usage
    ```python
    import cenplot

    chrom = "chm13_chr10:38568472-42561808"
    track_list, settings = cenplot.read_one_cen_tracks("tracks_example_api.toml", chrom=chrom)
    fig, axes, outfiles = cenplot.plot_one_cen(track_list.tracks, "plots", chrom, settings)
    ```
    """
    # Show chrom trimmed of spaces for logs and filenames.
    logging.info(f"Plotting {chrom}...")

    if not settings.xlim:
        # Get min and max position of all tracks for this cen.
        _, min_st_pos = get_min_max_track(tracks, typ="min")
        _, max_end_pos = get_min_max_track(tracks, typ="max", default_col="chrom_end")
    else:
        min_st_pos = settings.xlim[0]
        max_end_pos = settings.xlim[1]

    # # Scale height based on track length.
    # adj_height = height * (trk_max_end / max_end_pos)
    # height = height if adj_height == 0 else adj_height

    fig, axes, track_indices = create_subplots(
        tracks,
        settings,
    )
    if settings.legend_pos == LegendPosition.Left:
        track_col, legend_col = 1, 0
    else:
        track_col, legend_col = 0, 1

    track_labels: list[str] = []

    def get_track_label(chrom: str, track: Track, all_track_labels: list[str]) -> str:
        if not track.title:
            return ""
        try:
            fmt_track_label = track.title.format(chrom=chrom)
        except KeyError:
            fmt_track_label = track.title

        track_label = fmt_track_label.encode("ascii", "ignore").decode("unicode_escape")

        # Update track label for each overlap.
        if track.pos == TrackPosition.Overlap:
            try:
                track_label = f"{all_track_labels[-1]}\n{track_label}"
            except IndexError:
                pass

        return track_label

    num_hor_split = 0
    for idx, track in enumerate(tracks):
        track_row = track_indices[idx]
        track_label = get_track_label(chrom, track, track_labels)

        try:
            track_ax: Axes = axes[track_row, track_col]
        except IndexError:
            print(f"Cannot get track ({track_row, track_col}) for {track}.")
            continue
        try:
            legend_ax: Axes = axes[track_row, legend_col]
        except IndexError:
            legend_ax = None

        # Set xaxis limits
        track_ax.set_xlim(min_st_pos, max_end_pos)

        if track.opt == TrackType.Legend:
            draw_legend(track_ax, axes, track, tracks, track_row, track_col)
        elif track.opt == TrackType.Position:
            # Hide everything but x-axis
            format_ax(
                track_ax,
                grid=True,
                yticks=True,
                spines=("right", "left", "top"),
            )
        elif track.opt == TrackType.Spacer:
            # Hide everything.
            format_ax(
                track_ax,
                grid=True,
                xticks=True,
                yticks=True,
                spines=("right", "left", "top", "bottom"),
            )
        else:
            # Switch track option. {bar, label, ident, hor}
            # Add legend.
            if track.opt == TrackType.HOR or track.opt == TrackType.HORSplit:
                draw_fn = draw_hor
            elif track.opt == TrackType.HOROrt:
                draw_fn = draw_hor_ort
            elif track.opt == TrackType.Label:
                draw_fn = draw_label
            elif track.opt == TrackType.SelfIdent:
                draw_fn = draw_self_ident
            elif track.opt == TrackType.Bar:
                draw_fn = draw_bars

            else:
                raise ValueError("Invalid TrackType. Unreachable.")

            draw_fn(
                ax=track_ax,
                legend_ax=legend_ax,
                track=track,
                zorder=idx,
            )

        # Store label if more overlaps.
        track_labels.append(track_label)

        # Set labels for both x and y axis.
        set_both_labels(track_label, track_ax, track)

        if not legend_ax:
            continue

        # Make legend title invisible for HORs split after 1.
        if track.opt == TrackType.HORSplit:
            legend_ax_legend = legend_ax.get_legend()
            if legend_ax_legend and num_hor_split != 0:
                legend_title = legend_ax_legend.get_title()
                legend_title.set_alpha(0.0)
            num_hor_split += 1

        # Minimalize all legend cols except self-ident
        if track.opt != TrackType.SelfIdent or (
            track.opt == TrackType.SelfIdent and not track.options.legend
        ):
            format_ax(
                legend_ax,
                grid=True,
                xticks=True,
                xticklabel_fontsize=track.options.legend_fontsize,
                yticks=True,
                yticklabel_fontsize=track.options.legend_fontsize,
                spines=("right", "left", "top", "bottom"),
            )
        else:
            format_ax(
                legend_ax,
                grid=True,
                xticklabel_fontsize=track.options.legend_fontsize,
                yticklabel_fontsize=track.options.legend_fontsize,
                spines=("right", "top"),
            )

    # Add title
    if settings.title:
        title = settings.title.format(chrom=chrom)
        fig.suptitle(
            title, x=0.02, y=0.98, horizontalalignment="left", fontsize="xx-large"
        )

    os.makedirs(outdir, exist_ok=True)
    if isinstance(settings.format, str):
        output_format = [settings.format]
    else:
        output_format = settings.format
    # Pad between axes.
    fig.get_layout_engine().set(h_pad=settings.axis_h_pad)

    # PNG must always be plotted last.
    # Matplotlib modifies figure settings causing formatting errors in vectorized image formats (svg, pdf)
    png_output = "png" in output_format
    if png_output:
        output_format.remove("png")

    outfiles = []
    for fmt in output_format:
        outfile = os.path.join(outdir, f"{chrom}.{fmt}")
        fig.savefig(outfile, dpi=settings.dpi, transparent=settings.transparent)
        outfiles.append(outfile)

    if png_output:
        outfile = os.path.join(outdir, f"{chrom}.png")
        fig.savefig(
            outfile,
            dpi=settings.dpi,
            transparent=settings.transparent,
        )
        outfiles.append(outfile)

    return fig, axes, outfiles
