import numpy as np
import matplotlib.pyplot as plt

from typing import Any, Iterable, NamedTuple

from matplotlib.axes import Axes
from matplotlib.artist import Artist
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages

from .settings import PlotSettings
from ..utils import Unit
from ..track.types import LegendPosition, Track, TrackType, TrackPosition
from ..track.settings import DefaultTrackSettings


class ComparePlotAxCoords(NamedTuple):
    ref_coords: dict[int, tuple[int, int]]
    ref_legend_coords: dict[int, tuple[int, int]]
    qry_coords: dict[int, tuple[int, int]]
    qry_legend_coords: dict[int, tuple[int, int]]
    ident_coords: tuple[int, int]
    ident_hist_coords: tuple[int, int]


class Delta(NamedTuple):
    coord: tuple[int, int]
    dt: tuple[int, int]

    def adjust_by(self, adj: int) -> tuple[int, int]:
        return (self.coord[0] + self.dt[0] * adj, self.coord[1] + self.dt[1] * adj)


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
        squeeze=False,
        layout=settings.layout,
        **kwargs,
    )

    return fig, axes, track_indices


def create_cmp_subplots(
    tracks_ref: list[Track],
    tracks_qry: list[Track],
    settings_ref: PlotSettings,
    settings_qry: PlotSettings,
    pos_ref: str,
    pos_qry: str,
    prop_ident: float,
    **kwargs: Any,
) -> tuple[Figure, np.ndarray, ComparePlotAxCoords]:
    """
    Generate a subplot figure from ref and query `Track`s.

    # Args
    * `tracks_ref`
    * `tracks_qry`
    * `settings_ref`
        * Plot settings.
        * Legend pos is ignored.
    * `settings_qry`
        * Plot settings.
        * Legend pos is ignored.
    * `pos_qry`
    * `tracks_qry`
    * `prop_ident`
    * `kwargs`
        * Additional arguments passed to `plt.subplots`

    # Returns
    * Figure, its axes, and a mapping of the idx of `tracks` to the row idx of axes.
    """

    def get_track_ordering(
        tracks: list[Track],
    ) -> tuple[list[float], dict[int, int]]:
        track_idx = 0
        track_props = []
        track_indices = {}
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

        return track_props, track_indices

    ref_track_props, ref_track_indices = get_track_ordering(tracks_ref)
    qry_track_props, qry_track_indices = get_track_ordering(tracks_qry)

    # Determine correct proportions for entire figure.
    ref_legend_prop = settings_ref.legend_prop
    qry_legend_prop = settings_qry.legend_prop

    # 1. Ref top.     2. Ref bottom.
    #    Query left.     Query left. (DEFAULT)
    #      +----+-+      +-+
    #    +-+----+-+      +-+----+
    #    | | \  |        | | \  |
    #    | |  \ |h       | |  \ |h
    #    +-+----+        +-+----+-+
    #    +-+               +----+-+
    #
    # 3. Ref top.     4. Ref bottom.
    #    Query right.    Query right.
    #     +-+----+              +-+
    #     +-+----+-+       +----+-+
    #       | \  | |       | \  | |
    #      h|  \ | |      h|  \ | |
    #       +----+-+     +-+----+-+
    #            +-+     +-+----+
    if pos_ref == "top" and pos_qry == "left":
        ident_coord = (len(ref_track_props), len(qry_track_props))
        ident_hist_coord = (ident_coord[0], ident_coord[1] + 1)
        # Delta for annot row and legend from origin scaled by number of tracks.
        ref_delta = (Delta(ident_coord, (-1, 0)), Delta(ident_hist_coord, (-1, 0)))
        qry_delta = (
            Delta(ident_coord, (0, -1)),
            Delta((ident_coord[0] + 1, ident_coord[1]), (1, 0)),
        )
        height_ratios = [*ref_track_props, prop_ident, qry_legend_prop]
        width_ratios = [*qry_track_props, prop_ident, ref_legend_prop]
    elif pos_ref == "bottom" and pos_qry == "left":
        ident_coord = (1, len(qry_track_props))
        ident_hist_coord = (1, ident_coord[1] + 1)
        ref_delta = (Delta(ident_coord, (1, 0)), Delta(ident_hist_coord, (1, 0)))
        qry_delta = (
            Delta(ident_coord, (0, -1)),
            Delta((ident_coord[0] - 1, ident_coord[1]), (0, -1)),
        )
        height_ratios = [qry_legend_prop, prop_ident, *ref_track_props]
        width_ratios = [
            *qry_track_props,
            prop_ident,
            ref_legend_prop,
        ]
    elif pos_ref == "top" and pos_qry == "right":
        ident_coord = (len(ref_track_props), 1)
        ident_hist_coord = (ident_coord[0], 0)
        ref_delta = (Delta(ident_coord, (-1, 0)), Delta(ident_hist_coord, (-1, 0)))
        qry_delta = (
            Delta(ident_coord, (0, 1)),
            Delta((ident_coord[0] + 1, ident_coord[1]), (0, 1)),
        )
        height_ratios = [*ref_track_props, prop_ident, qry_legend_prop]
        width_ratios = [ref_legend_prop, prop_ident, *qry_track_props]
    elif pos_ref == "bottom" and pos_qry == "right":
        ident_coord = (1, 1)
        ident_hist_coord = (ident_coord[0], 0)
        ref_delta = (Delta(ident_coord, (1, 0)), Delta(ident_hist_coord, (1, 0)))
        qry_delta = (
            Delta(ident_coord, (0, 1)),
            Delta((ident_coord[0] - 1, ident_coord[1]), (0, 1)),
        )
        height_ratios = [qry_legend_prop, prop_ident, *ref_track_props]
        width_ratios = [ref_legend_prop, prop_ident, *qry_track_props]
    else:
        raise NotImplementedError()

    def get_new_coords(
        track_indices: dict[int, int], deltas: tuple[Delta, Delta]
    ) -> tuple[dict[int, tuple[int, int]], dict[int, tuple[int, int]]]:
        track_coords = {}
        track_legend_coords = {}
        for track_idx, track_pos in track_indices.items():
            dt_track, dt_track_legend = deltas
            # Ident track is first position so offset.
            coords_track = dt_track.adjust_by(track_pos + 1)
            coords_track_legend = dt_track_legend.adjust_by(track_pos + 1)
            track_coords[track_idx] = coords_track
            track_legend_coords[track_idx] = coords_track_legend
        return track_coords, track_legend_coords

    ref_track_coords, ref_track_legend_coords = get_new_coords(
        ref_track_indices, ref_delta
    )
    qry_track_coords, qry_track_legend_coords = get_new_coords(
        qry_track_indices, qry_delta
    )

    coords_ax = ComparePlotAxCoords(
        ref_coords=ref_track_coords,
        ref_legend_coords=ref_track_legend_coords,
        qry_coords=qry_track_coords,
        qry_legend_coords=qry_track_legend_coords,
        ident_coords=ident_coord,
        ident_hist_coords=ident_hist_coord,
    )

    fig, axes = plt.subplots(
        # Count number of tracks
        nrows=len(height_ratios),
        ncols=len(width_ratios),
        figsize=settings_ref.dim,
        height_ratios=height_ratios,
        width_ratios=width_ratios,
        layout=settings_ref.layout,
        **kwargs,
    )

    return fig, axes, coords_ax


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
    new_xtick_labels = []
    units = Unit(track.options.units_x)
    xmin, xmax = ax.get_xlim()
    xticks, xticklabels = list(ax.get_xticks()), ax.get_xticklabels()
    xticklabels_txt = set([lbl.get_text() for lbl in xticklabels])
    ax.xaxis.set_major_formatter("plain")

    for txt in xticklabels:
        x, _ = txt.get_position()
        # Convert units and round.
        new_x_txt = units.convert_value(x, 3)
        txt.set_text(new_x_txt)
        new_xtick_labels.append(txt)

    # Add first and last position.
    # Check if txt exists before adding.
    xmin_txt = str(units.convert_value(xmin, 1))
    xmax_txt = str(units.convert_value(xmax, 1))

    if xmin_txt not in xticklabels_txt:
        xticks.append(xmin)
        new_xtick_labels.append(xmin_txt)  # type: ignore[arg-type]
    if xmax_txt not in xticklabels_txt:
        new_xtick_labels.append(xmax_txt)  # type: ignore[arg-type]
        xticks.append(xmax)

    ax.set_xticks(xticks, new_xtick_labels, fontsize=track.options.fontsize)  # type: ignore[arg-type]
    ax.set_xlabel(
        f"Position ({units.capitalize()})", fontsize=track.options.title_fontsize
    )
    ax.set_xlim(xmin, xmax)


def draw_uniq_entry_legend(
    ax: Axes,
    track: Track,
    ref_ax: Axes | None = None,
    ncols: int | None = DefaultTrackSettings.legend_ncols,
    label_order: Iterable[str] | None = None,
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
    # Order by appearance or set order.
    handles, labels = ref_ax.get_legend_handles_labels()
    by_label: dict[Any, Artist] = dict(zip(labels, handles))
    if label_order:
        by_label = {
            label: by_label[label] for label in label_order if by_label.get(label)
        }
    else:
        by_label = dict(zip(labels, handles))

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
            fontsize=track.options.legend_fontsize,
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


def add_rect(
    ax: Axes,
    height: float,
    zorder: float,
    color: str | None = None,
    fill: bool = False,
    edgecolor: str | None = None,
):
    xmin, xmax = ax.get_xlim()
    rect = Rectangle(
        (xmin, 0),
        xmax - xmin,
        height,
        edgecolor=edgecolor,
        fill=fill,
        color=color,
        zorder=zorder,
    )
    ax.add_patch(rect)


def set_ylim(ax: Axes, track: Track) -> None:
    ylim_args = {}
    if hasattr(track.options, "ymin"):
        ymin = track.options.ymin
        ymin_add = track.options.ymin_add
        ymin_val = track.data["name"].min()
        if isinstance(ymin_val, (int, float)) and isinstance(ymin_add, float):
            if ymin == "min":
                ylim_args["ymin"] = ymin_val + (ymin_val * ymin_add)
            elif isinstance(ymin, (int, float)):
                ylim_args["ymin"] = ymin + (ymin_val * ymin_add)

    if hasattr(track.options, "ymax"):
        ymax = track.options.ymax
        ymax_add = track.options.ymax_add
        ymax_val = track.data["name"].max()
        if isinstance(ymax_val, (int, float)) and isinstance(ymax_add, float):
            if ymax == "max":
                ylim_args["ymax"] = ymax_val + (ymax_val * ymax_add)
            elif isinstance(ymax, (int, float)):
                ylim_args["ymax"] = ymax + (ymax_val * ymax_add)

    if ylim_args:
        ax.set_ylim(**ylim_args)  # type: ignore[arg-type]

    if track.options.add_end_yticks:
        ymin, ymax = ax.get_ylim()
        yticks, yticklabels = list(ax.get_yticks()), ax.get_yticklabels()

        # Add first and last position.
        if ymax not in yticks:
            yticks.append(ymax)
            yticklabels.append(str(round(ymax, 1)))  # type: ignore[arg-type]
        if ymin not in yticks:
            yticks.append(ymin)
            yticklabels.append(str(round(ymin, 1)))  # type: ignore[arg-type]

        ax.set_yticks(yticks, yticklabels, fontsize=track.options.fontsize)  # type: ignore[arg-type]
        ax.set_ylim(ymin, ymax)
