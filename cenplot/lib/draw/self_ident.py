import polars as pl

from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from intervaltree import Interval, IntervalTree

from .utils import format_ax
from ..track.types import Track
from ..defaults import IDENT_COLOR_RANGE


def draw_self_ident_hist(ax: Axes, track: Track, *, zorder: float = 1.0):
    """
    Draw self identity histogram plot on axis with the given `Track`.
    """
    legend_bins = track.options.legend_bins
    legend_xmin = track.options.legend_xmin
    legend_asp_ratio = track.options.legend_asp_ratio

    cmap = IntervalTree(
        Interval(rng[0], rng[1], color) for rng, color in IDENT_COLOR_RANGE.items()
    )
    cnts, values, bars = ax.hist(
        track.data["percent_identity_by_events"], bins=legend_bins, zorder=zorder
    )
    ax.set_xlim(legend_xmin, 100.0)
    ax.minorticks_on()
    ax.set_xlabel(
        "Mean nucleotide identity\nbetween pairwise intervals",
        fontsize=track.options.legend_title_fontsize,
    )
    ax.set_ylabel(
        "# of Intervals (thousands)", fontsize=track.options.legend_title_fontsize
    )

    # Ensure that legend is only a portion of the total height.
    # Otherwise, take up entire axis dim.
    ax.set_box_aspect(legend_asp_ratio)

    for _, value, bar in zip(cnts, values, bars):
        # Make value a non-null interval
        # ex. (1,1) -> (1, 1.000001)
        color = cmap.overlap(value, value + 0.00001)
        try:
            color = next(iter(color)).data
        except Exception:
            color = None
        bar.set_facecolor(color)


def draw_self_ident(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
) -> None:
    """
    Draw self identity plot on axis with the given `Track`.
    """
    hide_x = track.options.hide_x
    invert = track.options.invert
    legend = track.options.legend

    colors, verts = [], []
    spines = ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
    format_ax(
        ax,
        xticks=hide_x,
        xticklabel_fontsize=track.options.fontsize,
        yticks=True,
        yticklabel_fontsize=track.options.fontsize,
        spines=spines,
    )

    if invert:
        df_track = track.data.with_columns(y=-pl.col("y"))
    else:
        df_track = track.data

    for _, df_diam in df_track.group_by(["group"]):
        df_points = df_diam.select("x", "y")
        color = df_diam["color"].first()
        colors.append(color)
        verts.append(df_points)

    # https://stackoverflow.com/a/29000246
    polys = PolyCollection(verts, zorder=zorder)
    polys.set(array=None, facecolors=colors)
    ax.add_collection(polys)

    ax.set_ylim(df_track["y"].min(), df_track["y"].max())

    if legend_ax and legend:
        draw_self_ident_hist(legend_ax, track, zorder=zorder)
