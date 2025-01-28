import polars as pl

from matplotlib.axes import Axes
from matplotlib.collections import PolyCollection
from intervaltree import Interval, IntervalTree

from .utils import minimalize_ax
from ..track import Track
from ..defaults import IDENT_COLOR_RANGE


def draw_self_ident(
    ax: Axes,
    track: Track,
    *,
    zorder: float,
    hide_x: bool,
    flip_y: bool,
    bins: int = 300,
    legend_ax: Axes | None = None,
    legend_xmin: float | None = None,
    legend_aspect_ratio: float | None = None,
) -> None:
    colors, verts = [], []
    spines = ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
    minimalize_ax(ax, xticks=hide_x, yticks=True, spines=spines)

    if not legend_xmin:
        legend_xmin = 0.7

    if flip_y:
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

    if legend_ax:
        cmap = IntervalTree(
            Interval(rng[0], rng[1], color) for rng, color in IDENT_COLOR_RANGE.items()
        )
        cnts, values, bars = legend_ax.hist(
            track.data["percent_identity_by_events"], bins=bins, zorder=zorder
        )
        legend_ax.set_xlim(legend_xmin, 100.0)
        legend_ax.minorticks_on()
        legend_ax.set_xlabel("Mean nucleotide identity\nbetween pairwise intervals")
        legend_ax.set_ylabel("# of Intervals (thousands)")

        # Ensure that legend is only a portion of the total height.
        # Otherwise, take up entire axis dim.
        legend_ax.set_box_aspect(legend_aspect_ratio)

        for _, value, bar in zip(cnts, values, bars):
            # Make value a non-null interval
            # ex. (1,1) -> (1, 1.000001)
            color = cmap.overlap(value, value + 0.00001)
            try:
                color = next(iter(color)).data
            except Exception:
                color = None
            bar.set_facecolor(color)
