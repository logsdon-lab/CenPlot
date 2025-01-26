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
    legend_ax: Axes | None = None,
    zorder: float,
    hide_x: bool,
    flip_y: bool,
    bins: int = 300,
) -> None:
    colors, verts = [], []

    minimalize_ax(
        ax, xticks=hide_x, yticks=True, spines=("left", "right", "top", "bottom")
    )

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
        legend_ax.set_xlim(70.0, 100.0)
        legend_ax.set_xlabel("Average Nucleotide Identity (%)")
        legend_ax.set_ylabel("Number of Intervals (k)")

        for _, value, bar in zip(cnts, values, bars):
            color = cmap.overlap(value, value + 0.00001)
            try:
                color = next(iter(color)).data
            except Exception:
                color = None
            bar.set_facecolor(color)
