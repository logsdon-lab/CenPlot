from matplotlib.axes import Axes
from matplotlib.patches import Rectangle

from cenplot.lib.draw.strand import draw_strand

from .utils import add_border, draw_uniq_entry_legend, format_ax
from ..track.types import Track, TrackPosition


def draw_hor_ort(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
):
    """
    Draw HOR ort plot on axis with the given `Track`.
    """
    draw_strand(ax, track, zorder=zorder, legend_ax=legend_ax)


def draw_hor(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
):
    """
    Draw HOR plot on axis with the given `Track`.
    """
    hide_x = track.options.hide_x
    legend = track.options.legend
    border = track.options.border

    if track.pos != TrackPosition.Overlap:
        spines = (
            ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
        )
    else:
        spines = None

    format_ax(
        ax,
        xticks=hide_x,
        xticklabel_fontsize=track.options.fontsize,
        yticks=True,
        yticklabel_fontsize=track.options.fontsize,
        spines=spines,
    )

    ylim = ax.get_ylim()
    height = ylim[1] - ylim[0]

    if track.options.mode == "hor":
        colname = "name"
    else:
        colname = "mer"

    # Add HOR track.
    for row in track.data.iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        color = row["color"]
        rect = Rectangle(
            (start, 0),
            end + 1 - start,
            height,
            color=color,
            lw=0,
            label=row[colname],
            zorder=zorder,
        )
        ax.add_patch(rect)

    if border:
        # Ensure border is always on top.
        add_border(ax, height, zorder + 1.0)

    if legend_ax and legend:
        draw_uniq_entry_legend(
            legend_ax,
            track,
            ref_ax=ax,
            ncols=track.options.legend_ncols,
            loc="center left",
            alignment="left",
        )
