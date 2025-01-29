from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, FancyArrowPatch

from .utils import draw_uniq_entry_legend, format_ax
from ..track.types import Track
from ..defaults import MONOMER_COLORS, MONOMER_LEN


def draw_hor_ort(
    ax: Axes,
    track: Track,
    *,
    zorder: float,
    legend_ax: Axes | None = None,
):
    hide_x = track.options.hide_x
    fwd_color = (
        track.options.fwd_color if track.options.fwd_color else track.options.DEF_COLOR
    )
    rev_color = (
        track.options.rev_color if track.options.rev_color else track.options.DEF_COLOR
    )
    scale = track.options.scale
    legend = track.options.legend

    spines = ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
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

    for row in track.data.iter_rows(named=True):
        # sample arrow
        start = row["chrom_st"]
        end = row["chrom_end"]
        strand = row["strand"]
        length = end - start
        # Skip single orts
        if length < MONOMER_LEN:
            continue

        if strand == "-":
            tmp_start = start
            start = end
            end = tmp_start
            color = rev_color
        else:
            color = fwd_color

        arrow = FancyArrowPatch(
            (start, height * 0.5),
            (end, height * 0.5),
            mutation_scale=scale,
            color=color,
            clip_on=False,
            zorder=zorder,
        )
        ax.add_patch(arrow)

    if legend_ax and legend:
        draw_uniq_entry_legend(legend_ax, track, ref_ax=ax, loc="center left", ncols=3)


def draw_hor(
    ax: Axes,
    track: Track,
    *,
    zorder: float,
    legend_ax: Axes | None = None,
):
    hide_x = track.options.hide_x
    legend = track.options.legend

    spines = ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
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

    # Add HOR track.
    for row in track.data.with_row_index().iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        row_idx = row["index"]
        # TODO: Adjust zorder so mer order is respected.
        zorder_adj = zorder + (row_idx / track.data.shape[0])
        color = MONOMER_COLORS.get(row["mer"])
        rect = Rectangle(
            (start, 0),
            end + 1 - start,
            height,
            color=color,
            lw=0,
            label=row["mer"],
            zorder=zorder_adj,
        )
        ax.add_patch(rect)

    if legend_ax and legend:
        draw_uniq_entry_legend(legend_ax, track, ref_ax=ax, loc="center left", ncols=3)
