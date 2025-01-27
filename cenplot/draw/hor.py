import polars as pl
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, FancyArrowPatch

from .utils import draw_uniq_entry_legend, minimalize_ax
from ..defaults import MONOMER_COLORS, MONOMER_LEN


def draw_hor_ort(
    ax: Axes,
    df_stv_ort: pl.DataFrame,
    zorder: float,
    scale: float | None = None,
    fwd_color: str | None = None,
    rev_color: str | None = None,
):
    minimalize_ax(
        ax, xticks=True, yticks=True, spines=("right", "left", "top", "bottom")
    )

    if not fwd_color:
        fwd_color = "black"
    if not rev_color:
        rev_color = "black"
    if not scale:
        scale = 50

    for row in df_stv_ort.iter_rows(named=True):
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
            (start, 0),
            (end, 0),
            mutation_scale=scale,
            color=color,
            clip_on=False,
            zorder=zorder,
        )
        ax.add_patch(arrow)


def draw_hor(
    ax: Axes,
    df_stv: pl.DataFrame,
    *,
    hide_x: bool,
    zorder: float,
    legend_ax: Axes | None,
):
    minimalize_ax(
        ax, xticks=hide_x, yticks=True, spines=("right", "left", "top", "bottom")
    )

    # Add HOR track.
    for row in df_stv.with_row_index().iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        row_idx = row["index"]
        # TODO: Adjust zorder so mer order is respected.
        zorder_adj = zorder + (row_idx / df_stv.shape[0])
        color = MONOMER_COLORS.get(row["mer"])
        rect = Rectangle(
            (start, 0),
            end + 1 - start,
            1,
            color=color,
            lw=0,
            label=row["mer"],
            zorder=zorder_adj,
        )
        ax.add_patch(rect)

    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center left", ncols=3)
