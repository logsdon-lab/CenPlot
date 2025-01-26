import polars as pl
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, FancyArrow

from .utils import draw_uniq_entry_legend, minimalize_ax
from ..defaults import MONOMER_COLORS


def draw_hor_w_ort(
    ax: Axes,
    df_stv: pl.DataFrame,
    df_stv_ort: pl.DataFrame,
    *,
    hide_x: bool,
    ort: bool,
    ort_pos: str,
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

    # Add ort track.
    if ort_pos == "top":
        y_ort = 1.2
    else:
        y_ort = -0.1
    if ort:
        for row in df_stv_ort.iter_rows(named=True):
            # sample arrow
            start = row["chrom_st"]
            end = row["chrom_end"]
            strand = row["strand"]
            length = end - start
            # Skip single orts
            if length < 171:
                continue

            if strand == "-":
                start = end
                length = -length

            arrow = FancyArrow(
                start,
                y_ort,
                length,
                0,
                width=0.05,
                color="black",
                length_includes_head=True,
                # 10% of entire length
                head_length=0.05 * -length,
                clip_on=False,
                zorder=zorder,
            )
            ax.add_patch(arrow)

    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center right", ncols=3)
