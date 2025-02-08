import polars as pl

from typing import TextIO

from .bed9 import read_bed9
from .utils import map_value_colors
from ..defaults import MONOMER_COLORS, MONOMER_LEN
from ..track.settings import HORPlotSettings


def read_bed_hor(
    infile: str | TextIO,
    *,
    chrom: str | None = None,
    live_only: bool = True,
    mer_filter: int = HORPlotSettings.mer_filter,
    hor_filter: int | None = None,
    sort_col: str = "mer",
    sort_order: str = HORPlotSettings.sort_order,
    use_item_rgb: bool = HORPlotSettings.use_item_rgb,
) -> pl.DataFrame:
    """
    Read a HOR BED9 file with no header.

    # Args
    * `infile`
        * Input file or IO stream.
    * `chrom`
        * Chromsome in `chrom` column to filter for.
    * `live_only`
        * Filter for only live data.
        * Contains `L` in `name` column.
    * `mer_filter`
        * Filter for HORs with at least this many monomers.
    * `hor_filter`
        * Filter for HORs that occur at least this many times.
    * `sort_col`
        * Sort `pl.DataFrame` by `mer` or `hor_count`.
    * `sort_order`
        * Sort in ascending or descending order.
    * `use_item_rgb`
        * Use `item_rgb` column or generate random colors.

    # Returns
    * HOR `pl.DataFrame`
    """
    df = (
        read_bed9(infile, chrom=chrom)
        .lazy()
        .with_columns(
            length=pl.col("chrom_end") - pl.col("chrom_st"),
        )
        .with_columns(
            mer=(pl.col("length") / MONOMER_LEN).round().cast(pl.Int8).clip(1, 100)
        )
        .filter(
            pl.when(live_only).then(pl.col("name").str.contains("L")).otherwise(True)
            & (pl.col("mer") >= mer_filter)
        )
        .collect()
    )

    df = map_value_colors(
        df,
        map_col="mer",
        map_values=MONOMER_COLORS,
        use_item_rgb=use_item_rgb,
    )
    df = df.join(df.get_column("name").value_counts(name="hor_count"), on="name")

    if hor_filter:
        df = df.filter(pl.col("hor_count") >= hor_filter)

    if sort_col == "mer":
        df = df.sort(sort_col, descending=sort_order == HORPlotSettings.sort_order)
    else:
        df = df.sort("hor_count", descending=sort_order == HORPlotSettings.sort_order)

    return df
