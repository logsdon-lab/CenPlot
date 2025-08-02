import os
import logging
import polars as pl

from typing import Any, TextIO

from .bed9 import read_bed9
from .utils import map_value_colors
from ..defaults import MONOMER_COLORS, BED9_COLS
from ..track.settings import HORTrackSettings


def read_bed_hor(
    infile: str | TextIO,
    *,
    chrom: str | None = None,
    live_only: bool = True,
    mer_size: int = HORTrackSettings.mer_size,
    mer_filter: int = HORTrackSettings.mer_filter,
    hor_filter: int | None = None,
    sort_by: str = "mer",
    sort_order: str = HORTrackSettings.sort_order,
    sort_fill_missing: str | None = HORTrackSettings.split_fill_missing,
    sort_order_only: bool = False,
    color_map_file: str | None = None,
    use_item_rgb: bool = HORTrackSettings.use_item_rgb,
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
    * `mer_size`
        * Monomer size to calculate monomer number.
    * `mer_filter`
        * Filter for HORs with at least this many monomers.
    * `hor_filter`
        * Filter for HORs that occur at least this many times.
    * `color_map_file`
        * Convenience color map file for `mer` or `hor`.
        * Two-column TSV file with no header.
        * If `None`, use default color map.
    * `sort_by`
        * Sort `pl.DataFrame` by `mer`, `hor`, or `hor_count`.
        * Can be a path to a list of `mer` or `hor` names
    * `sort_order`
        * Sort in ascending or descending order.
    * `sort_fill_missing`
        * Fill in missing elements in defined sort order with this color.
    * `sort_order_only`
        * Convenience switch to keep only elements in defined sort order.
    * `use_item_rgb`
        * Use `item_rgb` column or generate random colors.

    # Returns
    * HOR `pl.DataFrame`
    """
    df = read_bed9(infile, chrom=chrom)

    if df.is_empty():
        return pl.DataFrame(schema=[*BED9_COLS, "mer", "length", "color", "hor_count"])

    df = (
        df.lazy()
        .with_columns(
            length=pl.col("chrom_end") - pl.col("chrom_st"),
        )
        .with_columns(
            mer=(pl.col("length") / mer_size).round().cast(pl.Int8).clip(1, 100)
        )
        .filter(
            pl.when(live_only).then(pl.col("name").str.contains("L")).otherwise(True)
            & (pl.col("mer") >= mer_filter)
        )
        .collect()
    )
    # Read color map.
    if color_map_file:
        color_map: dict[str, str] = {}
        with open(color_map_file, "rt") as fh:
            for line in fh.readlines():
                try:
                    name, color = line.strip().split()
                except Exception:
                    logging.error(f"Invalid color map. ({line})")
                    continue
                color_map[name] = color
    else:
        color_map = MONOMER_COLORS

    df = map_value_colors(
        df,
        map_col="mer",
        map_values=MONOMER_COLORS,
        use_item_rgb=use_item_rgb,
    )
    df = df.join(df.get_column("name").value_counts(name="hor_count"), on="name")

    if hor_filter:
        df = df.filter(pl.col("hor_count") >= hor_filter)

    if os.path.exists(sort_order):
        with open(sort_order, "rt") as fh:
            defined_sort_order = []
            for line in fh:
                line = line.strip()
                defined_sort_order.append(int(line) if sort_by == "mer" else line)
    else:
        defined_sort_order = None

    if sort_by == "mer":
        sort_col = "mer"
    elif sort_by == "name" and defined_sort_order:
        sort_col = "name"
    else:
        sort_col = "hor_count"

    if defined_sort_order:
        # Add missing elems in df not in sort order so all elements covered.
        all_elems = [
            *defined_sort_order,
            *set(df[sort_col]).difference(defined_sort_order),
        ]
        # Missing elements in sort order not in df
        missing_elems = set(defined_sort_order).difference(df[sort_col])

        # Fill in missing.
        if sort_fill_missing and missing_elems:
            row_template = df.row(0, named=True)
            min_st, max_end = df["chrom_st"].min(), df["chrom_end"].max()
            df_missing_element_rows = pl.DataFrame(
                [
                    {
                        **row_template,
                        "chrom_st": min_st,
                        "chrom_end": max_end,
                        "strand": ".",
                        "thick_st": min_st,
                        "thick_end": max_end,
                        "name": elem,
                        "mer": 0,
                        "length": 0,
                        "hor_count": 0,
                        "item_rgb": sort_fill_missing,
                        "color": sort_fill_missing,
                    }
                    for elem in missing_elems
                ],
                schema=df.schema,
            )

            df = pl.concat([df, df_missing_element_rows])
        # Only take elements in sort order.
        if sort_order_only:
            df = df.filter(pl.col(sort_col).is_in(defined_sort_order))
            all_elems = defined_sort_order

        df = df.cast({sort_col: pl.Enum(all_elems)}).sort(by=sort_col)
    else:
        df = df.sort(sort_col, descending=sort_order == HORTrackSettings.sort_order)

    return df


def read_bed_hor_from_settings(
    path: str, options: dict[str, Any], chrom: str | None = None
) -> pl.DataFrame:
    live_only = options.get("live_only", HORTrackSettings.live_only)
    mer_filter = options.get("mer_filter", HORTrackSettings.mer_filter)
    hor_filter = options.get("hor_filter", HORTrackSettings.hor_filter)
    use_item_rgb = options.get("use_item_rgb", HORTrackSettings.use_item_rgb)
    sort_order = options.get("sort_order", HORTrackSettings.sort_order)
    color_map_file = options.get("color_map_file", HORTrackSettings.color_map_file)
    mer_size = options.get("mer_size", HORTrackSettings.mer_size)
    split_fill_missing = options.get(
        "split_fill_missing", HORTrackSettings.split_fill_missing
    )
    split_sort_order_only = options.get(
        "split_sort_order_only", HORTrackSettings.split_sort_order_only
    )

    if options.get("mode", HORTrackSettings.mode) == "hor":
        split_colname = "name"
    else:
        split_colname = "mer"

    df_hor = read_bed_hor(
        path,
        chrom=chrom,
        mer_size=mer_size,
        sort_by=split_colname,
        sort_order=sort_order,
        sort_fill_missing=split_fill_missing,
        sort_order_only=split_sort_order_only,
        live_only=live_only,
        mer_filter=mer_filter,
        hor_filter=hor_filter,
        use_item_rgb=use_item_rgb,
        color_map_file=color_map_file,
    )
    colnames = (
        df_hor.unique([split_colname, "hor_count"], maintain_order=True)
        .slice(0, options.get("split_top_n"))
        .get_column(split_colname)
    )
    return df_hor.filter(pl.col(split_colname).is_in(colnames))
