import polars as pl

from typing import TextIO

from .utils import adj_by_ctg_coords, header_info
from ..defaults import BED9_COLS


def read_bed9(infile: str | TextIO, *, chrom: str | None = None) -> pl.DataFrame:
    """
    Read a BED9 file with no header.

    # Args
    * `infile`
        * Input file or IO stream.
    * `chrom`
        * Chromsome in `chrom` column to filter for. If contains coordinates, subset to those coordinates.

    # Returns
    * BED9 pl.DataFrame.
    """
    skip_rows, number_cols = header_info(infile)

    try:
        df = pl.scan_csv(
            infile,
            separator="\t",
            has_header=False,
            skip_rows=skip_rows,
            new_columns=BED9_COLS[0:number_cols],
        )
        try:
            chrom_no_coords, coords = chrom.rsplit(":", 1)
            chrom_st, chrom_end = [int(elem) for elem in coords.split("-")]
        except Exception:
            chrom_no_coords = None
            chrom_st, chrom_end = None, None

        if chrom:
            df = df.filter(
                pl.when(pl.col("chrom").is_in([chrom_no_coords]))
                .then(
                    (pl.col("chrom") == chrom_no_coords)
                    & (pl.col("chrom_st").is_between(chrom_st, chrom_end))
                    & (pl.col("chrom_end").is_between(chrom_st, chrom_end))
                )
                .when(pl.col("chrom").is_in([chrom]))
                .then(pl.col("chrom") == chrom)
                .otherwise(True)
            ).collect()

        df_adj = adj_by_ctg_coords(df, "chrom").sort(by="chrom_st")
    except pl.exceptions.NoDataError:
        df_adj = pl.DataFrame(schema=BED9_COLS)

    if "item_rgb" not in df_adj.columns:
        df_adj = df_adj.with_columns(item_rgb=pl.lit("0,0,0"))
    if "name" not in df_adj.columns:
        df_adj = df_adj.with_columns(name=pl.lit("-"))

    return df_adj
