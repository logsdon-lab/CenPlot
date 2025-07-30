import csv
import polars as pl

from typing import TextIO

from .utils import adj_by_ctg_coords
from ..defaults import BED9_COL_MAP


def read_bed9(infile: str | TextIO, *, chrom: str | None = None) -> pl.DataFrame:
    """
    Read a BED9 file with no header.

    # Args
    * `infile`
        * Input file or IO stream.
    * `chrom`
        * Chromsome in `chrom` column to filter for.

    # Returns
    * BED9 pl.DataFrame.
    """
    sniffer = csv.Sniffer()
    fname = infile if isinstance(infile, str) else infile.name
    with open(fname, "rt") as fh:
        header = next(fh)

    skip_rows = 1 if sniffer.has_header(header) else 0

    try:
        df = pl.read_csv(infile, separator="\t", has_header=False, skip_rows=skip_rows)
        df = df.rename(
            {col: val for col, val in BED9_COL_MAP.items() if col in df.columns}
        )
        df_adj = adj_by_ctg_coords(df, "chrom").sort(by="chrom_st")
    except Exception:
        df_adj = pl.DataFrame(schema=BED9_COL_MAP.values())

    if chrom:
        df_adj = df_adj.filter(pl.col("chrom") == chrom)
    if "item_rgb" not in df_adj.columns:
        df_adj = df_adj.with_columns(item_rgb=pl.lit("0,0,0"))
    if "name" not in df_adj.columns:
        df_adj = df_adj.with_columns(name=pl.lit("-"))

    return df_adj
