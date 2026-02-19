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

        def expr_chrom_coords(
            expr_no_coords: pl.Expr, expr_coords: pl.Expr, expr_otherwise: pl.Expr
        ) -> pl.Expr:
            return (
                pl.when(pl.col("chrom").eq(chrom_no_coords))
                .then(expr_no_coords)
                .when(pl.col("chrom").eq(chrom))
                .then(expr_coords)
                .otherwise(expr_otherwise)
            )

        # Chrom coordinates can be one of three states:
        # 1. chr1:0-10:0-5
        # 2. chr1:0-10
        # 3. chr1
        # We assume if an exact match for chrom_no_coords is found (1), the user wants to trim to some coordinates.
        # Coordinate are right split once.
        if chrom_no_coords and chrom_st and chrom_end:
            df_filtered = (
                df.filter(
                    expr_chrom_coords(
                        pl.col("chrom") == chrom_no_coords,
                        pl.col("chrom") == chrom,
                        True,
                    )
                )
                .with_columns(
                    chrom_st=expr_chrom_coords(
                        pl.col("chrom_st").clip(chrom_st, chrom_end),
                        pl.col("chrom_st"),
                        pl.col("chrom_st"),
                    ),
                    chrom_end=expr_chrom_coords(
                        pl.col("chrom_end").clip(chrom_st, chrom_end),
                        pl.col("chrom_end"),
                        pl.col("chrom_end"),
                    ),
                )
                # Remove null intervals created by clipping to boundaries
                .filter(
                    expr_chrom_coords(
                        ~(
                            (
                                pl.col("chrom_st").eq(chrom_st)
                                & pl.col("chrom_st").eq(chrom_end)
                            )
                            | (
                                pl.col("chrom_end").eq(chrom_st)
                                & pl.col("chrom_end").eq(chrom_end)
                            )
                        ),
                        True,
                        True,
                    )
                )
                .collect()
            )
        elif chrom:
            df_filtered = df.filter(pl.col("chrom") == chrom).collect()
        else:
            df_filtered = df.collect()

        df_adj = adj_by_ctg_coords(df_filtered, "chrom").sort(by="chrom_st")
    except pl.exceptions.NoDataError:
        df_adj = pl.DataFrame(schema=BED9_COLS)

    if "item_rgb" not in df_adj.columns:
        df_adj = df_adj.with_columns(item_rgb=pl.lit("0,0,0"))
    if "name" not in df_adj.columns:
        df_adj = df_adj.with_columns(name=pl.lit("-"))

    return df_adj
