import polars as pl

from .bed9 import read_bed9
from ..defaults import MONOMER_LEN


def read_bed_hor(
    infile: str,
    *,
    chrom: str | None = None,
    mer_order: str = "large",
    live_only: bool = True,
    mer_filter: int = 2,
    visualization: bool = False,
) -> pl.DataFrame:
    lf = (
        read_bed9(infile, chrom=chrom)
        .lazy()
        .with_columns(
            length=pl.col("chrom_end") - pl.col("chrom_st"),
        )
        .with_columns(
            mer=(pl.col("length") / MONOMER_LEN).round().cast(pl.Int8).clip(1, 100)
        )
    )

    if visualization:
        lf = lf.filter(
            pl.when(pl.col("chrom_name") == "chr10")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chr20")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chrY")
            .then(pl.col("mer") >= 30)
            .when(pl.col("chrom_name") == "chr17")
            .then(pl.col("mer") >= 4)
            .otherwise(True)
        ).with_columns(pl.col("length") + 2000)
    return (
        lf.filter(
            pl.when(live_only).then(pl.col("name").str.contains("L")).otherwise(True)
            & (pl.col("mer") >= mer_filter)
        )
        .sort("mer", descending=mer_order == "large")
        .cast({"mer": pl.String})
        .collect()
    )
