import polars as pl

from .bed9 import read_bed9
from ..defaults import MONOMER_LEN


def read_bed_hor(
    infile: str, *, chrom: str | None = None, mer_order: str = "large"
) -> pl.DataFrame:
    return (
        read_bed9(infile, chrom=chrom)
        .with_columns(
            length=pl.col("chrom_end") - pl.col("chrom_st"),
        )
        .with_columns(mer=(pl.col("length") / MONOMER_LEN).round().cast(pl.Int8))
        .filter(
            pl.when(pl.col("chrom_name") == "chr10")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chr20")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chrY")
            .then(pl.col("mer") >= 30)
            .when(pl.col("chrom_name") == "chr17")
            .then(pl.col("mer") >= 4)
            .otherwise(True)
        )
        .sort("mer", descending=mer_order == "large")
        .cast({"mer": pl.String})
        # Add 1000 bp for visualization purposes.
        .with_columns(pl.col("length") + 2000)
    )
