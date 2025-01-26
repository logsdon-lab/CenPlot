import polars as pl


def adj_by_ctg_coords(df: pl.DataFrame, colname: str) -> pl.DataFrame:
    final_adj_coords = {
        f"{colname}_st": pl.col(f"{colname}_st") - pl.col("ctg_st"),
        f"{colname}_end": pl.col(f"{colname}_end") - pl.col("ctg_st"),
    }
    return df.with_columns(
        chrom_name=pl.col(colname).str.extract(r"(chr[\dXY]+)").fill_null(""),
        # Use simplified coordinates if possible, otherwise, take everything.
        ctg_name=pl.col(colname).str.extract(r"^(.*?):|^(.*?)$"),
        ctg_st=pl.col(colname).str.extract(r":(\d+)-").cast(pl.Int64).fill_null(0),
        ctg_end=pl.col(colname).str.extract(r"-(\d+)$").cast(pl.Int64).fill_null(0),
    ).with_columns(
        # Use ctg name without coordinates. Avoids 0 and 1 offset issues.
        chrom=pl.col("ctg_name"),
        **final_adj_coords,
    )
