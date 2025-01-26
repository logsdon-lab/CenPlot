import math
import polars as pl

from ..defaults import BED_SELF_IDENT_COLS, IDENT_COLOR_RANGE


def read_bed_identity(infile: str, *, chrom: str | None = None) -> pl.DataFrame:
    df = pl.read_csv(
        infile, separator="\t", has_header=False, new_columns=BED_SELF_IDENT_COLS
    ).with_columns(
        query_name=pl.col("query").str.extract(r"^(.*?):|^(.*?)$"),
        chrom_name=pl.col("query").str.extract(r"(chr[\dXY]+)").fill_null(""),
    )
    if chrom:
        df = df.filter(pl.col("chrom_name") == chrom)

    # Build expr to filter range of colors.
    # TODO: Allow custom range.
    color_expr = None
    for rng, color in IDENT_COLOR_RANGE.items():
        if not isinstance(color_expr, pl.Expr):
            color_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
        else:
            color_expr = color_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))

    if isinstance(color_expr, pl.Expr):
        color_expr = color_expr.otherwise(None)
    else:
        color_expr = pl.lit(None)

    tri_side = math.sqrt(2) / 2
    df_tri = (
        df.lazy()
        .with_columns(color=color_expr)
        # Get window size.
        .with_columns(
            window=(pl.col("query_end") - pl.col("query_st")).max().over("query_name")
        )
        .with_columns(
            first_pos=pl.col("query_st") // pl.col("window"),
            second_pos=pl.col("ref_st") // pl.col("window"),
        )
        # x y coords of diamond
        .with_columns(
            x=pl.col("first_pos") + pl.col("second_pos"),
            y=-pl.col("first_pos") + pl.col("second_pos"),
        )
        .with_columns(
            scale=(pl.col("query_st").max() / pl.col("x").max()).over("query_name"),
            group=pl.int_range(pl.len()).over("query_name"),
        )
        .with_columns(
            window=pl.col("window") / pl.col("scale"),
        )
        # Rather than generate new dfs. Add new x,y as arrays per row.
        .with_columns(
            new_x=[tri_side, 0.0, -tri_side, 0.0],
            new_y=[0.0, tri_side, 0.0, -tri_side],
        )
        # Convert to array to save memory
        .with_columns(
            pl.col("new_x").list.to_array(4), pl.col("new_y").list.to_array(4)
        )
        # Rescale x and y.
        .with_columns(
            ((pl.col("new_x") * pl.col("window")) + pl.col("x")) * pl.col("scale"),
            ((pl.col("new_y") * pl.col("window")) + pl.col("y")) * pl.col("window"),
        )
        .select(
            "query_name",
            "new_x",
            "new_y",
            "color",
            "group",
            "percent_identity_by_events",
        )
        # arr to new rows
        .explode("new_x", "new_y")
        # Rename to filter later on.
        .rename({"query_name": "chrom", "new_x": "x", "new_y": "y"})
        .collect()
    )
    return df_tri
