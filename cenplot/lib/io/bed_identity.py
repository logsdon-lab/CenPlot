import csv
import logging
import math
import polars as pl

from typing import TextIO

from ..track.settings import LocalSelfIdentTrackSettings
from ..defaults import BED9_COLS, BED_SELF_IDENT_COLS, IDENT_COLORSCALE, Colorscale
from censtats.self_ident.cli import convert_2D_to_1D_ident, Dim


def read_ident_colorscale(
    colorscale: Colorscale | str | None,
) -> dict[tuple[float, float], str]:
    if not colorscale:
        return IDENT_COLORSCALE

    if isinstance(colorscale, dict):
        return colorscale

    ident_colorscale = {}
    with open(colorscale, "rt") as fh:
        for line in fh:
            st, end, color, *_ = line.strip().split("\t")
            fst = float(st)
            fend = float(end)
            ident_colorscale[(fst, fend)] = color
    return ident_colorscale


def read_bed_identity(
    infile: str | TextIO,
    *,
    chrom: str | None = None,
    mode: str = "2D",
    colorscale: Colorscale | str | None = None,
    band_size: int = LocalSelfIdentTrackSettings.band_size,
    ignore_band_size=LocalSelfIdentTrackSettings.ignore_band_size,
) -> tuple[pl.DataFrame, Colorscale]:
    """
    Read a self, sequence identity BED file generate by `ModDotPlot`.

    Requires the following columns
    * `query,query_st,query_end,ref,ref_st,ref_end,percent_identity_by_events`

    # Args
    * `infile`
        * File or IO stream.
    * `chrom`
        * Chromosome name in `query` column to filter for.
    * `mode`
        * 1D or 2D self-identity.
    * `band_size`
        * Number of windows to calculate average sequence identity over. Only applicable if mode is 1D.
    * `ignore_band_size`
        * Number of windows ignored along self-identity diagonal. Only applicable if mode is 1D.

    # Returns
    * Coordinates of colored polygons in 2D space.
    """
    sniffer = csv.Sniffer()
    fname = infile if isinstance(infile, str) else infile.name
    with open(fname, "rt") as fh:
        header = next(fh)

    skip_rows = 1 if sniffer.has_header(header) else 0

    df = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=BED_SELF_IDENT_COLS,
        skip_rows=skip_rows,
    )
    if chrom:
        df = df.filter(pl.col("query") == chrom)

    # Check mode. Set by dev not user.
    mode = Dim(mode)

    # Build expr to filter range of colors.
    color_expr = None
    rng_expr = None
    ident_colorscale = read_ident_colorscale(colorscale)
    for rng, color in ident_colorscale.items():
        if not isinstance(color_expr, pl.Expr):
            color_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
            rng_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(f"{rng[0]}-{rng[1]}"))
        else:
            color_expr = color_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
            rng_expr = rng_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(f"{rng[0]}-{rng[1]}"))

    if isinstance(color_expr, pl.Expr):
        color_expr = color_expr.otherwise(None)
    else:
        color_expr = pl.lit(None)
    if isinstance(rng_expr, pl.Expr):
        rng_expr = rng_expr.otherwise(None)
    else:
        rng_expr = pl.lit(None)

    if mode == Dim.ONE:
        df_window = (
            (df["query_end"] - df["query_st"])
            .value_counts(sort=True)
            .rename({"query_end": "window"})
        )
        if df_window.shape[0] > 1:
            logging.warning(f"Multiple windows detected. Taking largest.\n{df_window}")
        window = df_window.row(0, named=True)["window"] + 1
        df_local_ident = pl.DataFrame(
            convert_2D_to_1D_ident(df.iter_rows(), window, band_size, ignore_band_size),
            schema=[
                "chrom_st",
                "chrom_end",
                "percent_identity_by_events",
            ],
            orient="row",
        )
        query = df["query"][0]
        df_res = (
            df_local_ident.lazy()
            .with_columns(
                chrom=pl.lit(query),
                color=color_expr,
                name=rng_expr,
                score=pl.col("percent_identity_by_events"),
                strand=pl.lit("."),
                thick_st=pl.col("chrom_st"),
                thick_end=pl.col("chrom_end"),
                item_rgb=pl.lit("0,0,0"),
            )
            .select(*BED9_COLS, "color")
            .collect()
        )
    else:
        tri_side = math.sqrt(2) / 2
        df_res = (
            df.lazy()
            .with_columns(color=color_expr)
            # Get window size.
            .with_columns(
                window=(pl.col("query_end") - pl.col("query_st")).max().over("query")
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
                scale=(pl.col("query_st").max() / pl.col("x").max()).over("query"),
                group=pl.int_range(pl.len()).over("query"),
            )
            .with_columns(
                window=pl.col("window") / pl.col("scale"),
            )
            # Rather than generate new dfs. Add new x,y as arrays per row.
            .with_columns(
                new_x=[tri_side, 0.0, -tri_side, 0.0],
                new_y=[0.0, tri_side, 0.0, -tri_side],
            )
            # Rescale x and y.
            .with_columns(
                ((pl.col("new_x") * pl.col("window")) + pl.col("x")) * pl.col("scale"),
                ((pl.col("new_y") * pl.col("window")) + pl.col("y")) * pl.col("window"),
            )
            .select(
                "query",
                "new_x",
                "new_y",
                "color",
                "group",
                "percent_identity_by_events",
            )
            # arr to new rows
            .explode("new_x", "new_y")
            # Rename to filter later on.
            .rename({"query": "chrom", "new_x": "x", "new_y": "y"})
            .collect()
        )
    return df_res, ident_colorscale
