import os
import sys
import logging
import argparse
import matplotlib.colors
import polars as pl
from typing import Any, BinaryIO, TYPE_CHECKING

from cenplot import (
    read_tracks,
)
from cenplot.lib.draw.core import plot_comparison_tracks
from cenplot.lib.io.bed_identity import (
    expr_ident_to_color_n_range,
    read_ident_colorscale,
    read_bedpe,
)
from cenplot.lib.track.types import TrackType

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s.%(msecs)03d \033[32m%(levelname)s\033[0m [cenplot::%(name)s] %(message)s",
    datefmt="%Y-%m-%dT%H:%M:%S",
    handlers=[logging.StreamHandler(sys.stderr)],
)


def add_draw_cmp_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "draw_cmp",
        description="Draw comparison centromere tracks.",
    )
    ap.add_argument(
        "-i",
        "--ident",
        type=argparse.FileType("rb"),
        required=True,
        help="Identity matrix as BEDPE comparing reference and query.",
    )
    ap.add_argument(
        "-r",
        "--input_tracks_ref",
        required=True,
        type=argparse.FileType("rb"),
        help=(
            "Reference TOML or YAML file with headerless BED files to plot. "
            "Specify under tracks the following fields: {name, position, type, proportion, path, or options}."
        ),
    )
    ap.add_argument(
        "--pos_ref",
        default="bottom",
        type=str,
        choices=["top", "bottom"],
        help="Position of reference tracks.",
    )
    ap.add_argument(
        "-q",
        "--input_tracks_qry",
        required=True,
        type=argparse.FileType("rb"),
        help=(
            "Query TOML or YAML file with headerless BED files to plot. "
            "Specify under tracks the following fields: {name, position, type, proportion, path, or options}."
        ),
    )
    ap.add_argument(
        "--pos_qry",
        default="left",
        type=str,
        choices=["left", "right"],
        help="Position of query tracks.",
    )
    ap.add_argument(
        "--ident_colorscale",
        type=str,
        default=None,
        help="Identity matrix colorscale as 3-column TSV file (ident_st, ident_end, color).",
    )
    ap.add_argument(
        "--ident_prop",
        type=float,
        default=5.0,
        help="Identity matrix proporition for both height and width of plot.",
    )
    ap.add_argument(
        "-d",
        "--outdir",
        help="Output dir to plot multiple separate figures.",
        type=str,
        default=".",
    )

    return None


def draw_cmp(
    input_ident: BinaryIO,
    input_tracks_ref: BinaryIO,
    input_tracks_qry: BinaryIO,
    pos_tracks_ref: str,
    pos_tracks_qry: str,
    ident_colorscale: str,
    ident_prop: float,
    outdir: str,
):
    os.makedirs(outdir, exist_ok=True)

    tracks_ref, settings_ref = read_tracks(input_tracks_ref)
    tracks_qry, settings_qry = read_tracks(input_tracks_qry)

    df_bedpe = read_bedpe(input_ident, to_abs_coords=True)
    # Convert hex to rgb to plot.
    colorscale = {
        rng: matplotlib.colors.to_rgb(hx)
        for rng, hx in read_ident_colorscale(ident_colorscale).items()
    }
    color_expr, rng_expr = expr_ident_to_color_n_range(colorscale)
    window = (df_bedpe["query_end"] - df_bedpe["query_st"]).median()
    df_bedpe = df_bedpe.select(
        x=(pl.col("query_st") // window).cast(pl.UInt32),
        y=(pl.col("ref_st") // window).cast(pl.UInt32),
        color=color_expr,
        name=rng_expr,
    )
    # Convert coordinates to windows.
    for tracks in (tracks_qry.tracks, tracks_ref.tracks):
        for trk in tracks:
            if trk.opt == TrackType.SelfIdent:
                trk.data = trk.data.with_columns(pl.col("x") / window)
            else:
                trk.data = trk.data.with_columns(
                    pl.col("chrom_st") / window, pl.col("chrom_end") / window
                )

    plot_comparison_tracks(
        tracks_ref=tracks_ref.tracks,
        tracks_qry=tracks_qry.tracks,
        settings_ref=settings_ref,
        settings_qry=settings_qry,
        pos_ref=pos_tracks_ref,
        pos_qry=pos_tracks_qry,
        df_ident=df_bedpe,
        prop_ident=ident_prop,
        window=window,
    )

    logging.info("Done!")
