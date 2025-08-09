import os
import sys
import logging
import argparse
import multiprocessing

import polars as pl

from typing import Any, BinaryIO, TYPE_CHECKING
from concurrent.futures import ProcessPoolExecutor

from cenplot import (
    plot_one_cen,
    merge_plots,
    read_one_cen_tracks,
    Track,
    PlotSettings,
)

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


def get_draw_args(
    input_tracks: BinaryIO, chroms: list[str], share_xlim: bool, outdir: str
) -> list[tuple[list[Track], str, str, PlotSettings]]:
    inputs = []
    tracks_settings = [
        (chrom, *read_one_cen_tracks(input_tracks, chrom=chrom)) for chrom in chroms
    ]
    xmin_all, xmax_all = sys.maxsize, 0
    if share_xlim:
        for *_, settings in tracks_settings:
            if settings.xlim:
                xmin, xmax = settings.xlim
                xmin_all = min(xmin_all, xmin)
                xmax_all = max(xmax_all, xmax)

    for chrom, tracks_summary, plot_settings in tracks_settings:
        if share_xlim:
            plot_settings.xlim = (xmin_all, xmax_all)
        chrom_no_coords = chrom.rsplit(":", 1)[0]

        tracks = []
        for trk in tracks_summary.tracks:
            df = trk.data

            if isinstance(df, pl.DataFrame):
                try:
                    has_no_coords = chrom_no_coords in df["chrom"]
                except Exception:
                    has_no_coords = False

                if has_no_coords:
                    df = df.filter(pl.col("chrom") == chrom_no_coords)
                else:
                    df = df.filter(pl.col("chrom") == chrom)

            tracks.append(
                Track(
                    trk.title,
                    trk.pos,
                    trk.opt,
                    trk.prop,
                    df if isinstance(df, pl.DataFrame) else None,
                    trk.options,
                )
            )

        inputs.append(
            (
                tracks,
                outdir,
                chrom,
                plot_settings,
            )
        )
    return inputs


def add_draw_cli(parser: SubArgumentParser) -> None:
    ap = parser.add_parser(
        "draw",
        description="Draw centromere tracks.",
    )
    ap.add_argument(
        "-t",
        "--input_tracks",
        required=True,
        type=argparse.FileType("rb"),
        help=(
            "TOML or YAML file with headerless BED files to plot. "
            "Specify under tracks the following fields: {name, position, type, proportion, path, or options}."
        ),
    )
    ap.add_argument(
        "-c",
        "--chroms",
        nargs="+",
        help="Regions to plot in this order. Corresponds to 1st col in BED files. If contains ':' and coords, subsets to those coordinates.",
        required=True,
    )
    ap.add_argument(
        "-d",
        "--outdir",
        help="Output dir to plot multiple separate figures.",
        type=str,
        default=".",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        help="Output file merging all figures. Either pdf of png.",
        type=str,
        default=None,
    )
    ap.add_argument("--share_xlim", help="Share x-axis limits.", action="store_true")
    ap.add_argument("-p", "--processes", type=int, default=4, help="Processes to run.")

    return None


def draw(
    input_tracks: BinaryIO,
    chroms: list[str],
    outdir: str,
    outfile: str,
    share_xlim: bool,
    processes: int,
):
    draw_args = get_draw_args(
        input_tracks=input_tracks, chroms=chroms, share_xlim=share_xlim, outdir=outdir
    )

    os.makedirs(outdir, exist_ok=True)
    if processes == 1:
        plots = [plot_one_cen(*draw_arg) for draw_arg in draw_args]
    else:
        with ProcessPoolExecutor(
            max_workers=processes, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            futures = [
                (draw_arg[2], pool.submit(plot_one_cen, *draw_arg))
                for draw_arg in draw_args
            ]  # type: ignore[assignment]
            plots = []
            for chrom, future in futures:
                if future.exception():
                    logging.error(f"Failed to plot {chrom} ({future.exception()})")
                    continue
                plots.append(future.result())

    if outfile:
        logging.info(f"Merging {len(plots)} plots into {outfile}.")
        merge_plots(plots, outfile)

    logging.info("Done!")
