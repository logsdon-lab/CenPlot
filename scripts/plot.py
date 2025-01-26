import os
import argparse
import multiprocessing

import polars as pl

from concurrent.futures import ProcessPoolExecutor

from cenplot import plot_one_cen, merge_plots, read_all_tracks, Track


def get_inputs(
    args: argparse.Namespace,
) -> list[tuple[list[Track], str, str, str, int, int, int, int]]:
    tracks_summary = read_all_tracks(args.input_tracks)
    if args.chroms:
        all_chroms = args.chroms
    return [
        (
            [
                Track(
                    trk.name,
                    trk.pos,
                    trk.opt,
                    trk.prop,
                    trk.data.filter(pl.col("chrom") == chrom),
                    trk.options,
                )
                for trk in tracks_summary.tracks
            ],
            args.outdir,
            args.format,
            chrom,
            tracks_summary.min_pos,
            tracks_summary.max_pos,
            args.width,
            args.height,
        )
        for chrom in all_chroms
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-t",
        "--input_tracks",
        nargs="*",
        type=str,
        help=(
            "TOML file with headerless BED files to plot. "
            "Specify under tracks the following fields: {name, position, type, proportion, path, or options}. "
            "One of more TOML can be provided."
        ),
    )
    ap.add_argument(
        "-c",
        "--chroms",
        nargs="*",
        help="Names to plot in this order. Corresponds to 1st col in BED files.",
        default=None,
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
    ap.add_argument(
        "-f",
        "--format",
        help="Output format passed as ext to matplotlib. Without '.'",
        default="png",
    )
    ap.add_argument(
        "-w",
        "--width",
        type=float,
        help="Figure width in inches per centromere.",
        default=16.0,
    )
    ap.add_argument(
        "-ht",
        "--height",
        type=float,
        help="Figure height in inches per centromere.",
        default=8.0,
    )

    ap.add_argument("-p", "--processes", type=int, default=4, help="Processes to run.")
    args = ap.parse_args()

    tracks = get_inputs(args)

    os.makedirs(args.outdir, exist_ok=True)
    if args.processes == 1:
        plots = [plot_one_cen(*track) for track in tracks]
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            plots = pool.map(plot_one_cen, *zip(*tracks))  # type: ignore[assignment]

    if args.outfile:
        merge_plots(plots, args.outfile)


if __name__ == "__main__":
    raise SystemExit(main())
