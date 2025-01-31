import os
import sys
import tomllib

import polars as pl

from typing import Any, Generator
from censtats.length import hor_array_length

from .bed9 import read_bed9
from .bed_identity import read_bed_identity
from .bed_label import read_bed_label
from .bed_hor import read_bed_hor
from ..track.settings import (
    HORPlotSettings,
    HOROrtPlotSettings,
    SelfIdentPlotSettings,
    LabelPlotSettings,
    BarPlotSettings,
    PlotSettings,
)
from ..track.types import Track, TrackOption, TrackPosition, TrackList
from ..draw.settings import SinglePlotSettings


def read_one_track_info(
    track: dict[str, Any], *, chrom: str | None = None
) -> Generator[Track, None, None]:
    prop = track.get("proportion", 0.0)
    title = track.get("title")
    pos = track.get("position")
    opt = track.get("type")
    path: str | None = track.get("path")
    options: dict[str, Any] = track.get("options", {})

    try:
        track_pos = TrackPosition(pos)  # type: ignore[arg-type]
    except ValueError:
        print(
            f"Invalid plot position ({pos}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None
    try:
        track_opt = TrackOption(opt)  # type: ignore[arg-type]
    except ValueError:
        print(
            f"Invalid plot option ({opt}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None

    if not path:
        raise ValueError("Path to data required.")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data does not exist for track ({track})")

    if track_opt == TrackOption.HORSplit:
        mer_order = options.get("mer_order", HORPlotSettings.mer_order)
        df_track = read_bed_hor(path, chrom=chrom, mer_order=mer_order)
        if df_track.is_empty():
            print(
                f"Empty file or chrom not found for {track_opt} and {path}. Skipping",
                file=sys.stderr,
            )
            return None

        uniq_mers = df_track["mer"].unique()
        track_prop = prop / len(uniq_mers)
        if track_pos == TrackPosition.Overlap:
            print(
                f"Overlap not supported for {track_opt}. Using relative position.",
                file=sys.stderr,
            )

        for mer, df_mer_track in df_track.group_by(["mer"], maintain_order=True):
            mer = mer[0]
            # Add mer to name if formatted.
            try:
                mer_title = str(title).format(mer=mer) if title else ""
            except KeyError:
                mer_title = str(title) if title else ""
            # Disallow overlap.
            # Split proportion over uniq monomers.
            yield Track(
                mer_title,
                TrackPosition.Relative,
                TrackOption.HORSplit,
                track_prop,
                df_mer_track,
                HORPlotSettings(**options),
            )

        return None

    track_options: PlotSettings
    if track_opt == TrackOption.HOR:
        mer_order = options.get("mer_order", HORPlotSettings.mer_order)
        df_track = read_bed_hor(path, chrom=chrom, mer_order=mer_order)
        track_options = HORPlotSettings(**options)
    elif track_opt == TrackOption.HOROrt:
        _, df_track = hor_array_length(
            read_bed_hor(path, chrom=chrom), output_strand=True
        )
        # Merged HOR arrays must be 90% HORs.
        df_track = df_track.filter(pl.col("prop") > 0.9)
        track_options = HOROrtPlotSettings(**options)
    elif track_opt == TrackOption.SelfIdent:
        df_track = read_bed_identity(path, chrom=chrom)
        track_options = SelfIdentPlotSettings(**options)
    elif track_opt == TrackOption.Bar:
        df_track = read_bed9(path, chrom=chrom)
        track_options = BarPlotSettings(**options)
    else:
        df_track = read_bed_label(path, chrom=chrom)
        track_options = LabelPlotSettings(**options)

    yield Track(title, track_pos, track_opt, prop, df_track, track_options)


def get_min_max_track(
    tracks: list[Track], typ: str, default_col: str = "chrom_st"
) -> tuple[Track, int]:
    track = None
    if typ == "min":
        pos = sys.maxsize
    else:
        pos = 0

    for trk in tracks:
        if trk.opt == TrackOption.SelfIdent:
            col = "x"
        else:
            col = default_col
        if typ == "min":
            trk_min = trk.data.filter(pl.col(col) > 0)[col].min()
            if trk_min < pos:
                track = trk
                pos = trk_min
        else:
            trk_max = trk.data[col].max()
            if trk_max > pos:
                track = trk
                pos = trk_max
    if not track:
        raise ValueError("No tracks.")
    return track, pos


def read_one_cen_tracks(
    input_track: str, *, chrom: str | None = None
) -> tuple[TrackList, SinglePlotSettings]:
    """
    Read a `cenplot` track file optionally filtering for a chrom name.

    Args:
        input_track:
            Input track file.
        chrom:
            Chromosome name in 1st column (`chrom`) to filter for. ex. `chr4`

    Returns:
        List of tracks w/contained chroms and plot settings.
    """
    all_tracks = []
    chroms = set()
    with open(input_track, "rb") as fh:
        toml = tomllib.load(fh)
        settings: dict[str, Any] = toml.get("settings", {})
        format = settings.get("format", SinglePlotSettings.format)
        transparent = settings.get("transparent", SinglePlotSettings.transparent)
        dim = tuple(settings.get("dim", SinglePlotSettings.dim))
        dpi = settings.get("dpi", SinglePlotSettings.dpi)
        legend_pos = settings.get("legend_pos", SinglePlotSettings.legend_pos)
        legend_prop = settings.get("legend_prop", SinglePlotSettings.legend_prop)
        axis_h_pad = settings.get("axis_h_pad", SinglePlotSettings.axis_h_pad)
        layout = settings.get("layout", SinglePlotSettings.layout)

        tracks = toml.get("tracks", [])

        for track_info in tracks:
            for track in read_one_track_info(track_info, chrom=chrom):
                all_tracks.append(track)
                chroms.update(track.data["chrom"])

    _, min_st_pos = get_min_max_track(all_tracks, typ="min")
    _, max_end_pos = get_min_max_track(all_tracks, typ="max", default_col="chrom_end")
    tracklist = TrackList(all_tracks, chroms)
    plot_settings = SinglePlotSettings(
        format,
        transparent,
        dim,
        dpi,
        layout,
        legend_pos,
        legend_prop,
        axis_h_pad,
        shared_xlim=(
            tuple(settings.get("shared_xlim"))  # type: ignore[arg-type]
            if settings.get("shared_xlim")
            else (min_st_pos, max_end_pos)
        ),
    )
    return tracklist, plot_settings
