import os
import sys
import tomllib

import numpy as np
import polars as pl

from typing import Any, Generator
from censtats.length import hor_array_length
from matplotlib.colors import LinearSegmentedColormap, rgb2hex

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
from ..defaults import MONOMER_COLORS


def map_value_colors(
    df: pl.DataFrame,
    map_col: str | None = None,
    map_values: dict[Any, Any] | None = None,
    use_item_rgb: bool = False,
) -> pl.DataFrame:
    if "item_rgb" in df.columns and use_item_rgb:
        # Convert colors from rgb str -> rgb tuple -> hex
        df = df.with_columns(
            color=pl.col("item_rgb")
            .str.split(",")
            .list.eval(pl.element().cast(pl.Int16) / 255)
            .map_elements(lambda x: rgb2hex(x), return_dtype=pl.String)
        )
    elif map_col:
        if map_values:
            val_color_mapping = map_values
        else:
            unique_vals = df[map_col].unique(maintain_order=True)
            cmap = LinearSegmentedColormap.from_list(
                "rg", ["r", "w", "g"], N=len(unique_vals)
            )
            val_color_mapping = {
                val: rgb2hex(color)
                for val, color in zip(
                    unique_vals, cmap(np.linspace(0, 1, len(unique_vals)))
                )
            }
        df = df.with_columns(
            color=pl.col(map_col)
            .cast(pl.String)
            # If not in mapping, set to gray.
            .replace(val_color_mapping, default="#808080")
        )

    return df


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

    track_options: PlotSettings
    if track_opt == TrackOption.HORSplit:
        mer_order = options.get("mer_order", HORPlotSettings.mer_order)
        live_only = options.get("live_only", HORPlotSettings.live_only)
        mer_filter = options.get("mer_filter", HORPlotSettings.mer_filter)
        split_prop = options.get("split_prop", HORPlotSettings.split_prop)
        df_track = read_bed_hor(
            path,
            chrom=chrom,
            mer_order=mer_order,
            live_only=live_only,
            mer_filter=mer_filter,
        )
        if df_track.is_empty():
            print(
                f"Empty file or chrom not found for {track_opt} and {path}. Skipping",
                file=sys.stderr,
            )
            return None

        # Use item_rgb column otherwise, map name or mer to a color.
        use_item_rgb = options.get("use_item_rgb", HORPlotSettings.use_item_rgb)
        if options.get("mode", HORPlotSettings.mode) == "hor":
            split_colname = "name"
            df_track = map_value_colors(
                df_track,
                map_col=split_colname,
                use_item_rgb=use_item_rgb,
            )
        else:
            split_colname = "mer"
            df_track = map_value_colors(
                df_track,
                map_col=split_colname,
                map_values=MONOMER_COLORS,
                use_item_rgb=use_item_rgb,
            )

        srs_split_names = df_track[split_colname].unique()
        # Split proportion across tracks.
        if split_prop:
            track_prop = prop / len(srs_split_names)
        else:
            track_prop = prop

        if track_pos == TrackPosition.Overlap:
            print(
                f"Overlap not supported for {track_opt}. Using relative position.",
                file=sys.stderr,
            )

        for split, df_split_track in df_track.group_by(
            [split_colname], maintain_order=True
        ):
            split = split[0]
            # Add mer to name if formatted.
            try:
                mer_title = str(title).format(mer=split) if title else ""
            except KeyError:
                mer_title = str(title) if title else ""
            # Disallow overlap.
            # Split proportion over uniq monomers.
            yield Track(
                mer_title,
                TrackPosition.Relative,
                TrackOption.HORSplit,
                track_prop,
                df_split_track,
                HORPlotSettings(**options),
            )

        return None

    elif track_opt == TrackOption.HOR:
        mer_order = options.get("mer_order", HORPlotSettings.mer_order)
        live_only = options.get("live_only", HORPlotSettings.live_only)
        mer_filter = options.get("mer_filter", HORPlotSettings.mer_filter)

        # Use item_rgb column otherwise, map name or mer to a color.
        use_item_rgb = options.get("use_item_rgb", HORPlotSettings.use_item_rgb)
        df_track = read_bed_hor(
            path,
            chrom=chrom,
            mer_order=mer_order,
            live_only=live_only,
            mer_filter=mer_filter,
        )
        df_track = map_value_colors(
            df_track,
            map_col="mer",
            map_values=MONOMER_COLORS,
            use_item_rgb=use_item_rgb,
        )
        track_options = HORPlotSettings(**options)

        yield Track(title, track_pos, track_opt, prop, df_track, track_options)
        return None

    if track_opt == TrackOption.HOROrt:
        live_only = options.get("live_only", HOROrtPlotSettings.live_only)
        mer_filter = options.get("mer_filter", HOROrtPlotSettings.mer_filter)
        _, df_track = hor_array_length(
            read_bed_hor(
                path,
                chrom=chrom,
                live_only=live_only,
                mer_filter=mer_filter,
                visualization=False,
            ),
            output_strand=True,
        )
        track_options = HOROrtPlotSettings(**options)
    elif track_opt == TrackOption.SelfIdent:
        df_track = read_bed_identity(path, chrom=chrom)
        track_options = SelfIdentPlotSettings(**options)
    elif track_opt == TrackOption.Bar:
        df_track = read_bed9(path, chrom=chrom)
        track_options = BarPlotSettings(**options)
    else:
        use_item_rgb = options.get("use_item_rgb", LabelPlotSettings.use_item_rgb)
        df_track = read_bed_label(path, chrom=chrom)
        df_track = map_value_colors(
            df_track,
            map_col="name",
            use_item_rgb=use_item_rgb,
        )
        track_options = LabelPlotSettings(**options)

    df_track = map_value_colors(df_track)

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
        title = settings.get("title", SinglePlotSettings.title)
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
        title,
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
