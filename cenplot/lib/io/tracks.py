import os
import sys
import tomllib
import polars as pl


from typing import Any, Generator
from censtats.length import hor_array_length

from .utils import get_min_max_track, map_value_colors
from .bed9 import read_bed9
from .bed_identity import read_bed_identity
from .bed_label import read_bed_label
from .bed_hor import read_bed_hor
from ..track.settings import (
    HORPlotSettings,
    HOROrtPlotSettings,
    LegendPlotSettings,
    PositionPlotSettings,
    SelfIdentPlotSettings,
    LabelPlotSettings,
    BarPlotSettings,
    PlotSettings,
    SpacerPlotSettings,
)
from ..track.types import Track, TrackOption, TrackPosition, TrackList
from ..draw.settings import SinglePlotSettings


def split_hor_track(
    df_track: pl.DataFrame,
    track_pos: TrackPosition,
    track_opt: TrackOption,
    title: Any | None,
    prop: float,
    split_colname: str,
    split_prop: bool,
    options: dict[str, Any],
    chrom: str | None = None,
) -> Generator[Track, None, None]:
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

    plot_options = HORPlotSettings(**options)
    for split, df_split_track in df_track.group_by(
        [split_colname], maintain_order=True
    ):
        split = split[0]
        # Add mer to name if formatted.
        try:
            mer_title = str(title).format(**{split_colname: split}) if title else ""
        except KeyError:
            mer_title = str(title) if title else ""

        # Update legend title.
        if plot_options.legend_title and chrom:
            plot_options.legend_title = plot_options.legend_title.format(
                **{split_colname: split, "chrom": chrom}
            )

        # Disallow overlap.
        # Split proportion over uniq monomers.
        yield Track(
            mer_title,
            TrackPosition.Relative,
            TrackOption.HORSplit,
            track_prop,
            df_split_track,
            plot_options,
        )


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

    track_options: PlotSettings
    if track_opt == TrackOption.Position:
        track_options = PositionPlotSettings(**options)
        track_options.hide_x = False
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None
    elif track_opt == TrackOption.Legend:
        track_options = LegendPlotSettings(**options)
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None
    elif track_opt == TrackOption.Spacer:
        track_options = SpacerPlotSettings(**options)
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None

    if not path:
        raise ValueError("Path to data required.")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data does not exist for track ({track})")

    if track_opt == TrackOption.HORSplit:
        live_only = options.get("live_only", HORPlotSettings.live_only)
        mer_filter = options.get("mer_filter", HORPlotSettings.mer_filter)
        hor_filter = options.get("hor_filter", HORPlotSettings.hor_filter)
        split_prop = options.get("split_prop", HORPlotSettings.split_prop)
        use_item_rgb = options.get("use_item_rgb", HORPlotSettings.use_item_rgb)
        sort_order = options.get("sort_order", HORPlotSettings.sort_order)

        # Use item_rgb column otherwise, map name or mer to a color.
        if options.get("mode", HORPlotSettings.mode) == "hor":
            split_colname = "name"
        else:
            split_colname = "mer"

        df_track = read_bed_hor(
            path,
            chrom=chrom,
            sort_col=split_colname,
            sort_order=sort_order,
            live_only=live_only,
            mer_filter=mer_filter,
            hor_filter=hor_filter,
            use_item_rgb=use_item_rgb,
        )
        if df_track.is_empty():
            print(
                f"Empty file or chrom not found for {track_opt} and {path}. Skipping",
                file=sys.stderr,
            )
            return None

        yield from split_hor_track(
            df_track,
            track_pos,
            track_opt,
            title,
            prop,
            split_colname,
            split_prop,
            options,
            chrom=chrom,
        )
        return None

    elif track_opt == TrackOption.HOR:
        sort_order = options.get("sort_order", HORPlotSettings.sort_order)
        live_only = options.get("live_only", HORPlotSettings.live_only)
        mer_filter = options.get("mer_filter", HORPlotSettings.mer_filter)
        hor_filter = options.get("hor_filter", HORPlotSettings.hor_filter)

        # Use item_rgb column otherwise, map name or mer to a color.
        use_item_rgb = options.get("use_item_rgb", HORPlotSettings.use_item_rgb)
        df_track = read_bed_hor(
            path,
            chrom=chrom,
            sort_col="mer",
            sort_order=sort_order,
            live_only=live_only,
            mer_filter=mer_filter,
            hor_filter=hor_filter,
            use_item_rgb=use_item_rgb,
        )
        track_options = HORPlotSettings(**options)
        # Update legend title.
        if track_options.legend_title:
            track_options.legend_title = track_options.legend_title.format(chrom=chrom)

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
    # Update legend title.
    if track_options.legend_title:
        track_options.legend_title = track_options.legend_title.format(chrom=chrom)

    yield Track(title, track_pos, track_opt, prop, df_track, track_options)


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
                # Tracks legend and position have no data.
                if not isinstance(track.data, pl.DataFrame):
                    continue
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
        xlim=(
            tuple(settings.get("xlim"))  # type: ignore[arg-type]
            if settings.get("xlim")
            else (min_st_pos, max_end_pos)
        ),
    )
    return tracklist, plot_settings
