import os
import yaml
import tomllib
import logging
import polars as pl

from typing import Any, Generator, BinaryIO
from censtats.length import hor_array_length

from .utils import get_min_max_track, map_value_colors
from .bed9 import read_bed9
from .bed_identity import read_bed_identity
from .bed_label import read_bed_label
from .bed_hor import read_bed_hor, read_bed_hor_from_settings
from ..track.settings import (
    HORTrackSettings,
    HOROrtTrackSettings,
    LegendTrackSettings,
    PositionTrackSettings,
    SelfIdentTrackSettings,
    LabelTrackSettings,
    BarTrackSettings,
    TrackSettings,
    SpacerTrackSettings,
)
from ..track.types import Track, TrackType, TrackPosition, TrackList
from ..draw.settings import PlotSettings


def split_hor_track(
    df_track: pl.DataFrame,
    track_pos: TrackPosition,
    track_opt: TrackType,
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
        logging.error(
            f"Overlap not supported for {track_opt}. Using relative position.",
        )

    plot_options = HORTrackSettings(**options)
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
            TrackType.HORSplit,
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
        logging.error(f"Invalid plot position ({pos}) for {path}. Skipping.")
        return None
    try:
        track_opt = TrackType(opt)  # type: ignore[arg-type]
    except ValueError:
        logging.error(f"Invalid plot option ({opt}) for {path}. Skipping.")
        return None

    track_options: TrackSettings
    if track_opt == TrackType.Position:
        track_options = PositionTrackSettings(**options)
        track_options.hide_x = False
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None
    elif track_opt == TrackType.Legend:
        track_options = LegendTrackSettings(**options)
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None
    elif track_opt == TrackType.Spacer:
        track_options = SpacerTrackSettings(**options)
        yield Track(title, track_pos, track_opt, prop, None, track_options)
        return None

    if not path:
        raise ValueError("Path to data required.")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data does not exist for track ({track})")

    if track_opt == TrackType.HORSplit:
        df_track = read_bed_hor_from_settings(path, options, chrom)
        if df_track.is_empty():
            logging.error(
                f"Empty file or chrom not found for {track_opt} and {path}. Skipping"
            )
            return None
        if options.get("mode", HORTrackSettings.mode) == "hor":
            split_colname = "name"
        else:
            split_colname = "mer"
        split_prop = options.get("split_prop", HORTrackSettings.split_prop)
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

    elif track_opt == TrackType.HOR:
        df_track = read_bed_hor_from_settings(path, options, chrom)
        track_options = HORTrackSettings(**options)
        # Update legend title.
        if track_options.legend_title:
            track_options.legend_title = track_options.legend_title.format(chrom=chrom)

        yield Track(title, track_pos, track_opt, prop, df_track, track_options)
        return None

    if track_opt == TrackType.HOROrt:
        live_only = options.get("live_only", HOROrtTrackSettings.live_only)
        mer_filter = options.get("mer_filter", HOROrtTrackSettings.mer_filter)
        _, df_track = hor_array_length(
            read_bed_hor(
                path,
                chrom=chrom,
                live_only=live_only,
                mer_filter=mer_filter,
            ),
            output_strand=True,
        )
        track_options = HOROrtTrackSettings(**options)
    elif track_opt == TrackType.SelfIdent:
        df_track = read_bed_identity(path, chrom=chrom)
        track_options = SelfIdentTrackSettings(**options)
    elif track_opt == TrackType.Bar:
        df_track = read_bed9(path, chrom=chrom)
        track_options = BarTrackSettings(**options)
    else:
        use_item_rgb = options.get("use_item_rgb", LabelTrackSettings.use_item_rgb)
        df_track = read_bed_label(path, chrom=chrom)
        df_track = map_value_colors(
            df_track,
            map_col="name",
            use_item_rgb=use_item_rgb,
        )
        track_options = LabelTrackSettings(**options)

    df_track = map_value_colors(df_track)
    # Update legend title.
    if track_options.legend_title:
        track_options.legend_title = track_options.legend_title.format(chrom=chrom)

    yield Track(title, track_pos, track_opt, prop, df_track, track_options)


def read_one_cen_tracks(
    input_track: BinaryIO, *, chrom: str | None = None
) -> tuple[TrackList, PlotSettings]:
    """
    Read a `TOML` or `YAML` file of tracks to plot optionally filtering for a chrom name.

    Expected to have two items:
    * `[settings]`
        * See `cenplot.PlotSettings`
    * `[[tracks]]`
        * See one of the `cenplot.TrackSettings` for more details.

    Example:
    ```toml
    [settings]
    format = "png"
    transparent = true
    dim = [16.0, 8.0]
    dpi = 600
    ```

    ```yaml
    settings:
        format: "png"
        transparent: true
        dim: [16.0, 8.0]
        dpi: 600
    ```

    # Args:
    * input_track:
        * Input track `TOML` or `YAML` file.
    * chrom:
        * Chromosome name in 1st column (`chrom`) to filter for.
        * ex. `chr4`

    # Returns:
    * List of tracks w/contained chroms and plot settings.
    """
    all_tracks = []
    chroms = set()
    # Reset file position.
    input_track.seek(0)
    # Try TOML
    try:
        dict_settings = tomllib.load(input_track)
    except Exception:
        input_track.seek(0)
        # Then YAML
        try:
            dict_settings = yaml.safe_load(input_track)
        except Exception:
            raise TypeError("Invalid file type for settings.")

    settings: dict[str, Any] = dict_settings.get("settings", {})
    title = settings.get("title", PlotSettings.title)
    format = settings.get("format", PlotSettings.format)
    transparent = settings.get("transparent", PlotSettings.transparent)
    dim = tuple(settings.get("dim", PlotSettings.dim))
    dpi = settings.get("dpi", PlotSettings.dpi)
    legend_pos = settings.get("legend_pos", PlotSettings.legend_pos)
    legend_prop = settings.get("legend_prop", PlotSettings.legend_prop)
    axis_h_pad = settings.get("axis_h_pad", PlotSettings.axis_h_pad)
    layout = settings.get("layout", PlotSettings.layout)

    tracks = dict_settings.get("tracks", [])

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
    plot_settings = PlotSettings(
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
