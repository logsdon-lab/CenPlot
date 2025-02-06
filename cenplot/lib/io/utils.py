import sys
import numpy as np
import polars as pl

from typing import Any
from matplotlib.colors import ListedColormap, rgb2hex

from ..track.types import NO_DATA_TRACK_OPTS, Track, TrackOption


def map_value_colors(
    df: pl.DataFrame,
    map_col: str | None = None,
    map_values: dict[Any, Any] | None = None,
    use_item_rgb: bool = False,
) -> pl.DataFrame:
    if "item_rgb" in df.columns and use_item_rgb:
        # Convert colors from rgb str -> rgb tuple -> hex
        # If already hex, allow.
        df = df.with_columns(
            color=pl.when(pl.col("item_rgb").str.starts_with("#"))
            .then(pl.col("item_rgb"))
            .otherwise(
                pl.col("item_rgb")
                .str.split(",")
                .list.eval(pl.element().cast(pl.Int16) / 255)
                .map_elements(lambda x: rgb2hex(x), return_dtype=pl.String)
            )
        )
    elif map_col:
        if map_values:
            val_color_mapping = map_values
        else:
            unique_vals = df[map_col].unique(maintain_order=True)
            # Generate random number of colors.
            colors = np.random.rand(len(unique_vals), 3)
            cmap = ListedColormap(colors)
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


def adj_by_ctg_coords(df: pl.DataFrame, colname: str) -> pl.DataFrame:
    final_adj_coords = {
        f"{colname}_st": pl.col(f"{colname}_st") - pl.col("ctg_st"),
        f"{colname}_end": pl.col(f"{colname}_end") - pl.col("ctg_st"),
    }
    return df.with_columns(
        chrom_name=pl.col(colname).str.extract(r"(chr[\dXY]+)").fill_null(""),
        # Use simplified coordinates if possible, otherwise, take everything.
        ctg_st=pl.col(colname).str.extract(r":(\d+)-").cast(pl.Int64).fill_null(0),
        ctg_end=pl.col(colname).str.extract(r"-(\d+)$").cast(pl.Int64).fill_null(0),
    ).with_columns(**final_adj_coords)


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
        # Skip tracks which carry no data.
        elif trk.opt in NO_DATA_TRACK_OPTS:
            continue
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
