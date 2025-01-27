import polars as pl

from enum import StrEnum, auto
from typing import Any, NamedTuple

__all__ = ["Track", "TrackOption", "TrackPosition", "TrackList", "LegendPosition"]


class TrackPosition(StrEnum):
    Overlap = auto()
    Relative = auto()


class TrackOption(StrEnum):
    HOR = auto()
    HORSplit = auto()
    Label = auto()
    Value = auto()
    SelfIdent = auto()


class Track(NamedTuple):
    name: str
    pos: TrackPosition
    opt: TrackOption
    prop: float
    data: pl.DataFrame
    options: dict[str, Any]


class TrackList(NamedTuple):
    tracks: list[Track]
    chroms: set[str]
    min_pos: int
    max_pos: int


class LegendPosition(StrEnum):
    Left = auto()
    Right = auto()
