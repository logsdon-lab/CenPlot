import polars as pl

from enum import StrEnum, auto
from typing import NamedTuple

from cenplot.track.settings import PlotSettings


class TrackPosition(StrEnum):
    Overlap = auto()
    Relative = auto()


class TrackOption(StrEnum):
    HOR = auto()
    HORSplit = auto()
    HOROrt = auto()
    Label = auto()
    Bar = auto()
    SelfIdent = auto()


class Track(NamedTuple):
    name: str
    pos: TrackPosition
    opt: TrackOption
    prop: float
    data: pl.DataFrame
    options: PlotSettings


class TrackList(NamedTuple):
    tracks: list[Track]
    chroms: set[str]
    min_pos: int
    max_pos: int


class LegendPosition(StrEnum):
    Left = auto()
    Right = auto()
