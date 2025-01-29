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
    title: str | None
    """
    Title of track.
    ex. "{chrom}"
    ex. HOR monomers
    """
    pos: TrackPosition
    """
    Track position.
    """
    opt: TrackOption
    prop: float
    data: pl.DataFrame
    options: PlotSettings


class TrackList(NamedTuple):
    tracks: list[Track]
    chroms: set[str]


class LegendPosition(StrEnum):
    Left = auto()
    Right = auto()
