import polars as pl

from enum import StrEnum, auto
from typing import NamedTuple

from ..track.settings import TrackSettings


class TrackPosition(StrEnum):
    Overlap = auto()
    Relative = auto()


class TrackType(StrEnum):
    """
    Track options.
    * Input track data is expected to be headerless.
    """

    HOR = auto()
    """
    An alpha-satellite higher order repeat (HOR) track with HORs by monomer number overlapping.

    Expected format:
    * [`BED9`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        * `name` as HOR variant
            * ex. `S4CYH1L.44-1`
    """
    HORSplit = auto()
    """
    A split alpha-satellite higher order repeat (HOR) track with each type of HOR as a single track.
    * `mer` or the number of monomers within the HOR.
    * `hor` or HOR variant.

    Expected format:
    * [`BED9`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        * `name` as HOR variant
            * ex. `S4CYH1L.44-1`
    """
    HOROrt = auto()
    """
    An alpha-satellite higher order repeat (HOR) orientation track.
    * This is calculate with default settings via the [`censtats`](https://github.com/logsdon-lab/CenStats) library.

    Expected format:
    * [`BED9`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        * `name` as HOR variant
            * ex. `S4CYH1L.44-1`
        * `strand` as `+` or `-`
    """
    Label = auto()
    """
    A label track. Elements in the `name` column are displayed as colored rectangles.

    Expected format:
    * [`BED4-9`](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        * `name` as any string value.
    """
    Bar = auto()
    """
    A bar plot track. Elements in the `name` column are displayed as bars.

    Expected format:
    * `BED9`
        * `name` as any numeric value.
    """
    SelfIdent = auto()
    """
    A self, sequence identity heatmap track displayed as a triangle.
    * Similar to plots from [`ModDotPlot`](https://github.com/marbl/ModDotPlot)

    Expected format:
    * `BED*`
        * Identity bedfile produced by `ModDotPlot` without a header.

    |query|query_st|query_end|reference|reference_st|reference_end|percent_identity_by_events|
    |-|-|-|-|-|-|-|
    |x|1|5000|x|1|5000|100.0|

    """
    Position = auto()
    """
    Position track.
    * Displays the x-axis position as well as a label.

    Expected format:
    * None
    """

    Legend = auto()
    """
    Legend track. Displays the legend of a specified track.
    * NOTE: This does not work with `TrackType.HORSplit`

    Expected format:
    * None
    """

    Spacer = auto()
    """
    Spacer track. Empty space.

    Expected format:
    * None
    """


class Track(NamedTuple):
    """
    A centromere track.
    """

    title: str | None
    """
    Title of track.
    * ex. "{chrom}"
    * ex. "HOR monomers"
    """
    pos: TrackPosition
    """
    Track position.
    """
    opt: TrackType
    """
    Track option.
    """
    prop: float
    """
    Proportion of track in final figure.
    """
    data: pl.DataFrame
    """
    Track data.
    """
    options: TrackSettings
    """
    Plot settings.
    """


class TrackList(NamedTuple):
    """
    Track list.
    """

    tracks: list[Track]
    """
    Tracks.
    """
    chroms: set[str]
    """
    Chromosomes found with `tracks`.
    """


class LegendPosition(StrEnum):
    Left = auto()
    Right = auto()


NO_DATA_TRACK_OPTS = {TrackType.Legend, TrackType.Position, TrackType.Spacer}
