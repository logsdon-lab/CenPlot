from dataclasses import dataclass
from typing import Literal


@dataclass
class DefaultPlotSettings:
    """
    Default plot options settings.
    """

    title: bool = True
    """
    Add `name` as title of y-axis label.
    """

    chrom_as_title: bool = False
    """
    Use chrom name as title.
    """

    legend: bool = True
    """
    Show the legend.
    """
    hide_x: bool = True
    """
    Hide the x-axis ticks, ticklabels, and spines.
    """


@dataclass
class SelfIdentPlotSettings(DefaultPlotSettings):
    """
    Self-identity heatmap triangle plot options.
    """

    invert: bool = True
    """
    Invert the self identity triangle.
    """
    legend_bins: int = 300
    """
    Number of bins for `perc_identity_by_events` in the legend.
    """
    legend_xmin: float = 70.0
    """
    Legend x-min coordinate. Used to constrain x-axis limits.
    """
    legend_asp_ratio: float | None = 1.0
    """
    Aspect ratio of legend. If `None`, takes up entire axis.
    """


@dataclass
class LabelPlotSettings(DefaultPlotSettings):
    """
    Label plot options.
    """

    DEF_COLOR = "black"
    """
    Default color for label.
    """

    color: str | None = None
    """
    Label color. Used if no color is provided in `item_rgb` column.
    """

    alpha: float = 1.0
    """
    Label alpha.
    """


@dataclass
class BarPlotSettings(DefaultPlotSettings):
    """
    Bar plot options.
    """

    DEF_COLOR = "black"
    """
    Default color for bar plot.
    """

    color: str | None = None
    """
    Color of bars. If `None`, uses `item_rgb` column colors.
    """

    alpha: float = 1.0
    """
    Alpha of bars.
    """


@dataclass
class HOROrtPlotSettings(DefaultPlotSettings):
    """
    Higher order repeat orientation arrow plot options.
    """

    DEF_COLOR = "black"
    """
    Default color for arrows.
    """
    scale: float = 50
    """
    Scale arrow attributes by this factor as well as length.
    """
    merge: int = 100_000
    """
    Merge orientations by this number of bases.
    """
    fwd_color: str | None = None
    """
    Color of `+` arrows.
    """
    rev_color: str | None = None
    """
    Color of `-` arrows.
    """


@dataclass
class HORPlotSettings(DefaultPlotSettings):
    """
    Higher order repeat plot options.
    """

    mer_order: Literal["large", "small"] = "large"
    """
    Plot HORs with `{mer_order}` monomers on top.
    """


PlotSettings = (
    HORPlotSettings
    | HOROrtPlotSettings
    | SelfIdentPlotSettings
    | BarPlotSettings
    | LabelPlotSettings
)
"""
Type annotation for all possible settings for the various plot types.
"""
