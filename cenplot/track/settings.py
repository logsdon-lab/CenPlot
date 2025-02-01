from dataclasses import dataclass
from typing import Literal


@dataclass
class DefaultPlotSettings:
    """
    Default plot options settings.
    """

    fontsize: float | str = "medium"
    """
    Font size for track text.
    """

    title_fontsize: float | str = "x-large"
    """
    Font size for track title.
    """

    legend: bool = True
    """
    Show the legend.
    """

    legend_fontsize: str = "medium"
    """
    Legend font size.
    """

    legend_title: str | None = None
    """
    Set legend title.
    * ex. "HOR monomers for {chrom}"
    """

    legend_title_fontsize: str = "x-large"
    """
    Legend title font size.
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
    fwd_color: str | None = None
    """
    Color of `+` arrows.
    """
    rev_color: str | None = None
    """
    Color of `-` arrows.
    """
    live_only: bool = True
    """
    Only plot live HORs.
    """
    mer_filter: int = 2
    """
    Filter HORs that have at least 2 monomers.
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
    mode: Literal["mer", "hor"] = "mer"
    """
    Plot HORs with `mer` or `hor`.
    """
    live_only: bool = True
    """
    Only plot live HORs.
    """
    mer_filter: int = 2
    """
    Filter HORs that have at least 2 monomers.
    """
    use_item_rgb: bool = False
    """
    Use `item_rgb` column for color.
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
