from typing import Literal
from dataclasses import dataclass
from ..track.types import LegendPosition

OutputFormat = Literal["png", "pdf", "svg"]


@dataclass
class PlotSettings:
    """
    Plot settings for a single plot.
    """

    title: str | None = None
    """
    Figure title. ex. "{chrom}"
    """
    format: list[OutputFormat] | OutputFormat = "png"
    """
    Output format(s). Either `"pdf"`, `"png"`, or `"svg"`.
    """
    transparent: bool = True
    """
    Output a transparent image.
    """
    dim: tuple[float, float] = (20.0, 12.0)
    """
    The dimensions of each plot.
    """
    dpi: int = 600
    """
    Set the plot DPI per plot.
    """
    layout: str = "tight"
    """
    Layout engine option for matplotlib. See https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.figure.html#matplotlib.pyplot.figure.
    """
    legend_pos: LegendPosition = LegendPosition.Right
    """
    Legend position as `LegendPosition`. Either `LegendPosition.Right` or `LegendPosition.Left`.
    """
    legend_prop: float = 0.2
    """
    Legend proportion of plot.
    """
    axis_h_pad: float = 0.2
    """
    Apply a height padding to each axis.
    """
    xlim: tuple[int, int] | None = None
    """
    Set x-axis limit across all plots.
    * `None` - Use the min and max position across all tracks.
    * `tuple[float, float]` - Use provided coordinates as min and max position.
    """
