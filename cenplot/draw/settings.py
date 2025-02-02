from typing import Literal
from dataclasses import dataclass
from ..track.types import LegendPosition


@dataclass
class SinglePlotSettings:
    title: str | None = None
    """
    Figure title. ex. "{chrom}"
    """

    format: Literal["png", "pdf"] = "png"
    """
    Output format. Either `"pdf"` or `"png"`.
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
    shared_xlim: tuple[int, int] | None = None
    """
    Share x-axis limit across all plots.
    * `None` - Use the min and max position across all tracks.
    * `tuple[float, float]` - Use provided coordinates as min and max position.
    """
