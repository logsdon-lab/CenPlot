from typing import Literal
from dataclasses import dataclass
from ..track.types import LegendPosition


@dataclass
class SinglePlotSettings:
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
    legend_pos: LegendPosition = LegendPosition.Right
    """
    Legend position as `LegendPosition`. Either `LegendPosition.Right` or `LegendPosition.Left`.
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
