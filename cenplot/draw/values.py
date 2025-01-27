from matplotlib.axes import Axes

from .utils import draw_uniq_entry_legend, minimalize_ax
from ..track import Track


def draw_values(
    ax: Axes,
    track: Track,
    *,
    color: str | None = None,
    alpha: float | None = None,
    legend_ax: Axes | None = None,
    zorder: float,
    hide_x: bool,
) -> None:
    plot_options = {"color": "black", "zorder": zorder}
    if color:
        plot_options["color"] = color
    if alpha:
        plot_options["alpha"] = alpha

    # Add line
    ax.plot(track.data["chrom_st"], track.data["name"], **plot_options)
    # Fill in-between
    ax.fill_between(track.data["chrom_st"], track.data["name"], **plot_options)
    # Trim plot to margins
    ax.margins(x=0, y=0)
    # remove spines
    minimalize_ax(ax, xticks=hide_x, spines=("right", "top"))

    # Limit spine range.
    # TODO: Remove ticks not within bounds.
    ax.spines["bottom"].set_bounds(0, track.data["chrom_end"].max())

    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center left", ncols=3)
