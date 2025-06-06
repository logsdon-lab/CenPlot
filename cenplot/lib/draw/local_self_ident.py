from matplotlib.axes import Axes

from cenplot.lib.draw.label import draw_label

from ..track.types import Track


def draw_local_self_ident(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
) -> None:
    """
    Draw local, self identity plot on axis with the given `Track`.
    """
    draw_label(ax, track, zorder=zorder, legend_ax=legend_ax)
