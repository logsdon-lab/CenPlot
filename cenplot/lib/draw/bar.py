from matplotlib.axes import Axes


from .utils import draw_uniq_entry_legend, format_ax
from ..track.types import Track, TrackPosition


def draw_bars(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
) -> None:
    """
    Draw bar plot on axis with the given `Track`.
    """
    hide_x = track.options.hide_x
    color = track.options.color
    alpha = track.options.alpha
    legend = track.options.legend
    ymin = track.options.ymin
    ymax = track.options.ymax
    label = track.options.label

    if track.pos != TrackPosition.Overlap:
        spines = ("right", "top")
    else:
        spines = None

    format_ax(
        ax,
        xticks=hide_x,
        xticklabel_fontsize=track.options.fontsize,
        yticklabel_fontsize=track.options.fontsize,
        spines=spines,
    )

    plot_options = {"zorder": zorder, "alpha": alpha}
    if color:
        plot_options["color"] = color
    elif "color" in track.data.columns:
        plot_options["color"] = track.data["color"]
    else:
        plot_options["color"] = track.options.DEF_COLOR

    # Add bar
    ax.bar(
        track.data["chrom_st"],
        track.data["name"],
        track.data["chrom_end"] - track.data["chrom_st"],
        label=label,
        **plot_options,
    )
    # Trim plot to margins
    ax.margins(x=0, y=0)
    ax.set_ylim(ymin=ymin, ymax=ymax)

    if legend_ax and legend:
        draw_uniq_entry_legend(
            legend_ax,
            track,
            ref_ax=ax,
            ncols=track.options.legend_ncols,
            loc="center left",
            alignment="left",
        )
