import polars as pl
from matplotlib.axes import Axes


from .utils import draw_uniq_entry_legend, format_ax, set_ylim
from ..track.types import Track, TrackPosition


def draw_line(
    ax: Axes,
    track: Track,
    *,
    zorder: float = 1.0,
    legend_ax: Axes | None = None,
) -> None:
    """
    Draw line plot on axis with the given `Track`.
    """
    hide_x = track.options.hide_x
    color = track.options.color
    alpha = track.options.alpha
    legend = track.options.legend
    label = track.options.label
    linestyle = track.options.linestyle
    linewidth = track.options.linewidth
    marker = track.options.marker
    markersize = track.options.markersize

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

    if linestyle:
        plot_options["linestyle"] = linestyle
    if linewidth:
        plot_options["linewidth"] = linewidth

    # Fill between cannot add markers
    if not track.options.fill:
        plot_options["marker"] = marker
        if markersize:
            plot_options["markersize"] = markersize

    if track.options.position == "midpoint":
        df = track.data.with_columns(
            chrom_st=pl.col("chrom_st") + (pl.col("chrom_end") - pl.col("chrom_st")) / 2
        )
    else:
        df = track.data

    if track.options.log_scale:
        ax.set_yscale("log")

    # Add bar
    if track.options.fill:
        ax.fill_between(
            df["chrom_st"],
            df["name"],
            0,
            label=label,
            **plot_options,
        )
    else:
        ax.plot(
            df["chrom_st"],
            df["name"],
            label=label,
            **plot_options,
        )

    # Trim plot to margins
    ax.margins(x=0, y=0)

    set_ylim(ax, track)

    if legend_ax and legend:
        draw_uniq_entry_legend(
            legend_ax,
            track,
            ref_ax=ax,
            ncols=track.options.legend_ncols,
            loc="center left",
            alignment="left",
        )
