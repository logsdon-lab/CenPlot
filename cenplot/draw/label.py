import ast

from typing import Any
from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Rectangle

from .utils import draw_uniq_entry_legend, format_ax
from ..track.types import Track


def draw_label(
    ax: Axes,
    track: Track,
    *,
    zorder: float,
    legend_ax: Axes | None = None,
) -> None:
    hide_x = track.options.hide_x
    color = track.options.color
    alpha = track.options.alpha
    legend = track.options.legend

    patch_options: dict[str, Any] = {"zorder": zorder}
    patch_options["alpha"] = alpha

    # Convert colors from rgb str -> rgb tuple -> hex
    track_color_mapping = {
        name: rgb2hex([c / 255 for c in ast.literal_eval(rgb)])
        for name, rgb in track.data.select("name", "item_rgb").unique().rows()
    }

    spines = ("right", "left", "top", "bottom") if hide_x else ("right", "left", "top")
    format_ax(
        ax,
        xticks=hide_x,
        xticklabel_fontsize=track.options.fontsize,
        yticks=True,
        yticklabel_fontsize=track.options.fontsize,
        spines=spines,
    )

    ylim = ax.get_ylim()
    height = ylim[1] - ylim[0]

    for row in track.data.iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        lbl_color = track_color_mapping.get(row["name"])

        if row["name"] == "-" or not row["name"]:
            labels = {}
        else:
            labels = {"label": row["name"]}
        # Allow override.
        if color:
            patch_options["color"] = color
        else:
            patch_options["color"] = lbl_color

        rect = Rectangle(
            (start, 0),
            end + 1 - start,
            height,
            lw=0,
            **labels,
            **patch_options,
        )
        ax.add_patch(rect)

    # Draw legend.
    if legend_ax and legend:
        draw_uniq_entry_legend(legend_ax, track, ref_ax=ax, loc="center left", ncols=3)
