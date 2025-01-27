import ast

from typing import Any
from matplotlib.axes import Axes
from matplotlib.colors import rgb2hex
from matplotlib.patches import Rectangle

from .utils import draw_uniq_entry_legend, minimalize_ax
from ..track import Track


def draw_label(
    ax: Axes,
    track: Track,
    *,
    height: int = 12,
    color: str | None = None,
    alpha: float | None = None,
    legend_ax: Axes | None = None,
    hide_x: bool,
    zorder: float,
) -> None:
    patch_options: dict[str, Any] = {"zorder": zorder}
    if alpha:
        patch_options["alpha"] = alpha

    # Convert colors from rgb str -> rgb tuple -> hex
    track_color_mapping = {
        name: rgb2hex([c / 255 for c in ast.literal_eval(rgb)])
        for name, rgb in track.data.select("name", "item_rgb").unique().rows()
    }

    minimalize_ax(
        ax, xticks=hide_x, yticks=True, spines=("left", "right", "top", "bottom")
    )

    for row in track.data.iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        # G
        lbl_color = track_color_mapping.get(row["name"], "red")

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
    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center left", ncols=3)
