from typing import Literal
from dataclasses import dataclass

from ..defaults import Colorscale


@dataclass
class DefaultTrackSettings:
    """
    Default plot options settings.
    """

    fontsize: float | str | None = "medium"
    """
    Font size for track text.
    """

    title_fontsize: float | str | None = "large"
    """
    Font size for track title.
    """

    legend: bool = True
    """
    Show the legend.
    """

    legend_ncols: int | None = None
    """
    Number of columns for legend entries.
    """

    legend_fontsize: float | str | None = "medium"
    """
    Legend font size.
    """

    legend_title: str | None = None
    """
    Set legend title.
    * ex. "HOR monomers for {chrom}"
    """

    legend_title_fontsize: str = "large"
    """
    Legend title font size.
    """

    legend_title_only: bool = False
    """
    Hide all legend elements except titile.
    """

    legend_label_order: list[str] | None = None
    """
    Legend label order.
    """

    hide_x: bool = True
    """
    Hide the x-axis ticks, ticklabels, and spines.
    """

    units_x: Literal["bp", "kbp", "mbp"] = "mbp"
    """
    Set x-axis units.
    """


@dataclass
class SelfIdentTrackSettings(DefaultTrackSettings):
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
    colorscale: Colorscale | str | None = None
    """
    Colorscale for identity as TSV file.
    * Format: `[start, end, color]`
        * Color is a `str` representing a color name or hexcode.
        * See https://matplotlib.org/stable/users/explain/colors/colors.html
    * ex. `0\t90\tblue`
    """


@dataclass
class LabelTrackSettings(DefaultTrackSettings):
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

    use_item_rgb: bool = True
    """
    Use `item_rgb` column if provided. Otherwise, generate a random color for each value in column `name`.
    """

    alpha: float = 1.0
    """
    Label alpha.
    """

    shape: Literal["rect", "tri"] = "rect"
    """
    Shape to draw.
    * `"tri"` Always pointed down.
    """

    edgecolor: str | None = None
    """
    Edge color for each label.
    """

    bg_border: bool = False
    """
    Add black border containing all added labels.
    """


@dataclass
class LocalSelfIdentTrackSettings(LabelTrackSettings):
    """
    Local self-identity plot options.
    """

    colorscale: Colorscale | str | None = None
    """
    Colorscale for identity as TSV file.
    * Format: `[start, end, color]`
        * Color is a `str` representing a color name or hexcode.
        * See https://matplotlib.org/stable/users/explain/colors/colors.html
    * ex. `0\t90\tblue`
    """
    band_size: int = 5
    """
    Number of windows to calculate average sequence identity over.
    """
    ignore_band_size: int = 2
    """
    Number of windows ignored along self-identity diagonal.
    """


@dataclass
class BarTrackSettings(DefaultTrackSettings):
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

    ymin: int | Literal["min"] = 0
    """
    Minimum y-value.
    * Static value
    * 'min' for minimum value in data.
    """

    ymin_add: float = 0.0
    """
    Add some percent of y-axis minimum to y-axis limit.
    * ex. -0.05 subtracts 5% of min value so points aren't cutoff in plot.
    """

    ymax: int | Literal["max"] | None = None
    """
    Maximum y-value.
    * Static value
    * 'max' for maximum value in data.
    """

    ymax_add: float = 0.0
    """
    Add some percent of y-axis maximum to y-axis limit.
    * ex. 0.05 adds 5% of max value so points aren't cutoff in plot.
    """

    label: str | None = None
    """
    Label to add to legend.
    """

    add_end_yticks: bool = True
    """
    Add y-ticks showing beginning and end of data range.
    """


@dataclass
class LineTrackSettings(BarTrackSettings):
    """
    Line plot options.
    """

    position: Literal["start", "midpoint"] = "start"
    """
    Draw position at start or midpoint of interval.
    """
    fill: bool = False
    """
    Fill under line.
    """
    linestyle: str = "solid"
    """
    Line style. See https://matplotlib.org/stable/gallery/lines_bars_and_markers/linestyles.html.
    """
    linewidth: int | None = None
    """
    Line width.
    """
    marker: str | None = None
    """
    Marker shape. See https://matplotlib.org/stable/api/markers_api.html#module-matplotlib.markers,
    """
    markersize: int | None = None
    """
    Marker size.
    """
    log_scale: bool = False
    """
    Use log-scale for plot.
    """


@dataclass
class StrandTrackSettings(DefaultTrackSettings):
    """
    Strand arrow plot options.
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
    use_item_rgb: bool = False
    """
    Use `item_rgb` column if provided. Otherwise, use `fwd_color` and `rev_color`.
    """


@dataclass
class HOROrtTrackSettings(StrandTrackSettings):
    """
    Higher order repeat orientation arrow plot options.
    """

    live_only: bool = True
    """
    Only plot live HORs.
    """
    mer_filter: int = 2
    """
    Filter HORs that have at least 2 monomers.
    """
    arr_opt_bp_merge_units: int | None = 256
    """
    Merge HOR units into HOR blocks within this number of base pairs.
    """
    arr_opt_bp_merge_blks: int | None = 8000
    """
    Merge HOR blocks into HOR arrays within this number of bases pairs.
    """
    arr_opt_min_blk_hor_units: int | None = 2
    """
    Grouped stv rows must have at least `n` HOR units unbroken.
    """
    arr_opt_min_arr_hor_units: int | None = 10
    """
    hor_len_Require that a HOR array have at least `n` HOR units.
    """
    arr_opt_min_arr_len: int | None = 30_000
    """
    Require that a HOR array is this size in bp.
    """
    arr_opt_min_arr_prop: float | None = 0.9
    """
    Require that a HOR array has at least this proportion of HORs by length.
    """


@dataclass
class HORTrackSettings(DefaultTrackSettings):
    """
    Higher order repeat plot options.
    """

    sort_order: str = "descending"
    """
    Plot HORs by `{mode}` in `{sort_order}` order.

    Either:
    * `ascending`
    * `descending`
    * Or a path to a single column file specifying the order of elements of `mode`. Only for split.

    Mode:
    * If `{mer}`, sort by `mer` number
    * If `{hor}`, sort by `hor` frequency.
    """
    mode: Literal["mer", "hor"] = "mer"
    """
    Plot HORs with `mer` or `hor`.
    """
    live_only: bool = True
    """
    Only plot live HORs. Filters only for rows with `L` character in `name` column.
    """
    mer_size: int = 171
    """
    Monomer size to calculate number of monomers for mer_filter.
    """
    mer_filter: int = 2
    """
    Filter HORs that have less than `mer_filter` monomers.
    """
    hor_filter: int = 5
    """
    Filter HORs that occur less than `hor_filter` times.
    """
    color_map_file: str | None = None
    """
    Monomer color map TSV file. Two column headerless file that has `mode` to `color` mapping.
    """
    use_item_rgb: bool = False
    """
    Use `item_rgb` column for color. If omitted, use default mode color map or `color_map`.
    """
    split_prop: bool = False
    """
    If split, divide proportion evenly across each split track.
    """
    split_top_n: int | None = None
    """
    If split, show top n HORs for a given mode.
    """

    split_fill_missing: str | None = None
    """
    If split and defined sort order provided, fill in missing with this color. Otherwise, display random HOR variant.
    * Useful to maintain order across multiple plots.
    """

    split_sort_order_only: bool = False
    """
    If split and defined sort order provided, only show HORs within defined list.
    """

    bg_border: bool = False
    """
    Add black border containing all added labels.
    """

    bg_color: str | None = None
    """
    Background color for track.
    """


@dataclass
class LegendTrackSettings(DefaultTrackSettings):
    index: int | None = None
    """
    Index of plot to get legend of.
    """


@dataclass
class PositionTrackSettings(DefaultTrackSettings):
    pass


@dataclass
class SpacerTrackSettings(DefaultTrackSettings):
    pass


TrackSettings = (
    HORTrackSettings
    | HOROrtTrackSettings
    | SelfIdentTrackSettings
    | LocalSelfIdentTrackSettings
    | BarTrackSettings
    | LabelTrackSettings
    | LegendTrackSettings
    | PositionTrackSettings
    | SpacerTrackSettings
    | StrandTrackSettings
    | LineTrackSettings
)
"""
Type annotation for all possible settings for the various plot types.
"""
