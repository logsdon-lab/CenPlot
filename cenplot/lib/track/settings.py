from dataclasses import dataclass
from typing import Literal


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

    ymin: int = 0
    """
    Minimum y-value.
    """

    ymax: int | None = None
    """
    Maximum y-value.
    """

    label: str | None = None
    """
    Label to add to legend.
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
    bp_merge_units: int | None = 256
    """
    Merge HOR units into HOR blocks within this number of base pairs.
    """
    bp_merge_blks: int | None = 8000
    """
    Merge HOR blocks into HOR arrays within this number of bases pairs.
    """
    min_blk_hor_units: int | None = 2
    """
    Grouped stv rows must have at least `n` HOR units unbroken.
    """
    min_arr_hor_units: int | None = 10
    """
    Require that an HOR array have at least `n` HOR units.
    """
    min_arr_len: int | None = 30_000
    """
    Require that an HOR array is this size in bp.
    """
    min_arr_prop: float | None = 0.9
    """
    Require that an HOR array has at least this proportion of HORs by length.
    """


@dataclass
class HORTrackSettings(DefaultTrackSettings):
    """
    Higher order repeat plot options.
    """

    sort_order: Literal["ascending", "descending"] = "descending"
    """
    Plot HORs by `{mode}` in `{sort_order}` order.

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
