import ast
import math
import multiprocessing
import os
import sys
import argparse
import itertools
import tomllib

from enum import StrEnum, auto
from collections import defaultdict
from typing import Any, Generator, NamedTuple
from concurrent.futures import ProcessPoolExecutor

import matplotlib
import matplotlib.axes
import matplotlib.colors
import matplotlib.figure
import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from intervaltree import Interval, IntervalTree
from matplotlib.collections import PolyCollection
from matplotlib.backends.backend_pdf import PdfPages

MONOMER_LEN = 170
MONOMER_COLORS = {
    "1": "#A8275C",
    "10": "#9AC78A",
    "11": "#A53D63",
    "12": "#3997C6",
    "13": "#29A3CE",
    "14": "#5EB2A7",
    "15": "#38A49B",
    "16": "#45B4CE",
    "17": "#A53D63",
    "18": "#AA1B63",
    "19": "#3F66A0",
    "2": "#D66C54",
    "20": "#BFDD97",
    "21": "#C0D875",
    "22": "#E5E57A",
    "24": "#B75361",
    "26": "#F9E193",
    "3": "#C6625D",
    "30": "#E5D1A1",
    "32": "#A1B5E5",
    "34": "#9F68A5",
    "35": "#81B25B",
    "4": "#F4DC78",
    "5": "#7EC0B3",
    "6": "#A8275C",
    "7": "#8CC49F",
    "8": "#893F89",
    "9": "#6565AA",
}
BED9_COLS = [
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "strand",
    "thick_st",
    "thick_end",
    "item_rgb",
]
BED9_COL_MAP = dict(
    zip(
        [f"column_{i}" for i in range(1, len(BED9_COLS) + 1)],
        BED9_COLS,
    )
)
BED_SELF_IDENT_COLS = [
    "query",
    "query_st",
    "query_end",
    "ref",
    "ref_st",
    "ref_end",
    "percent_identity_by_events",
]


IDENT_CUTOFF = 97.5
IDENT_INCREMENT = 0.25
IDENT_COLORS = [
    "#4b3991",
    "#2974af",
    "#4a9da8",
    "#57b894",
    "#9dd893",
    "#e1f686",
    "#ffffb2",
    "#fdda79",
    "#fb9e4f",
    "#ee5634",
    "#c9273e",
    "#8a0033",
]

ident_range_ends = tuple(
    # Increment by IDENT_INCREMENT from IDENT_CUTOFF up to 100%
    i + IDENT_CUTOFF
    for i in itertools.accumulate(IDENT_INCREMENT for _ in range(10))
)
IDENT_RANGE = [
    (0, 90),
    (90, 97.5),
    *zip((IDENT_CUTOFF, *ident_range_ends[:-1]), ident_range_ends),
]


IDENT_COLOR_RANGE = dict(zip(IDENT_RANGE, IDENT_COLORS))
WIDTH = 14


class TrackPosition(StrEnum):
    Overlap = auto()
    Relative = auto()


class TrackOption(StrEnum):
    HOR = auto()
    HORSplit = auto()
    Label = auto()
    Value = auto()
    SelfIdent = auto()


class Track(NamedTuple):
    name: str
    pos: TrackPosition
    opt: TrackOption
    prop: float
    data: pl.DataFrame
    options: dict[str, Any]


class Unit(StrEnum):
    Bp = auto()
    Kbp = auto()
    Mbp = auto()

    def convert_value(self, value: int | float, round_to: int = 3) -> float:
        if self == Unit.Bp:
            new_value = value
        elif self == Unit.Kbp:
            new_value = value / 1_000
        else:
            new_value = value / 1_000_000

        return round(new_value, round_to)


def read_bed_hor(
    infile: str, *, chrom: str | None = None, mer_order: str = "large"
) -> pl.DataFrame:
    return (
        read_bed9(infile, chrom=chrom)
        .with_columns(
            length=pl.col("chrom_end") - pl.col("chrom_st"),
        )
        .with_columns(mer=(pl.col("length") / MONOMER_LEN).round().cast(pl.Int8))
        .filter(
            pl.when(pl.col("chrom_name") == "chr10")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chr20")
            .then(pl.col("mer") >= 5)
            .when(pl.col("chrom_name") == "chrY")
            .then(pl.col("mer") >= 30)
            .when(pl.col("chrom_name") == "chr17")
            .then(pl.col("mer") >= 4)
            .otherwise(True)
        )
        .sort("mer", descending=mer_order == "large")
        .cast({"mer": pl.String})
        # Add 1000 bp for visualization purposes.
        .with_columns(pl.col("length") + 2000)
    )


def get_stv_mon_ort(df_stv: pl.DataFrame, *, dst_merge: int) -> pl.DataFrame:
    stv_itrees: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    for grp, df_grp in df_stv.group_by(["chrom", "strand"]):
        chrom, strand = grp
        itree = IntervalTree.from_tuples(
            (row["chrom_st"] - dst_merge, row["chrom_end"] + dst_merge, strand)
            for row in df_grp.select("chrom_st", "chrom_end", "chrom").iter_rows(
                named=True
            )
        )
        # TODO: Change this to use custom merge function.
        itree.merge_overlaps(strict=False, data_reducer=lambda x, _: x)

        stv_itrees[chrom] = stv_itrees[chrom].union(
            # Restore original coords.
            IntervalTree(
                Interval(itv.begin + dst_merge, itv.end - dst_merge, itv.data)
                for itv in itree
            )
        )

    return pl.DataFrame(
        (
            (chrom, itv.begin, itv.end, itv.data)
            for chrom, itree in stv_itrees.items()
            for itv in itree.iter()
        ),
        orient="row",
        schema=["chrom", "chrom_st", "chrom_end", "strand"],
    ).sort(by="chrom_st")


def adj_by_ctg_coords(df: pl.DataFrame, colname: str) -> pl.DataFrame:
    final_adj_coords = {
        f"{colname}_st": pl.col(f"{colname}_st") - pl.col("ctg_st"),
        f"{colname}_end": pl.col(f"{colname}_end") - pl.col("ctg_st"),
    }
    return df.with_columns(
        chrom_name=pl.col(colname).str.extract(r"(chr[\dXY]+)").fill_null(""),
        # Use simplified coordinates if possible, otherwise, take everything.
        ctg_name=pl.col(colname).str.extract(r"^(.*?):|^(.*?)$"),
        ctg_st=pl.col(colname).str.extract(r":(\d+)-").cast(pl.Int64).fill_null(0),
        ctg_end=pl.col(colname).str.extract(r"-(\d+)$").cast(pl.Int64).fill_null(0),
    ).with_columns(
        # Use ctg name without coordinates. Avoids 0 and 1 offset issues.
        chrom=pl.col("ctg_name"),
        **final_adj_coords,
    )


def read_bed9(infile: str, *, chrom: str | None = None) -> pl.DataFrame:
    df = pl.read_csv(infile, separator="\t", has_header=False)
    df = df.rename({col: val for col, val in BED9_COL_MAP.items() if col in df.columns})
    df_adj = adj_by_ctg_coords(df, "chrom").sort(by="chrom_st")

    if chrom:
        df_adj = df_adj.filter(pl.col("chrom_name") == chrom)
    if "item_rgb" not in df.columns:
        df_adj = df_adj.with_columns(item_rgb=pl.lit("0,0,0"))
    if "name" not in df.columns:
        df_adj = df_adj.with_columns(name=pl.lit("-"))

    return df_adj


def read_bed_identity(infile: str, *, chrom: str | None = None) -> pl.DataFrame:
    df = pl.read_csv(
        infile, separator="\t", has_header=False, new_columns=BED_SELF_IDENT_COLS
    ).with_columns(
        query_name=pl.col("query").str.extract(r"^(.*?):|^(.*?)$"),
        chrom_name=pl.col("query").str.extract(r"(chr[\dXY]+)").fill_null(""),
    )
    if chrom:
        df = df.filter(pl.col("chrom_name") == chrom)

    # Build expr to filter range of colors.
    # TODO: Allow custom range.
    color_expr = None
    for rng, color in IDENT_COLOR_RANGE.items():
        if not isinstance(color_expr, pl.Expr):
            color_expr = pl.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))
        else:
            color_expr = color_expr.when(
                pl.col("percent_identity_by_events").is_between(rng[0], rng[1])
            ).then(pl.lit(color))

    if isinstance(color_expr, pl.Expr):
        color_expr = color_expr.otherwise(None)
    else:
        color_expr = pl.lit(None)

    tri_side = math.sqrt(2) / 2
    df_tri = (
        df.lazy()
        .with_columns(color=color_expr)
        # Get window size.
        .with_columns(
            window=(pl.col("query_end") - pl.col("query_st")).max().over("query_name")
        )
        .with_columns(
            first_pos=pl.col("query_st") // pl.col("window"),
            second_pos=pl.col("ref_st") // pl.col("window"),
        )
        # x y coords of diamond
        .with_columns(
            x=pl.col("first_pos") + pl.col("second_pos"),
            y=-pl.col("first_pos") + pl.col("second_pos"),
        )
        .with_columns(
            scale=(pl.col("query_st").max() / pl.col("x").max()).over("query_name"),
            group=pl.int_range(pl.len()).over("query_name"),
        )
        .with_columns(
            window=pl.col("window") / pl.col("scale"),
        )
        # Rather than generate new dfs. Add new x,y as arrays per row.
        .with_columns(
            new_x=[tri_side, 0.0, -tri_side, 0.0],
            new_y=[0.0, tri_side, 0.0, -tri_side],
        )
        # Convert to array to save memory
        .with_columns(
            pl.col("new_x").list.to_array(4), pl.col("new_y").list.to_array(4)
        )
        # Rescale x and y.
        .with_columns(
            ((pl.col("new_x") * pl.col("window")) + pl.col("x")) * pl.col("scale"),
            ((pl.col("new_y") * pl.col("window")) + pl.col("y")) * pl.col("window"),
        )
        .select(
            "query_name",
            "new_x",
            "new_y",
            "color",
            "group",
            "percent_identity_by_events",
        )
        # arr to new rows
        .explode("new_x", "new_y")
        # Rename to filter later on.
        .rename({"query_name": "chrom", "new_x": "x", "new_y": "y"})
        .collect()
    )
    return df_tri


def read_bed_label(infile: str, *, chrom: str | None = None) -> pl.DataFrame:
    df_track = read_bed9(infile, chrom=chrom)

    # Order facets by descending length. This prevents larger annotations from blocking others.
    fct_name_order = (
        df_track.group_by(["name"])
        .agg(len=(pl.col("chrom_end") - pl.col("chrom_st")).sum())
        .sort(by="len", descending=True)
        .get_column("name")
    )
    return df_track.cast({"name": pl.Enum(fct_name_order)})


def read_one_track_info(
    track: dict[str, Any], *, chrom: str | None = None
) -> Generator[Track, None, None]:
    prop = track.get("proportion", 0.0)
    name = track.get("name")
    pos = track.get("position")
    opt = track.get("type")
    path: str | None = track.get("path")
    options: dict[str, Any] = track.get("options", {})

    if not name:
        raise ValueError(f"Name not provided for track ({track}).")

    try:
        track_pos = TrackPosition(pos)  # type: ignore[arg-type]
    except ValueError:
        print(
            f"Invalid plot position ({pos}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None
    try:
        track_opt = TrackOption(opt)  # type: ignore[arg-type]
    except ValueError:
        print(
            f"Invalid plot option ({opt}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None

    if not path:
        raise ValueError("Path to data required.")

    if not os.path.exists(path):
        raise FileNotFoundError(f"Data does not exist for track ({track})")

    if opt == TrackOption.HORSplit:
        mer_order = options.get("mer_order", "large")
        df_track = read_bed_hor(path, chrom=chrom, mer_order=mer_order)
        uniq_mers = df_track["mer"].unique()
        track_prop = prop / len(uniq_mers)
        if track_pos == TrackPosition.Overlap:
            print(
                f"Overlap not supported for {opt}. Using relative position.",
                file=sys.stderr,
            )

        for mer, df_mer_track in df_track.group_by(["mer"], maintain_order=True):
            mer = mer[0]
            # Add (mer) to name.
            # Disallow overlap.
            # Split proportion over uniq monomers.
            yield Track(
                f"{name} ({mer})",
                TrackPosition.Relative,
                TrackOption.HOR,
                track_prop,
                df_mer_track,
                options,
            )

        return None

    if opt == TrackOption.HOR:
        mer_order = options.get("mer_order", "large")
        df_track = read_bed_hor(path, chrom=chrom, mer_order=mer_order)
    elif opt == TrackOption.SelfIdent:
        df_track = read_bed_identity(path, chrom=chrom)
    elif opt == TrackOption.Value:
        df_track = read_bed9(path, chrom=chrom)
    else:
        df_track = read_bed_label(path, chrom=chrom)

    yield Track(name, track_pos, track_opt, prop, df_track, options)


def get_track_info(
    input_tracks: list[str], *, chrom: str | None = None
) -> tuple[list[Track], set[str]]:
    dfs = []
    chroms = set()
    for input_track in input_tracks:
        with open(input_track, "rb") as fh:
            tracks = tomllib.load(fh).get("tracks", {})["plots"]
            for track_info in tracks:
                for track in read_one_track_info(track_info, chrom=chrom):
                    dfs.append(track)
                    chroms.update(track.data["chrom"])

    return dfs, chroms


def draw_uniq_entry_legend(
    ax: matplotlib.axes.Axes, ref_ax: matplotlib.axes.Axes | None = None, **kwargs: Any
) -> None:
    ref_ax = ref_ax if ref_ax else ax
    handles, labels = ref_ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), **kwargs)


def minimalize_ax(
    ax: matplotlib.axes.Axes,
    *,
    grid=False,
    xticks: bool = False,
    yticks: bool = False,
    spines: tuple[str, ...] | None = None,
) -> None:
    if grid:
        ax.grid(False)
    if xticks:
        ax.set_xticks([], [])
    if yticks:
        ax.set_yticks([], [])
    if spines:
        for spine in spines:
            ax.spines[spine].set_visible(False)


def draw_label(
    ax: matplotlib.axes.Axes,
    track: Track,
    *,
    height: int = 12,
    color: str | None = None,
    alpha: float | None = None,
    legend_ax: matplotlib.axes.Axes | None = None,
    hide_x: bool,
    zorder: float,
) -> None:
    patch_options: dict[str, Any] = {"zorder": zorder}
    if alpha:
        patch_options["alpha"] = alpha

    # Convert colors from rgb str -> rgb tuple -> hex
    track_color_mapping = {
        name: matplotlib.colors.rgb2hex([c / 255 for c in ast.literal_eval(rgb)])
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

        rect = mpatches.Rectangle(
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
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center right", ncols=3)


def draw_self_ident(
    ax: matplotlib.axes.Axes,
    track: Track,
    *,
    legend_ax: matplotlib.axes.Axes | None = None,
    zorder: float,
    hide_x: bool,
    flip_y: bool,
    bins: int = 300,
) -> None:
    colors, verts = [], []

    minimalize_ax(
        ax, xticks=hide_x, yticks=True, spines=("left", "right", "top", "bottom")
    )

    if flip_y:
        df_track = track.data.with_columns(y=-pl.col("y"))
    else:
        df_track = track.data

    for _, df_diam in df_track.group_by(["group"]):
        df_points = df_diam.select("x", "y")
        color = df_diam["color"].first()
        colors.append(color)
        verts.append(df_points)

    # https://stackoverflow.com/a/29000246
    polys = PolyCollection(verts, zorder=zorder)
    polys.set(array=None, facecolors=colors)
    ax.add_collection(polys)

    ax.set_ylim(df_track["y"].min(), df_track["y"].max())

    if legend_ax:
        cmap = IntervalTree(
            Interval(rng[0], rng[1], color) for rng, color in IDENT_COLOR_RANGE.items()
        )
        cnts, values, bars = legend_ax.hist(
            track.data["percent_identity_by_events"], bins=bins, zorder=zorder
        )
        legend_ax.set_xlim(70.0, 100.0)
        legend_ax.set_xlabel("Average Nucleotide Identity (%)")
        legend_ax.set_ylabel("Number of Intervals (k)")

        for _, value, bar in zip(cnts, values, bars):
            color = cmap.overlap(value, value + 0.00001)
            try:
                color = next(iter(color)).data
            except Exception:
                color = None
            bar.set_facecolor(color)


def draw_values(
    ax: matplotlib.axes.Axes,
    track: Track,
    *,
    color: str | None = None,
    alpha: float | None = None,
    legend_ax: matplotlib.axes.Axes | None = None,
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
        draw_uniq_entry_legend(
            legend_ax, ref_ax=ax, loc="center right", bbox_to_anchor=[0.5, 0], ncols=3
        )


def draw_hor_w_ort(
    ax: matplotlib.axes.Axes,
    df_stv: pl.DataFrame,
    df_stv_ort: pl.DataFrame,
    *,
    hide_x: bool,
    ort: bool,
    ort_pos: str,
    zorder: float,
    legend_ax: matplotlib.axes.Axes | None,
):
    minimalize_ax(
        ax, xticks=hide_x, yticks=True, spines=("right", "left", "top", "bottom")
    )

    # Add HOR track.
    for row in df_stv.with_row_index().iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        row_idx = row["index"]
        # TODO: Adjust zorder so mer order is respected.
        zorder_adj = zorder + (row_idx / df_stv.shape[0])
        color = MONOMER_COLORS.get(row["mer"])
        rect = mpatches.Rectangle(
            (start, 0),
            end + 1 - start,
            1,
            color=color,
            lw=0,
            label=row["mer"],
            zorder=zorder_adj,
        )
        ax.add_patch(rect)

    # Add ort track.
    if ort_pos == "top":
        y_ort = 1.2
    else:
        y_ort = -0.1
    if ort:
        for row in df_stv_ort.iter_rows(named=True):
            # sample arrow
            start = row["chrom_st"]
            end = row["chrom_end"]
            strand = row["strand"]
            length = end - start
            # Skip single orts
            if length < 171:
                continue

            if strand == "-":
                start = end
                length = -length

            arrow = mpatches.FancyArrow(
                start,
                y_ort,
                length,
                0,
                width=0.05,
                color="black",
                length_includes_head=True,
                # 10% of entire length
                head_length=0.05 * -length,
                clip_on=False,
                zorder=zorder,
            )
            ax.add_patch(arrow)

    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center right", ncols=3)


def create_subplots(
    dfs_track: list[Track], width: float, height: float, **kwargs: Any
) -> tuple[matplotlib.figure.Figure, np.ndarray, dict[str, int]]:
    track_props = []
    track_indices = {}
    tracks_added = set()
    requires_second_col = False

    track_idx = 0
    for track in dfs_track:
        if track.name in tracks_added:
            print(f"Skipping duplicate {track.name} track.", file=sys.stderr)
            continue

        # Store index.
        # Only increment index if takes up a subplot axis.
        if track.pos == TrackPosition.Relative:
            track_indices[track.name] = track_idx
            track_idx += 1
            track_props.append(track.prop)
        # For each unique HOR monomer number, create a new track.
        # Divide the proportion of the image allocated between each mer track.
        elif track.opt == TrackOption.HORSplit:
            uniq_mers = track.data["mer"].unique()
            track_prop = track.prop / len(uniq_mers)
            for mer in uniq_mers:
                track_indices[f"{track.name}_{mer}"] = track_idx
                track_props.append(track_prop)
                track_idx += 1
        else:
            track_indices[track.name] = track_idx - 1

        if not requires_second_col and track.options.get("legend"):
            requires_second_col = True

        tracks_added.add(track.name)

    # Adjust columns and width ratio.
    num_cols = 2 if requires_second_col else 1
    width_ratios = (0.8, 0.2) if requires_second_col else [1.0]
    fig, axes = plt.subplots(
        # Count number of tracks
        len(track_props),
        num_cols,
        figsize=(width, height),
        height_ratios=track_props,
        width_ratios=width_ratios,
        # Always return 2D ndarray
        squeeze=0,
        **kwargs,
    )

    return fig, axes, track_indices


def plot_one_cen(
    dfs_track: list[Track],
    outdir: str,
    outfmt: str,
    chrom: str,
    min_st_pos: int,
    max_end_pos: int,
    width: float,
    height: float,
) -> tuple[matplotlib.figure.Figure, np.ndarray, str]:
    print(f"Plotting {chrom}...", file=sys.stderr)

    # # Get min and max position of all tracks for this cen.
    #  = sys.maxsize
    # trk_maxtrk_min_st_end = 0
    # for trk in dfs_track:
    #     if trk.opt == TrackOption.SelfIdent:
    #         continue
    #     try:
    #         trk_min_st = min(trk.data["chrom_st"].min(), trk_min_st)
    #         trk_max_end = max(trk.data["chrom_end"].max(), trk_max_end)
    #     except TypeError:
    #         continue

    # # Scale height based on track length.
    # adj_height = height * (trk_max_end / max_end_pos)
    # height = height if adj_height == 0 else adj_height

    fig, axes, track_indices = create_subplots(
        dfs_track, width, height, tight_layout=True
    )
    track_col, legend_col = 0, 1

    track_labels: list[str] = []

    def get_track_label(track: Track, all_track_labels: list[str]) -> tuple[str, str]:
        esc_track_label = track.name.encode("unicode_escape").decode("utf-8")
        track_label = track.name.encode("ascii", "ignore").decode("unicode_escape")

        # Update track label for each overlap.
        if track.pos == TrackPosition.Overlap:
            try:
                track_label = f"{all_track_labels[-1]}\n{track_label}"
            except IndexError:
                pass

        return esc_track_label, track_label

    for zorder, track in enumerate(dfs_track):
        track_row = track_indices[track.name]
        esc_track_label, track_label = get_track_label(track, track_labels)

        try:
            track_ax: matplotlib.axes.Axes = axes[track_row, track_col]
        except IndexError:
            print(
                f"Cannot get track ({track_row, track_col}) for {esc_track_label} with {track.pos} position."
            )
            continue
        try:
            legend_ax: matplotlib.axes.Axes = axes[track_row, legend_col]
        except IndexError:
            legend_ax = None

        # Set xaxis limits
        track_ax.set_xlim(min_st_pos, max_end_pos)

        # Switch to line if different track option. {value, label, ident}
        if track.opt == TrackOption.HOR:
            ort_dst_merge = track.options.get("ort_dst_merge", 1)
            df_stv_ort = get_stv_mon_ort(track.data, dst_merge=ort_dst_merge)
            draw_hor_w_ort(
                ax=track_ax,
                legend_ax=legend_ax if track.options.get("legend") else None,
                df_stv=track.data,
                df_stv_ort=df_stv_ort,
                zorder=zorder,
                hide_x=track.options.get("hide_x", False),
                ort=track.options.get("ort", True),
                ort_pos=track.options.get("ort_pos", "top"),
            )

        elif track.opt == TrackOption.Label:
            draw_label(
                track_ax,
                track,
                color=track.options.get("color"),
                alpha=track.options.get("alpha"),
                legend_ax=legend_ax if track.options.get("legend") else None,
                hide_x=track.options.get("hide_x", False),
                zorder=zorder,
            )

        elif track.opt == TrackOption.SelfIdent:
            draw_self_ident(
                track_ax,
                track,
                legend_ax=legend_ax if track.options.get("legend") else None,
                hide_x=track.options.get("hide_x", False),
                flip_y=track.options.get("flip_y", True),
                zorder=zorder,
            )

        elif track.opt == TrackOption.Value:
            draw_values(
                track_ax,
                track,
                color=track.options.get("color"),
                alpha=track.options.get("alpha"),
                zorder=zorder,
                hide_x=track.options.get("hide_x", False),
            )

        # Store label if more overlaps.
        track_labels.append(track_label)

        # Set label.
        # Allow chrom as title or name.
        track_label = chrom if track.options.get("chrom_as_title") else track_label
        if track.options.get("title", True):
            track_ax.set_ylabel(
                track_label,
                rotation="horizontal",
                ha="right",
                va="center",
                ma="center",
            )

        if not legend_ax:
            continue

        if track.opt != TrackOption.SelfIdent or (
            track.opt == TrackOption.SelfIdent and not track.options.get("legend")
        ):
            # Minimalize all legend cols
            minimalize_ax(
                legend_ax,
                grid=True,
                xticks=True,
                yticks=True,
                spines=("right", "left", "top", "bottom"),
            )
        else:
            minimalize_ax(
                legend_ax,
                grid=True,
                spines=("right", "top"),
            )

    # Add title
    # fig.suptitle(chrom)
    outfile = os.path.join(outdir, f"{chrom}.{outfmt}")
    plt.savefig(outfile, dpi=600, transparent=True)

    return fig, axes, outfile


def get_min_max_track(
    tracks: list[Track], typ: str, default_col: str = "chrom_st"
) -> tuple[Track, int]:
    track = None
    if typ == "min":
        pos = sys.maxsize
    else:
        pos = 0

    for trk in tracks:
        if trk.opt == TrackOption.SelfIdent:
            col = "x"
        else:
            col = default_col
        if typ == "min":
            trk_min = trk.data.filter(pl.col(col) > 0)[col].min()
            if trk_min < pos:
                track = trk
                pos = trk_min
        else:
            trk_max = trk.data[col].max()
            if trk_max > pos:
                track = trk
                pos = trk_max
    if not track:
        raise ValueError("No tracks.")
    return track, pos


def read_standardized_tracks(
    input_tracks: list[str], chrom: str | None = None
) -> tuple[list[Track], set[str], tuple[int, int]]:
    tracks, all_chroms = get_track_info(input_tracks, chrom=chrom)

    _, min_st_pos = get_min_max_track(tracks, typ="min")
    _, max_end_pos = get_min_max_track(tracks, typ="max", default_col="chrom_end")

    return tracks, all_chroms, (min_st_pos, max_end_pos)


def merge_plots(
    figures: list[tuple[matplotlib.figure.Figure, np.ndarray, str]], outfile: str
) -> None:
    if outfile.endswith(".pdf"):
        with PdfPages(outfile) as pdf:
            for fig, _, _ in figures:
                pdf.savefig(fig)
    else:
        merged_images = np.concatenate([plt.imread(file) for _, _, file in figures])
        plt.imsave(outfile, merged_images)


def get_inputs(
    args: argparse.Namespace,
) -> list[tuple[list[Track], str, str, str, int, int, int, int]]:
    tracks, all_chroms, (min_st, max_end) = read_standardized_tracks(args.input_tracks)
    if args.chroms:
        all_chroms = args.chroms
    return [
        (
            [
                Track(
                    trk.name,
                    trk.pos,
                    trk.opt,
                    trk.prop,
                    trk.data.filter(pl.col("chrom") == chrom),
                    trk.options,
                )
                for trk in tracks
            ],
            args.outdir,
            args.format,
            chrom,
            min_st,
            max_end,
            args.width,
            args.height,
        )
        for chrom in all_chroms
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-t",
        "--input_tracks",
        nargs="*",
        type=str,
        help=(
            "TOML file with headerless BED files to plot. "
            "Specify under tracks.plots the following fields: {name, position, type, proportion, path, or options}. "
            "One of more TOML can be provided."
        ),
    )
    ap.add_argument(
        "-c",
        "--chroms",
        nargs="*",
        help="Names to plot in this order. Corresponds to 1st col in BED files.",
        default=None,
    )
    ap.add_argument(
        "-d",
        "--outdir",
        help="Output dir to plot multiple separate figures.",
        type=str,
        default=".",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        help="Output file merging all figures. Either pdf of png.",
        type=str,
        default=None,
    )
    ap.add_argument(
        "-f",
        "--format",
        help="Output format passed as ext to matplotlib. Without '.'",
        default="png",
    )
    ap.add_argument(
        "-w",
        "--width",
        type=float,
        help="Figure width in inches per centromere.",
        default=16.0,
    )
    ap.add_argument(
        "-ht",
        "--height",
        type=float,
        help="Figure height in inches per centromere.",
        default=8.0,
    )

    ap.add_argument("-p", "--processes", type=int, default=4, help="Processes to run.")
    args = ap.parse_args()

    tracks = get_inputs(args)

    os.makedirs(args.outdir, exist_ok=True)
    if args.processes == 1:
        plots = [plot_one_cen(*track) for track in tracks]
    else:
        with ProcessPoolExecutor(
            max_workers=args.processes, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            plots = pool.map(plot_one_cen, *zip(*tracks))  # type: ignore[assignment]

    if args.outfile:
        merge_plots(plots, args.outfile)


if __name__ == "__main__":
    raise SystemExit(main())
