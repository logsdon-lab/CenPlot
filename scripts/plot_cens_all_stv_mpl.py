import ast
import math
import sys
import argparse
import itertools

from enum import StrEnum, auto
from collections import defaultdict, deque
from typing import Any, NamedTuple
# from concurrent.futures import ProcessPoolExecutor

import matplotlib
import matplotlib.axes
import matplotlib.axis
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import matplotlib.colors
import polars as pl
from intervaltree import Interval, IntervalTree

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
    Top = auto()
    Bottom = auto()
    Center = auto()


class TrackOption(StrEnum):
    Label = auto()
    Value = auto()
    SelfIdent = auto()


class Track(NamedTuple):
    name: str
    pos: TrackPosition
    opt: TrackOption
    prop: float
    data: pl.DataFrame


def read_stv(
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


def get_stv_mon_ort(df_stv: pl.DataFrame, *, dst_merge: int = 100_000) -> pl.DataFrame:
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


def read_bed_identity(infile: str, *, chrom: str | None = None):
    df = pl.read_csv(
        infile, separator="\t", has_header=False, new_columns=BED_SELF_IDENT_COLS
    ).with_columns(
        query_name=pl.col("query").str.extract(r"^(.*?):|^(.*?)$"),
        chrom_name=pl.col("query").str.extract(r"(chr[\dXY]+)").fill_null(""),
    )
    if chrom:
        df = df.filter(pl.col("chrom_name") == chrom)

    # Build expr to filter range of colors.
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

    tri_side = math.sqrt(2) / 2
    df_tri = (
        df.lazy()
        .with_columns(color=color_expr.otherwise(None))
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
        .select("query_name", "new_x", "new_y", "color", "group", "percent_identity_by_events")
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
    name: str, pos: str, opt: str, prop: str, path: str, *, chrom: str | None = None
) -> Track | None:
    prop = float(prop)

    try:
        pos = TrackPosition(pos)
    except ValueError:
        print(
            f"Invalid plot position ({pos}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None
    try:
        opt = TrackOption(opt)
    except ValueError:
        print(
            f"Invalid plot option ({pos}) for {path}. Skipping.",
            file=sys.stderr,
        )
        return None

    if opt == TrackOption.SelfIdent:
        df_track = read_bed_identity(path, chrom=chrom)
    elif opt == TrackOption.Value:
        df_track = read_bed9(path, chrom=chrom)
    else:
        df_track = read_bed_label(path, chrom=chrom)

    return Track(name, pos, opt, prop, df_track)


def get_track_info(input_tracks: list[str], *, chrom: str | None = None) -> list[Track]:
    dfs = []
    for input_track in input_tracks:
        track_info = input_track.split(",")
        num_fields = len(track_info)
        # Allow csv.
        if num_fields == 1 and input_track.endswith(".csv"):
            with open(input_track, "rt") as fh:
                for line in fh.readlines():
                    line = line.strip().split(",")
                    track = read_one_track_info(*line, chrom=chrom)
                    if not track:
                        continue
                    dfs.append(track)
        else:
            track = read_one_track_info(*track_info, chrom=chrom)
            if not track:
                continue
            dfs.append(track)

    return dfs


def standardize_coords(
    df_stv: pl.DataFrame, dfs_track: list[Track], df_min_st: pl.DataFrame
) -> tuple[pl.DataFrame, list[Track]]:
    # Get min across all groups.
    new_min = df_min_st["chrom_st"].min()
    df_new_stv = df_stv.join(
        df_min_st.group_by("chrom").agg(dst_diff=pl.col("chrom_st").min() - new_min),
        on="chrom",
    ).with_columns(
        chrom_st=pl.col("chrom_st") - pl.col("dst_diff") - new_min,
        chrom_end=pl.col("chrom_end") - pl.col("dst_diff") - new_min,
    )
    new_dfs_track = []
    for track in dfs_track:
        if track.opt == TrackOption.SelfIdent:
            df_adj_track = track.data
        else:
            df_adj_track = track.data.join(
                df_min_st.group_by("chrom").agg(
                    dst_diff=pl.col("chrom_st").min() - new_min
                ),
                on="chrom",
            ).with_columns(
                chrom_st=pl.col("chrom_st") - pl.col("dst_diff") - new_min,
                chrom_end=pl.col("chrom_end") - pl.col("dst_diff") - new_min,
            )
        new_dfs_track.append(
            Track(track.name, track.pos, track.opt, track.prop, df_adj_track)
        )
    return df_new_stv, new_dfs_track


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
    ax: matplotlib.axes.Axes, legend_ax: matplotlib.axes.Axes, track: Track
) -> None:
    # Convert colors from rgb str -> rgb tuple -> hex
    track_color_mapping = {
        name: matplotlib.colors.rgb2hex([c / 255 for c in ast.literal_eval(rgb)])
        for name, rgb in track.data.select("name", "item_rgb").unique().rows()
    }

    minimalize_ax(
        ax, xticks=True, yticks=True, spines=("left", "right", "top", "bottom")
    )

    for row in track.data.iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        color = track_color_mapping.get(row["name"], "red")
        if row["name"] == "-" or not row["name"]:
            labels = {}
        else:
            labels = {"label": row["name"]}
        rect = mpatches.Rectangle(
            (start, 0), end + 1 - start, 12, color=color, lw=0, **labels
        )
        ax.add_patch(rect)

    # Draw legend.
    draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center right", ncols=3)


def draw_self_ident(
    ax: matplotlib.axes.Axes,
    track: Track,
    *,
    legend_ax: matplotlib.axes.Axes | None = None,
    is_below_hor_track: bool,
) -> None:
    colors, verts = [], []

    minimalize_ax(ax, yticks=True, spines=("left", "right", "top"))

    # Reorient data so correct orient to HOR track.
    if is_below_hor_track:
        df_track = track.data.with_columns(y=-pl.col("y"))
    else:
        df_track = track.data
        ax.tick_params(top=False, labeltop=False, bottom=True, labelbottom=True)

    for _, df_diam in df_track.group_by(["group"]):
        df_points = df_diam.select("x", "y")
        color = df_diam["color"].first()
        colors.append(color)
        verts.append(df_points)

    # https://stackoverflow.com/a/29000246
    polys = PolyCollection(verts)
    polys.set(array=None, facecolors=colors)
    ax.add_collection(polys)

    ax.set_xlim(df_track["x"].min(), df_track["x"].max())
    ax.set_ylim(df_track["y"].min(), df_track["y"].max())

    if legend_ax:
        cmap = IntervalTree(Interval(rng[0], rng[1], color) for rng, color in IDENT_COLOR_RANGE.items())
        cnts, values, bars = legend_ax.hist(track.data["percent_identity_by_events"], bins=300)
        legend_ax.set_xlim(70.0, 100.0)
        legend_ax.set_xlabel("Average Nucleotide Identity (%)")
        legend_ax.set_ylabel("Number of Intervals (k)")

        for _, value, bar in zip(cnts, values, bars):
            color = cmap.overlap(value, value+0.00001)
            try:
                color = next(iter(color)).data
            except Exception:
                color = None
            bar.set_facecolor(color)
        

def draw_values(
    ax: matplotlib.axes.Axes,
    track: Track,
    *,
    legend_ax: matplotlib.axes.Axes | None = None,
) -> None:
    # Add line
    ax.plot(track.data["chrom_st"], track.data["name"], color="black")
    # Fill in-between
    ax.fill_between(track.data["chrom_st"], track.data["name"], color="black")
    # Trim plot to margins
    ax.margins(x=0, y=0)
    # remove spines
    minimalize_ax(ax, spines=("right", "top"))

    # Limit spine range.
    # TODO: Remove ticks not within bounds.
    ax.spines["bottom"].set_bounds(0, track.data["chrom_end"].max())

    if legend_ax:
        draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center", bbox_to_anchor=[0.5, 0], ncols=3)


def draw_hor_w_ort(
    ax: matplotlib.axes.Axes,
    legend_ax: matplotlib.axes.Axes,
    df_stv: pl.DataFrame,
    df_stv_ort: pl.DataFrame,
):
    minimalize_ax(ax, yticks=True, spines=("right", "left", "top", "bottom"))

    # Add HOR track.
    for row in df_stv.iter_rows(named=True):
        start = row["chrom_st"]
        end = row["chrom_end"]
        color = MONOMER_COLORS.get(row["mer"])
        rect = mpatches.Rectangle(
            (start, 0), end + 1 - start, 12, color=color, lw=0, label=row["mer"]
        )
        ax.add_patch(rect)

    # Add ort track.
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
            -0.1,
            length,
            0,
            width=0.05,
            color="black",
            length_includes_head=True,
            # 10% of entire length
            head_length=0.1 * -length,
            clip_on=False,
        )
        ax.add_patch(arrow)

    draw_uniq_entry_legend(legend_ax, ref_ax=ax, loc="center", bbox_to_anchor=[0.5, 0], ncols=3)


def plot_one_cen(
    df_stv: pl.DataFrame,
    df_stv_ort: pl.DataFrame,
    dfs_track: list[Track],
    chrom: str,
    min_st_pos: int,
    max_end_pos: int,
    width: float,
    height: float,
) -> tuple[Any, matplotlib.axes.Axes]:
    print(f"Plotting {chrom}...", file=sys.stderr)

    track_props = deque([("hor", 0.0)])
    tracks_added = set()
    track_prop_sum = 0.0
    requires_second_col = True
    for track in dfs_track:
        if track.name in tracks_added:
            print(f"Skipping duplicate {track.name} track.", file=sys.stderr)
            continue
        if track.pos == TrackPosition.Overlap:
            pass
        elif track.pos == TrackPosition.Top:
            track_props.appendleft((track.name, track.prop))
            track_prop_sum += track.prop
        elif track.pos == TrackPosition.Bottom:
            track_props.append((track.name, track.prop))
            track_prop_sum += track.prop

        if not requires_second_col and (
            track.opt == TrackOption.Label or track.opt == TrackOption.SelfIdent
        ):
            requires_second_col = True

        tracks_added.add(track.name)

    hor_track_prop = 1.0 - track_prop_sum
    if hor_track_prop < 0.0:
        raise ValueError(
            f"HOR track prop ({hor_track_prop}) is less than 0. Total other track prop: {track_prop_sum}"
        )

    track_info = {track: (i, prop) for i, (track, prop) in enumerate(track_props)}

    track_props = []
    track_indices = {}
    for name, (idx, prop) in sorted(track_info.items(), key=lambda x: x[1][0]):
        if name == "hor":
            track_props.append(hor_track_prop)
        else:
            track_props.append(prop)
        track_indices[name] = idx

    # Adjust columns and width ratio.
    num_cols = 2 if requires_second_col else 1
    width_ratios = (0.8, 0.2) if requires_second_col else 1.0
    track_col = 0
    legend_col = 1

    fig, axes = plt.subplots(
        # Count number of tracks
        len(track_info),
        num_cols,
        figsize=(width, height),
        height_ratios=track_props,
        width_ratios=width_ratios,
    )

    # Get HOR axis.
    hor_track_row = track_indices["hor"]
    ax: matplotlib.axes.Axes = axes[hor_track_row, track_col]
    hor_legend_ax = axes[hor_track_row, legend_col]
    # Set bounds of HOR axis.
    ax.set_xlim(min_st_pos, max_end_pos)

    # HOR
    draw_hor_w_ort(ax=ax, legend_ax=hor_legend_ax, df_stv=df_stv, df_stv_ort=df_stv_ort)

    # Minimalize all legend cols
    for track in dfs_track:
        if track.pos == TrackPosition.Overlap:
            track_ax: matplotlib.axes.Axes = ax
            legend_ax: matplotlib.axes.Axes = axes[hor_track_row, legend_col]
        else:
            track_row = track_indices[track.name]
            track_ax = axes[track_row, track_col]
            legend_ax = axes[track_row, legend_col]

        minimize_legend_opts = {"xticks": True, "yticks": True, "spines": ("right", "left", "top", "bottom")}
        # Switch to line if different track option. {value, label, ident}
        if track.opt == TrackOption.Label:
            draw_label(track_ax, legend_ax, track)

        elif track.opt == TrackOption.SelfIdent:
            draw_self_ident(
                track_ax, track, legend_ax=legend_ax, is_below_hor_track=track_row > hor_track_row
            )
            # Don't remove legend axis elements
            minimize_legend_opts = {"spines": ("right", "top")}

        elif track.opt == TrackOption.Value:
            draw_values(track_ax, track)

        minimalize_ax(
            legend_ax,
            grid=True,
            **minimize_legend_opts,
        )
        # Set label.
        # Encode to ascii and decode to remove escaped characters
        track_ax.set_ylabel(
            track.name.encode("ascii", "ignore").decode("unicode_escape"),
            rotation="horizontal",
            ha="right",
            va='center',
            ma='center',
        )
        # Share axis with main hor axis.
        track_ax.sharex(ax)

    plt.tight_layout()
    plt.savefig("test.png")
    plt.clf()
    breakpoint()

    return fig, ax


def main():
    ap = argparse.ArgumentParser()

    ap.add_argument(
        "-i",
        "--input_stv",
        type=argparse.FileType("rb"),
        help="HOR stv BED9 file",
    )
    ap.add_argument(
        "-t",
        "--input_tracks",
        nargs="*",
        type=str,
        metavar="(name),(overlap|top|bottom),(label|value|self-ident),(prop),(path)",
        help=(
            "Additional headerless BED9 files to plot. Specify track name, track position, track type, proportion of height, and path. "
            "One of more CSV files with the same field can be provided."
        ),
    )
    ap.add_argument("-s", "--sort", help="Sort based on another list.", default=None)
    ap.add_argument("-c", "--chrom", help="Chromosome to plot.", default=None)
    ap.add_argument("-o", "--output", help="Output figure.")
    ap.add_argument(
        "-d",
        "--output_dir",
        help="Output dir to plot multiple separate figures.",
        default=None,
    )
    ap.add_argument(
        "--width",
        type=float,
        help="Figure width in inches per centromere.",
        default=20.0,
    )
    ap.add_argument(
        "--height",
        type=float,
        help="Figure height in inches per centromere.",
        default=12.0,
    )
    ap.add_argument(
        "-m",
        "--mer_order",
        default="large",
        choices=["small", "large"],
        help="HOR monomer order. Small or larger monomers on top.",
    )
    # ap.add_argument("-p", "--processes", type=int, default=4, help="Processes to run.")
    args = ap.parse_args()

    df_stv = read_stv(args.input_stv, chrom=args.chrom)
    dfs_track = get_track_info(args.input_tracks, chrom=args.chrom)
    STV_TRACK = Track("stv", TrackPosition.Center, TrackOption.Label, 0.0, df_stv)

    # Self ident track is already standardized.
    check_tracks = [trk for trk in dfs_track if trk.opt != TrackOption.SelfIdent]
    # Get a reference to the dataframe with the lowest starting position.
    *_, df_track_min_st = min(
        (STV_TRACK, *check_tracks), key=lambda trk: trk.data["chrom_st"].min()
    )

    # Standardize coords
    df_stv, dfs_track = standardize_coords(df_stv, dfs_track, df_track_min_st)
    *_, df_track_min_st = min(
        (STV_TRACK, *check_tracks), key=lambda trk: trk.data["chrom_st"].min()
    )
    *_, df_track_max_end = max(
        (STV_TRACK, *check_tracks), key=lambda trk: trk.data["chrom_end"].max()
    )
    min_st_pos = df_track_min_st["chrom_st"].min()
    max_end_pos = df_track_max_end["chrom_end"].max()

    df_stv_ort = get_stv_mon_ort(df_stv)

    # with ProcessPoolExecutor(max_workers=args.processes) as pool:
    #     filtered_data = [
    #         (
    #             df_stv.filter(pl.col("chrom") == chrom),
    #             df_stv_ort.filter(pl.col("chrom") == chrom),
    #             dfs_track,
    #             chrom,
    #             min_st_pos,
    #             max_end_pos
    #         )
    #         for chrom in df_stv["chrom"].unique()
    #     ]
    #     print("Finished filtering data...", file=sys.stderr)
    #     plots = pool.map(
    #         plot_one_cen,
    #         *zip(*filtered_data),
    #     )

    plots = [
        plot_one_cen(
            df_stv.filter(pl.col("chrom") == chrom),
            df_stv_ort.filter(pl.col("chrom") == chrom),
            [
                Track(
                    trk.name,
                    trk.pos,
                    trk.opt,
                    trk.prop,
                    trk.data.filter(pl.col("chrom") == chrom),
                )
                for trk in dfs_track
            ],
            chrom,
            min_st_pos,
            max_end_pos,
            args.width,
            args.height,
        )
        for chrom in df_stv["chrom"].unique()
    ]

    # TODO: Maybe split per page.
    plt.tight_layout()


if __name__ == "__main__":
    raise SystemExit(main())
