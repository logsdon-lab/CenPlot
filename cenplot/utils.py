from enum import StrEnum, auto
import polars as pl

from collections import defaultdict

from intervaltree import Interval, IntervalTree


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
