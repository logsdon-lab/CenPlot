import itertools


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

IDENT_RANGE_ENDS = tuple(
    # Increment by IDENT_INCREMENT from IDENT_CUTOFF up to 100%
    i + IDENT_CUTOFF
    for i in itertools.accumulate(IDENT_INCREMENT for _ in range(10))
)
IDENT_RANGE = [
    (0, 90),
    (90, 97.5),
    *zip((IDENT_CUTOFF, *IDENT_RANGE_ENDS[:-1]), IDENT_RANGE_ENDS),
]

IDENT_COLOR_RANGE = dict(zip(IDENT_RANGE, IDENT_COLORS))
