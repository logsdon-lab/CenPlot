import itertools


MONOMER_COLORS = {
    "1": "#A8275C",
    "10": "#9AC78A",
    "11": "#CC8FC1",
    "12": "#3997C6",
    "13": "#8882C4",
    "14": "#8ABDD6",
    "15": "#096858",
    "16": "#45B4CE",
    "17": "#AFA7D8",
    "18": "#A874B5",
    "19": "#3F66A0",
    "2": "#D66C54",
    "20": "#BFDD97",
    "21": "#AF5D87",
    "22": "#E5E57A",
    "24": "#ED975D",
    "25": "#CC99FF",
    "26": "#F9E193",
    "27": "#99CCFF",
    "29": "#004C99",
    "3": "#93430C",
    "30": "#E5D1A1",
    "31": "#FF66B2",
    "32": "#A1B5E5",
    "34": "#9F68A5",
    "35": "#81B25B",
    "36": "#009999",
    "37": "#66FFFF",
    "4": "#F4DC78",
    "42": "#FF9999",
    "44": "#6666FF",
    "46": "#9933FF",
    "5": "#7EC0B3",
    "6": "#B23F73",
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
