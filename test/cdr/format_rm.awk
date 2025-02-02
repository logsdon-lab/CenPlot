#!/usr/bin/awk -f

# Format repeatmasker bed.
BEGIN {
    OFS="\t"
} {
    strand=($6 == "+") ? "+" : "-";
    if ($4 ~ "ALR" ) {
        color="244,219,112";
        name="asat";
    } else if ($4 ~ "BSR") {
        color="213,177,206";
        name="bsat";
    } else if ($4 ~ "GSAT") {
        color="133,117,170";
        name="gsat";
    } else if ($4 ~ "SAR") {
        color="70,150,139";
        name="hsat1A";
    } else if ($4 ~ "HSATII") {
        color="49,53,100";
        name="hsat2";
    } else if ($4 ~ "HSATI") {
        color="70,150,139";
        name="hsat1B";
    } else if ($4 ~ "(CATTC)n" || $4 ~ "(GAATG)n") {
        color="127,159,184";
        name="hsat3";
    } else {
        color="209,211,212"
        name="ct";
    };
    print $1, $2, $3, name, 0, strand, $2, $3, color
}
