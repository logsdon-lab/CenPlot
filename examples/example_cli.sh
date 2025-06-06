#!/bin/bash

set -euo pipefail

cenplot draw \
-t tracks_hor.toml \
-c "chm13_chr10:38568472-42561808" \
-p 4 \
-d plots

cenplot draw \
-t tracks_selfident.toml \
-c "HG00731_chrY_haplotype2-0000041:9700692-11101963" \
-p 4 \
-d plots

cenplot draw \
-t tracks_bar_label.toml \
-c "haplotype1-0000003" \
-p 4 \
-d plots

cenplot draw \
-t tracks_local_selfident.toml \
-c "HG00096_chr1_haplotype1-0000018" \
-p 4 \
-d plots
