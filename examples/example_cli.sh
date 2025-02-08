#!/bin/bash

set -euo pipefail

cenplot draw \
-t tracks_hor.toml \
-p 4 \
-d plots

cenplot draw \
-t tracks_selfident.toml \
-p 4 \
-d plots

cenplot draw \
-t tracks_bar_label.toml \
-p 4 \
-d plots
