#!/bin/bash

set -euo pipefail

cenplot draw \
-t tracks_example_api.toml \
-p 4 \
-d hor \
-o "hor/merged_hor.png"
