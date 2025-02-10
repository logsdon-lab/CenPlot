r"""
[![PyPI - Version](https://img.shields.io/pypi/v/cenplot)](https://pypi.org/project/cenplot/)
[![CI](https://github.com/logsdon-lab/cenplot/actions/workflows/main.yaml/badge.svg)](https://github.com/logsdon-lab/cenplot/actions/workflows/main.yaml)
[![docs](https://github.com/logsdon-lab/cenplot/actions/workflows/docs.yaml/badge.svg)](https://github.com/logsdon-lab/cenplot/actions/workflows/docs.yaml)

A library for building centromere figures.

<figure float="left">
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/example_multiple.png" width="100%">
</figure>

# Quickstart

.. include:: ../docs/quickstart.md


# Overview
Configuration comes in the form of `TOML` files with two fields, `[settings]` and `[[tracks]]`.
```toml
[settings]
format = "png"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "relative"

[[tracks]]
title = "Sequence Composition"
position = "relative"
```

`[[settings]]` determines figure level settings while `[[tracks]]` determines track level settings.
* To view all of the possible options for `[[settings]]`, see `cenplot.PlotSettings`
* To view all of the possible options for `[[tracks]]`, see one of `cenplot.TrackSettings`

## Track Order
Order is determined by placement of tracks. Here the `"Alpha-satellite HOR monomers"` comes before the `"Sequence Composition"` track.
```toml
[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "relative"

[[tracks]]
title = "Sequence Composition"
position = "relative"
```

<figure float="left">
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/simple_hor_top.png" width="100%">
</figure>

Reversing this will plot `"Sequence Composition"` before `"Alpha-satellite HOR monomers"`

```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "relative"
```

<figure float="left">
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/simple_hor_bottom.png" width="100%">
</figure>

## Overlap
Tracks can be overlapped with the `position` or `cenplot.TrackPosition` setting.

```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "overlap"
```

<figure float="left">
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/simple_hor_overlap.png" width="100%">
</figure>

The preceding track is overlapped and the legend elements are merged.

## Track Types and Data
Track types, or `cenplot.TrackType`s, are specified via the `type` parameter.
```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
type = "label"
path = "rm.bed"
```

Each type will expect different BED files in the `path` option.
* For example, the option `TrackType.SelfIdent` expects the following values.

|query|query_st|query_end|reference|reference_st|reference_end|percent_identity_by_events|
|-|-|-|-|-|-|-|
|x|1|5000|x|1|5000|100.0|

When using the `Python` API, each will have an associated `read_*` function (ex. `cenplot.read_bed_identity`).
* Using `cenplot.read_one_cen_tracks` is preferred.

> [!NOTE] If input BED files have contigs with coordinates in their name, the coordinates are expected to be in absolute coordinates.

Absolute coordinates
|chrom|chrom_st|chrom_end|
|-|-|-|
|chm13:100-200|105|130|

## Proportion
Each track must account for some proportion of the total plot dimensions.
* The plot dimensions are specified with `cenplot.PlotSettings.dim`

Here, with a total proportion of `0.2`, each track will take up `50%` of the total plot dimensions.
```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
type = "label"
proportion = 0.1
path = "rm.bed"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "relative"
type = "hor"
proportion = 0.1
path = "stv_row.bed"
```

When the position is `cenplot.TrackPosition.Overlap`, the proportion is ignored.
```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
type = "label"
proportion = 0.1
path = "rm.bed"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "overlap"
type = "hor"
path = "stv_row.bed"
```

## Options
Options for specific `cenplot.TrackType` types can be specified in `options`.
* See `cenplot.TrackSettings`

```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
proportion = 0.5
type = "label"
path = "rm.bed"
# Both need to be false to keep x
options = { hide_x = false }

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "overlap"
type = "hor"
path = "stv_row.bed"
# Change mode to showing HOR variant and reduce legend number of cols.
options = { hide_x = false, mode = "hor", legend_ncols = 2 }
```

<figure float="left">
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/simple_hor_track_options.png" width="100%">
</figure>

## Examples
Examples of both the CLI and Python API can be found in the root of `cenplot`'s project directory under `examples/` or `test/`

---
"""

import logging

from .lib.draw import (
    draw_hor,
    draw_hor_ort,
    draw_label,
    draw_self_ident,
    draw_bars,
    draw_legend,
    draw_self_ident_hist,
    plot_one_cen,
    merge_plots,
    PlotSettings,
)
from .lib.io import (
    read_bed9,
    read_bed_hor,
    read_bed_identity,
    read_bed_label,
    read_one_cen_tracks,
)
from .lib.track import (
    Track,
    TrackType,
    TrackPosition,
    TrackList,
    LegendPosition,
    TrackSettings,
    SelfIdentTrackSettings,
    HORTrackSettings,
    HOROrtTrackSettings,
    BarTrackSettings,
    LabelTrackSettings,
    PositionTrackSettings,
    LegendTrackSettings,
    SpacerTrackSettings,
)

__author__ = "Keith Oshima (oshimak@pennmedicine.upenn.edu)"
__license__ = "MIT"
__all__ = [
    "plot_one_cen",
    "merge_plots",
    "draw_hor",
    "draw_hor_ort",
    "draw_label",
    "draw_self_ident",
    "draw_self_ident_hist",
    "draw_bars",
    "draw_legend",
    "read_bed9",
    "read_bed_hor",
    "read_bed_identity",
    "read_bed_label",
    "read_one_cen_tracks",
    "Track",
    "TrackType",
    "TrackPosition",
    "TrackList",
    "LegendPosition",
    "PlotSettings",
    "TrackSettings",
    "SelfIdentTrackSettings",
    "HORTrackSettings",
    "HOROrtTrackSettings",
    "BarTrackSettings",
    "LabelTrackSettings",
    "PositionTrackSettings",
    "LegendTrackSettings",
    "SpacerTrackSettings",
]

logging.getLogger(__name__).addHandler(logging.NullHandler())
