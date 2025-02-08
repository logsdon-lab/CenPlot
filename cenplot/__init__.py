r"""
A library for building centromere figures via tracks.

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
* To view all of the possible options for `[[settings]]`, see `cenplot.SinglePlotSettings`
* To view all of the possible options for `[[tracks]]`, see one of `cenplot.PlotSettings`

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
    <img align="middle" src="https://raw.githubusercontent.com/logsdon-lab/cenplot/refs/heads/main/docs/simple_overlap.png" width="100%">
</figure>

The preceding track is overlapped and the legend elements are merged.

## Proportion
Each track must account for some proportion of the total plot dimensions.
* The plot dimensions are specified with `cenplot.SinglePlotSettings.dim`

Here, with a total proportion of `0.2`, each track will take up `50%` of the total plot dimensions.
```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
proportion = 0.1
path = "rm.bed"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "relative"
proportion = 0.1
path = "stv_row.bed"
```

When the position is `cenplot.TrackPosition.Overlap`, the proportion is ignored.
```toml
[[tracks]]
title = "Sequence Composition"
position = "relative"
proportion = 0.1
path = "rm.bed"

[[tracks]]
title = "Alpha-satellite HOR monomers"
position = "overlap"
path = "stv_row.bed"
```
---
"""

import logging

from .lib.draw import (
    draw_hor,
    draw_hor_ort,
    draw_label,
    draw_self_ident,
    draw_bars,
    plot_one_cen,
    merge_plots,
    SinglePlotSettings,
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
    TrackOption,
    TrackPosition,
    TrackList,
    LegendPosition,
    PlotSettings,
    SelfIdentPlotSettings,
    HORPlotSettings,
    HOROrtPlotSettings,
    BarPlotSettings,
    LabelPlotSettings,
    PositionPlotSettings,
    LegendPlotSettings,
    SpacerPlotSettings,
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
    "draw_bars",
    "read_bed9",
    "read_bed_hor",
    "read_bed_identity",
    "read_bed_label",
    "read_one_cen_tracks",
    "Track",
    "TrackOption",
    "TrackPosition",
    "TrackList",
    "LegendPosition",
    "PlotSettings",
    "SelfIdentPlotSettings",
    "HORPlotSettings",
    "HOROrtPlotSettings",
    "BarPlotSettings",
    "LabelPlotSettings",
    "SinglePlotSettings",
    "PositionPlotSettings",
    "LegendPlotSettings",
    "SpacerPlotSettings",
]

logging.getLogger(__name__).addHandler(logging.NullHandler())
