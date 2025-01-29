# `cenplot`
Library for producing centromere figures.

![](docs/example_multiple.png)

> WIP

### Plot Settings
General settings for output plots.

#### `single`
```toml
[settings]
format = "png"
transparent = true
dim = [16.0, 12.0]
dpi = 600
legend_pos = "right"
axis_h_pad = 0.2
```

#### `multiple`
```toml
```

### Tracks
Tracks can be provided in the format of a `TOML` file under `[[tracks]]`.

#### `title`
The title of a given track. This is added as a label to each track.
* The chrom can be added by using format string formatting. ex. `{chrom}`

#### `position`
The given position of a track. Either `relative` or `overlap`.

##### `relative`
Positions the track in relative order within the tracks file.

Here the `CDR` track comes before the `HOR` track.
```toml
[[tracks]]
title = "CDR"
position = "relative"

[[tracks]]
title = "HOR"
position = "relative"
```

##### `overlap`
Positions the track to overlap the previous track within the tracks file.

```toml
[[tracks]]
title = "CDR"
position = "relative"

[[tracks]]
title = "HOR"
position = "overlap"
```

#### `type`
The track type.

##### `HOR`
Higher Order Repeat track. Orientation of each HOR is added as an arrow.

##### `HORSplit`
Higher Order Repeat track. Split by number of monomers in a HOR.

##### `Label`
Label track. Each label is plotted as a bar on a single track.

##### `Bar`
Bar plot for each interval in the BED file.

##### `SelfIdent`
Self sequence identity plot. See [`ModDotPlot`](https://github.com/marbl/ModDotPlot).


#### `options`
Additional plot options. Dependent on [`type`](#type)

|type|option|description|default|
|-|-|-|-|
|`all`|`legend`|Display the legend.|`False`|
|`all`|`legend_title`|Title for legend.|`None`|
|`all`|`hide_x`|Hide the x-axis label and ticks|`False`|
|`all`|`fontsize`|Fontsize for axis elements.|`"medium"`|
|`all`|`title_fontsize`|Fontsize for axis title.|`"medium"`|
|`all`|`legend_fontsize`|Fontsize for legend elements|`"x-large"`|
|`all`|`legend_title_fontsize`|Fontsize for legend title.|`"x-large"`|
|`hor`|`mer_order`|Display this HORs with `x` monomers on top.|`"large"`|
|`horort`|`scale`|Scaling factor for arrow by length.|`50`|
|`horort`|`merge`|Merge same stranded monomers by this number of bases.|`100000`|
|`horort`|`fwd_color`|Color `+` monomers this color.|`"black"`|
|`horort`|`rev_color`|Color `-` monomers this color.|`"black"`|
|`label`|`color`|Label color. Used if no color is provided in `item_rgb` column.|`None`|
|`label`|`alpha`|Label alpha.|`1.0`|
|`bar`|`color`|Bar color. Used if no color is provided in `item_rgb` column.|`None`|
|`bar`|`alpha`|Bar alpha.|`1.0`|
|`selfident`|`invert`|Invert the self-identity triangle.|`True`|
|`selfident`|`legend_bins`| Number of bins for `perc_identity_by_events` in the legend.|`300`|
|`selfident`|`legend_xmin`| Legend x-min coordinate. Used to constrain x-axis limits.|`70.0`|
|`selfident`|`legend_asp_ratio`|  Aspect ratio of legend. If `None`, takes up entire axis.|`1.0`|

#### Example:
```toml
[[tracks]]
name = "CDR"
position = "relative"
type = "label"
proportion = 0.025
path = "test/chrY/cdrs.bed"

[[tracks]]
name = "Approximate\nNucleotide\nIdentity"
position = "relative"
type = "selfident"
proportion = 0.7
path = "test/chrY/ident.bed"
options = { invert = true }
```

```bash
make venv && make build && make install
python scripts/plot.py -t test/tracks_multiple.toml -ht 12
```

### TODO:
* [ ] Monomer order
* [x] Split monomers
* [x] Refactor to library.
* [x] Options to struct.
* [ ] Examples
* [ ] Tests
* [ ] Merge images.
