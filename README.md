# `cenplot`
Library for producing centromere figures.

![](docs/example_multiple.png)

> WIP

### Tracks
Tracks can be provided in the format of a `TOML` file under `[[tracks]]`.

#### `name`
The name of a given track. This is added as a label to each track and must be unique.

#### `position`
The given position of a track. Either `relative` or `overlap`.

##### `relative`
Positions the track in relative order within the tracks file.

Here the `CDR` track comes before the `HOR` track.
```toml
[[tracks]]
name = "CDR"
position = "relative"

[[tracks]]
name = "HOR"
position = "relative"
```

##### `overlap`
Positions the track to overlap the previous track within the tracks file.

```toml
[[tracks]]
name = "CDR"
position = "relative"

[[tracks]]
name = "HOR"
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

##### `Value`
Value track plotted as a line and area plot.

##### `SelfIdent`
Self sequence identity plot. See [`ModDotPlot`](https://github.com/marbl/ModDotPlot).


#### `options`
Additional plot options. Dependent on [`type`](#type)

|type|option|description|default|
|-|-|-|-|
|`all`|`legend`|Display the legend.|`false`|
|`all`|`title`|Display the title|`false`|
|`all`|`hide_x`|Hide the x-axis label and ticks|`false`|
|`hor`|`mer_order`|Display this HORs with `x` monomers on top.|`"large"`|
|`horort`|`scale`|Scaling factor for arrow by length.|`50`|
|`horort`|`merge`|Merge same stranded monomers by this number of bases.|`100000`|
|`horort`|`fwd_color`|Color `+` monomers this color.|`"black"`|
|`horort`|`rev_color`||Color `-` monomers this color.|`"black"`|
|`selfident`|`invert`|Invert the self-identity triangle.|`true`|

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
* [ ] Options to struct.
* [ ] Examples
* [ ] Tests
* [ ] Merge images.
