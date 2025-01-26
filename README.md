# `cenplot`
Library for producing centromere figures.

> WIP

### Test

```toml
[tracks]
[[tracks.plots]]
name = "CDR"
position = "relative"
type = "label"
proportion = 0.025
path = "test/chrY_cdrs.bed"

[[tracks.plots]]
name = "Approximate\nNucleotide\nIdentity"
position = "relative"
type = "selfident"
proportion = 0.7
path = "test/chrY_ident.bed"
options = { flip_y = true }
```

```bash
python scripts/plot_cens_all_stv_mpl.py -t test/tracks.toml
```

### TODO:
* [] Monomer order
* [] Split monomers
* [] Refactor to library.
* [] Examples
* [] Tests
* [] Merge images.
