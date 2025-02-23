### Getting Started
Install the package from `pypi`.
```bash
pip install cenplot
```

## CLI
Generating a split HOR tracks using the `cenplot draw` command.
```bash
# examples/example_cli.sh
cenplot draw \
-t tracks_hor.toml \
-c "chm13_chr10:38568472-42561808" \
-p 4 \
-d plots \
-o "plot/merged_image.png"
```

## Python API
The same HOR track can be created with a few lines of code.
```python
# examples/example_api.py
from cenplot import plot_one_cen, read_one_cen_tracks

chrom = "chm13_chr10:38568472-42561808"
track_list, settings = read_one_cen_tracks("tracks_hor.toml", chrom=chrom)
fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
```

### Development
Requires `Git LFS` to pull test files.
Create a `venv`, build `cenplot`, and install it. Also, generate the docs.
```bash
git lfs install && git lfs pull
make dev && make build && make install
pdoc ./cenplot -o docs/
```
