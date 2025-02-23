import pytest
import subprocess
import tempfile


# Not perfect. Just need to check nothing crashes.
@pytest.mark.parametrize(
    ["track_file", "ctg_file"],
    [
        # CDR track.
        ("test/tracks_cdr.toml", "test/cdr/HG02953_cdr.bed"),
        # Multiple tracks. ModDotPlot
        ("test/tracks_multiple.toml", "test/chrY/cdrs.bed"),
        ("test/tracks_chr1.toml", "test/chr1/cdrs.bed"),
        # Simple cases.
        ("test/tracks_simple.toml", "test/chrY/cdrs.bed"),
        ("test/tracks_simple_overlap.toml", "test/chrY/cdrs.bed"),
        # Split HOR
        ("test/tracks_split_hor.toml", "test/chr1/cdrs.bed"),
        # YAML
        ("test/tracks_simple.yaml", "test/chrY/cdrs.bed"),
    ],
)
def test_cli_draw(track_file, ctg_file):
    prefix = str(hash(track_file + ctg_file))
    with (
        tempfile.TemporaryDirectory(prefix=prefix) as tmp_dir,
        open(ctg_file, "rt") as fh,
    ):
        # Only read the first line.
        ctg = fh.readlines()[0].split("\t")[0]
        _ = subprocess.run(
            [
                "python",
                "-m",
                "cenplot.main",
                "draw",
                "-t",
                track_file,
                "-c",
                ctg,
                "-d",
                tmp_dir,
                "-p",
                "1",
            ],
            check=True,
        )
