import os
from cenplot import plot_tracks, read_tracks


def hor():
    chrom = "chm13_chr10:38568472-42561808"
    tracks = os.path.abspath("tracks_hor.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes


def self_ident():
    chrom = "HG00731_chrY_haplotype2-0000041:9700692-11101963"
    tracks = os.path.abspath("tracks_selfident.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes


def bar_label():
    chrom = "haplotype1-0000003"
    tracks = os.path.abspath("tracks_bar_label.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes


def local_self_ident():
    chrom = "HG00096_chr1_haplotype1-0000018"
    tracks = os.path.abspath("tracks_local_selfident.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes


def local_self_ident_subset():
    chrom = "HG00096_chr1_haplotype1-0000018:1000000-3000000"
    tracks = os.path.abspath("tracks_local_selfident.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes


def strand():
    chrom = "chm13_chr1:121119216-127324115"
    tracks = os.path.abspath("tracks_strand.toml")
    with open(tracks, "rb") as fh:
        track_list, settings = read_tracks(fh, chrom=chrom)
    fig, axes, _ = plot_tracks(track_list.tracks, settings)
    return fig, axes
