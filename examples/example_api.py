from cenplot import plot_one_cen, read_one_cen_tracks


def hor():
    chrom = "chm13_chr10:38568472-42561808"
    tracks = "tracks_selfident.toml"
    with open(tracks, "rb") as fh:
        track_list, settings = read_one_cen_tracks(fh, chrom=chrom)
    fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
    return fig, axes, outfile


def self_ident():
    chrom = "HG00731_chrY_haplotype2-0000041:9700692-11101963"
    tracks = "tracks_selfident.toml"
    with open(tracks, "rb") as fh:
        track_list, settings = read_one_cen_tracks(fh, chrom=chrom)
    fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
    return fig, axes, outfile


def bar_label():
    chrom = "haplotype1-0000003"
    tracks = "tracks_bar_label.toml"
    with open(tracks, "rb") as fh:
        track_list, settings = read_one_cen_tracks(fh, chrom=chrom)
    fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
    return fig, axes, outfile


def strand():
    chrom = "chm13_chr1:121119216-127324115"
    tracks = "tracks_strand.toml"
    with open(tracks, "rb") as fh:
        track_list, settings = read_one_cen_tracks(fh, chrom=chrom)
    fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
    return fig, axes, outfile
