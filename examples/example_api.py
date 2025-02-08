from cenplot import plot_one_cen, read_one_cen_tracks


def hor():
    chrom = "chm13_chr10:38568472-42561808"
    tracks = "tracks_example_api.toml"
    track_list, settings = read_one_cen_tracks(tracks, chrom=chrom)
    fig, axes, outfile = plot_one_cen(track_list.tracks, "plots", chrom, settings)
    return fig, axes, outfile
