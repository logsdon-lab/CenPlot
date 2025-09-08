import argparse
from matplotlib import rcParams
from .cli.draw import add_draw_cli, draw
from .cli.draw_comparison import add_draw_cmp_cli, draw_cmp

rcParams["pdf.use14corefonts"] = True
rcParams["text.usetex"] = False


def main() -> int:
    ap = argparse.ArgumentParser(description="Centromere ploting library.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_draw_cli(sub_ap)
    add_draw_cmp_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "draw":
        return draw(
            args.input_tracks,
            args.chroms,
            args.outdir,
            args.outfile,
            args.share_xlim,
            args.processes,
        )
    elif args.cmd == "draw_cmp":
        return draw_cmp(
            args.ident,
            args.input_tracks_ref,
            args.input_tracks_qry,
            args.pos_ref,
            args.pos_qry,
            args.ident_colorscale,
            args.ident_prop,
            args.outdir,
        )
    else:
        raise ValueError(f"Not a valid command ({args.cmd})")


if __name__ == "__main__":
    raise SystemExit(main())
