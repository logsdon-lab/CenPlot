import argparse
from matplotlib import rcParams
from .cli.draw import add_draw_cli, draw

rcParams["pdf.use14corefonts"] = True
rcParams["text.usetex"] = False


def main() -> int:
    ap = argparse.ArgumentParser(description="Centromere ploting library.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_draw_cli(sub_ap)

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
    else:
        raise ValueError(f"Not a valid command ({args.cmd})")


if __name__ == "__main__":
    raise SystemExit(main())
