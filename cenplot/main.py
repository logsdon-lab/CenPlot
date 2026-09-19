import argparse
import matplotlib.pyplot as plt

from matplotlib import rcParams
from .cli.draw import add_draw_cli, draw

rcParams["pdf.use14corefonts"] = True
rcParams["text.usetex"] = False

plt.switch_backend("agg")


def main() -> int:
    ap = argparse.ArgumentParser(description="Centromere ploting library.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_draw_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "draw":
        draw(
            args.input_tracks,
            args.chroms,
            args.outdir,
            args.outfile,
            args.share_xlim,
            args.processes,
        )
        return 0
    else:
        raise ValueError(f"Not a valid command ({args.cmd})")


if __name__ == "__main__":
    raise SystemExit(main())
