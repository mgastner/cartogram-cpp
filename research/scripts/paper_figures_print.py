#!/usr/bin/env python3
"""Draw the paper's figures at print size, from the outputs of paper_runs.py.

Usage: paper_figures_print.py OUT_DIR [--text-width-pt 482] [--python PYTHON]

The figures of paper_runs.py --figures are drawn for the screen: scaled to
the text width of the ICA template (17 cm, about 482 pt), their labels would
print at 3-4 pt. This script draws the same figures with smaller panels, so
that the labels print at about 7 pt or more when a figure spans the text
width, and writes them into OUT_DIR/figs_print. For every figure, it prints
the width of the PDF and the smallest font size after scaling to the text
width.

Figures: the graticule and the distortion field of every run with 32 modes
along the longer side, and the convergence of the USA and China runs.
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).parent

# Panel widths in pixels, chosen so that the smallest label prints at about
# 7 pt at a text width of 482 pt
GRATICULE_PANEL = 480
DISTORTION_PANEL = 200
CONVERGENCE_PANEL = 260


def smallest_printed_font(svg, text_width_pt):
    s = svg.read_text()
    m = re.search(r'<svg[^>]*\swidth="([0-9.]+)"', s)
    width_px = float(m.group(1))
    sizes = [float(v) for v in re.findall(r"font-size=[\"']([0-9.]+)[\"']", s)]
    # rsvg-convert maps 1 px to 0.75 pt; the figure is then scaled to the text width
    return width_px * 0.75, min(sizes) * text_width_pt / width_px


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--text-width-pt", type=float, default=482.0)
    parser.add_argument("--python", default=sys.executable)
    args = parser.parse_args()

    figs = args.out_dir / "figs_print"
    figs.mkdir(exist_ok=True)
    written = []
    for run in sorted(args.out_dir.iterdir()):
        if not run.is_dir() or not re.fullmatch(r"[a-z]+32", run.name):
            continue
        plot_json = next(run.glob("*_conformality_plot.json"))
        for script, suffix, width in [
            ("plot_conformality.py", "graticule", GRATICULE_PANEL),
            ("plot_distortion_field.py", "distortion", DISTORTION_PANEL),
        ]:
            out = figs / f"{run.name}_{suffix}.svg"
            subprocess.run(
                [args.python, str(SCRIPTS / script), str(plot_json), "-o", str(out), "--width", str(width)],
                check=True, stdout=subprocess.DEVNULL,
            )
            written.append(out)
    conv = [args.out_dir / n for n in ("usa16", "usa32", "usa64", "china16", "china32", "china64")
            if (args.out_dir / n).is_dir()]
    if conv:
        out = figs / "convergence_usa_china.svg"
        subprocess.run(
            [args.python, str(SCRIPTS / "plot_convergence.py"), *map(str, conv), "-o", str(out),
             "--panel-width", str(CONVERGENCE_PANEL)],
            check=True, stdout=subprocess.DEVNULL,
        )
        written.append(out)

    print(f"| Figure | PDF width (pt) | Smallest label at {args.text_width_pt:g} pt text width |")
    print("|---|---|---|")
    for svg in written:
        width_pt, font_pt = smallest_printed_font(svg, args.text_width_pt)
        print(f"| {svg.stem} | {width_pt:.0f} | {font_pt:.1f} pt |")
    print(f"\nwrote the figures into {figs}")


if __name__ == "__main__":
    main()
