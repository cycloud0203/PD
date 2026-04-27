#!/usr/bin/env python3
"""
Floorplan result visualizer.

Reads a block description file (.block) together with its placement report
(.rpt) produced by the floorplanner and draws the resulting layout using
matplotlib.

Report (.rpt) format expected:
  line 1 : total cost (alpha * A + (1 - alpha) * W)
  line 2 : total wirelength
  line 3 : total chip area
  line 4 : chip_width chip_height
  line 5 : runtime (seconds)
  line 6+: <block_name> <x1> <y1> <x2> <y2>

Block (.block) format expected:
  Outline: <W> <H>
  NumBlocks: <N>
  NumTerminals: <T>
  <name> <w> <h>                   (N lines, blocks)
  <name> terminal <x> <y>          (T lines, fixed IO pads)

Usage:
  Single report:
    python scripts/visualize.py <report.rpt> [--block <file.block>]
                                [--nets <file.nets>] [--out <image>]
                                [--no-show] [--no-labels] [--no-terminals]

  Batch (visualize every .rpt under a directory, PNG per report):
    python scripts/visualize.py --all [--input-dir output]
                                [--out-dir output/viz]
                                [--nets auto|none] [--no-labels] [--no-terminals]

If --block is omitted the script tries to locate the matching .block file in
./input_pa2 based on the report stem (e.g. output/apte.rpt     -> input_pa2/apte.block,
output/ami33_0.5.rpt -> input_pa2/ami33.block).
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------
@dataclass
class PlacedBlock:
    name: str
    x1: float
    y1: float
    x2: float
    y2: float

    @property
    def w(self) -> float:
        return self.x2 - self.x1

    @property
    def h(self) -> float:
        return self.y2 - self.y1

    @property
    def cx(self) -> float:
        return 0.5 * (self.x1 + self.x2)

    @property
    def cy(self) -> float:
        return 0.5 * (self.y1 + self.y2)


@dataclass
class Report:
    cost: float
    wirelength: float
    area: float
    chip_w: float
    chip_h: float
    runtime: float
    blocks: List[PlacedBlock] = field(default_factory=list)


@dataclass
class BlockInfo:
    outline_w: float
    outline_h: float
    block_sizes: Dict[str, Tuple[float, float]] = field(default_factory=dict)
    terminals: Dict[str, Tuple[float, float]] = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------
def parse_report(path: str) -> Report:
    with open(path, "r", encoding="utf-8") as f:
        raw = [ln.strip() for ln in f if ln.strip() != ""]

    if len(raw) < 5:
        raise ValueError(
            "Report '%s' has only %d non-empty lines; expected >= 5." % (path, len(raw))
        )

    cost = float(raw[0])
    wirelength = float(raw[1])
    area = float(raw[2])
    w_h = raw[3].split()
    if len(w_h) < 2:
        raise ValueError("Line 4 of '%s' must contain chip width and height." % path)
    chip_w = float(w_h[0])
    chip_h = float(w_h[1])
    runtime = float(raw[4].split()[0])

    blocks: List[PlacedBlock] = []
    for ln in raw[5:]:
        parts = ln.split()
        if len(parts) < 5:
            continue
        name = parts[0]
        x1, y1, x2, y2 = (float(v) for v in parts[1:5])
        if x2 < x1:
            x1, x2 = x2, x1
        if y2 < y1:
            y1, y2 = y2, y1
        blocks.append(PlacedBlock(name, x1, y1, x2, y2))

    return Report(cost, wirelength, area, chip_w, chip_h, runtime, blocks)


def parse_block(path: str) -> BlockInfo:
    info = BlockInfo(0.0, 0.0)
    with open(path, "r", encoding="utf-8") as f:
        lines = [ln.rstrip("\n") for ln in f]

    n_blocks = 0
    n_terms = 0
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        low = s.lower()
        if low.startswith("outline"):
            parts = s.replace(":", " ").split()
            info.outline_w = float(parts[1])
            info.outline_h = float(parts[2])
        elif low.startswith("numblocks"):
            n_blocks = int(s.split(":")[1])
        elif low.startswith("numterminals"):
            n_terms = int(s.split(":")[1])
            break

    remaining = [ln.strip() for ln in lines if ln.strip() != ""]
    # Drop the first 3 header entries (Outline, NumBlocks, NumTerminals).
    header_seen = 0
    body: List[str] = []
    for ln in remaining:
        low = ln.lower()
        if header_seen < 3 and (
            low.startswith("outline")
            or low.startswith("numblocks")
            or low.startswith("numterminals")
        ):
            header_seen += 1
            continue
        body.append(ln)

    # First n_blocks body lines are soft blocks, next n_terms are terminals.
    for ln in body[:n_blocks]:
        parts = ln.split()
        if len(parts) >= 3:
            info.block_sizes[parts[0]] = (float(parts[1]), float(parts[2]))

    for ln in body[n_blocks : n_blocks + n_terms]:
        parts = ln.split()
        if len(parts) >= 4 and parts[1].lower() == "terminal":
            info.terminals[parts[0]] = (float(parts[2]), float(parts[3]))
        elif len(parts) >= 3:
            # Fallback: "<name> <x> <y>".
            info.terminals[parts[0]] = (float(parts[-2]), float(parts[-1]))

    return info


def parse_nets(path: str) -> List[List[str]]:
    """Parse a .nets file into a list of nets (each net is a list of pin names)."""
    nets: List[List[str]] = []
    current: List[str] = []
    expected = 0
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            s = ln.strip()
            if not s:
                continue
            low = s.lower()
            if low.startswith("numnets"):
                continue
            if low.startswith("netdegree"):
                if current:
                    nets.append(current)
                current = []
                try:
                    expected = int(s.split(":")[1])
                except (IndexError, ValueError):
                    expected = 0
                continue
            current.append(s)
            if expected and len(current) == expected:
                nets.append(current)
                current = []
                expected = 0
    if current:
        nets.append(current)
    return nets


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
def _color_for(idx: int) -> Tuple[float, float, float]:
    # Matplotlib 'tab20' provides 20 distinguishable colors; cycle after that.
    cmap = plt.get_cmap("tab20")
    return cmap(idx % 20)


def plot_floorplan(
    report: Report,
    info: BlockInfo,
    nets: List[List[str]],
    out_path: str,
    show: bool,
    show_labels: bool,
    show_terminals: bool,
    title_suffix: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 10))

    blk_min_x = min((b.x1 for b in report.blocks), default=0.0)
    blk_min_y = min((b.y1 for b in report.blocks), default=0.0)
    blk_max_x = max((b.x2 for b in report.blocks), default=0.0)
    blk_max_y = max((b.y2 for b in report.blocks), default=0.0)
    bbox_w = max(0.0, blk_max_x - blk_min_x)
    bbox_h = max(0.0, blk_max_y - blk_min_y)

    # Extent of the drawing area.
    max_x = max(info.outline_w, blk_max_x)
    max_y = max(info.outline_h, blk_max_y)
    min_x = min(0.0, blk_min_x)
    min_y = min(0.0, blk_min_y)
    if show_terminals and info.terminals:
        tx = [p[0] for p in info.terminals.values()]
        ty = [p[1] for p in info.terminals.values()]
        min_x = min(min_x, min(tx))
        min_y = min(min_y, min(ty))
        max_x = max(max_x, max(tx))
        max_y = max(max_y, max(ty))

    # Fixed outline constraint from the .block file.
    if info.outline_w > 0 and info.outline_h > 0:
        ax.add_patch(
            mpatches.Rectangle(
                (0, 0),
                info.outline_w,
                info.outline_h,
                linewidth=2.0,
                edgecolor="red",
                facecolor="none",
                linestyle="--",
                label="Outline (%g x %g)" % (info.outline_w, info.outline_h),
            )
        )

    # Actual bounding box of the placement.
    ax.add_patch(
        mpatches.Rectangle(
            (blk_min_x, blk_min_y),
            bbox_w,
            bbox_h,
            linewidth=1.2,
            edgecolor="black",
            facecolor="none",
            linestyle="-",
            label="Chip bbox (%g x %g)" % (bbox_w, bbox_h),
        )
    )

    # Blocks.
    for idx, blk in enumerate(report.blocks):
        color = _color_for(idx)
        ax.add_patch(
            mpatches.Rectangle(
                (blk.x1, blk.y1),
                blk.w,
                blk.h,
                linewidth=0.8,
                edgecolor="black",
                facecolor=color,
                alpha=0.55,
            )
        )
        if show_labels:
            ax.text(
                blk.cx,
                blk.cy,
                blk.name,
                ha="center",
                va="center",
                fontsize=max(6, min(11, 0.6 * min(blk.w, blk.h) / max(1.0, max_x) * 60)),
                color="black",
            )

    # Terminals.
    if show_terminals and info.terminals:
        tx = [p[0] for p in info.terminals.values()]
        ty = [p[1] for p in info.terminals.values()]
        ax.scatter(tx, ty, s=10, c="black", marker="o", label="Terminals", zorder=3)

    # Nets (HPWL polylines through centers and terminals).
    if nets:
        pin_pos: Dict[str, Tuple[float, float]] = {
            b.name: (b.cx, b.cy) for b in report.blocks
        }
        pin_pos.update(info.terminals)
        segments: List[List[Tuple[float, float]]] = []
        for net in nets:
            pts = [pin_pos[p] for p in net if p in pin_pos]
            if len(pts) < 2:
                continue
            # Draw the HPWL bounding box edges for visual clarity.
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            bx1, bx2 = min(xs), max(xs)
            by1, by2 = min(ys), max(ys)
            segments.append([(bx1, by1), (bx2, by1)])
            segments.append([(bx2, by1), (bx2, by2)])
            segments.append([(bx2, by2), (bx1, by2)])
            segments.append([(bx1, by2), (bx1, by1)])
        if segments:
            lc = LineCollection(
                segments, colors="steelblue", linewidths=0.3, alpha=0.25, zorder=2
            )
            ax.add_collection(lc)

    # Padding for display.
    pad_x = 0.03 * max(max_x - min_x, 1.0)
    pad_y = 0.03 * max(max_y - min_y, 1.0)
    ax.set_xlim(min_x - pad_x, max_x + pad_x)
    ax.set_ylim(min_y - pad_y, max_y + pad_y)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    fit = "N/A"
    if info.outline_w > 0 and info.outline_h > 0:
        fit = (
            "fit"
            if blk_min_x >= 0.0 and blk_min_y >= 0.0 and blk_max_x <= info.outline_w and blk_max_y <= info.outline_h
            else "OVERFLOW"
        )
    title = (
        "Floorplan: %s\n"
        "cost=%.4f  WL=%.0f  A=%.0f  bbox=%gx%g  time=%.3fs  [%s]"
        % (
            title_suffix,
            report.cost,
            report.wirelength,
            report.area,
            report.chip_w,
            report.chip_h,
            report.runtime,
            fit,
        )
    )
    ax.set_title(title, fontsize=11)
    ax.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax.grid(True, linestyle=":", linewidth=0.3, alpha=0.5)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print("Saved visualization to %s" % out_path)
    if show:
        plt.show()
    plt.close(fig)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _guess_block(rpt_path: str) -> str:
    stem = os.path.splitext(os.path.basename(rpt_path))[0]
    # Strip trailing suffixes like "_0.5", "_test", etc., by taking the first
    # underscore-separated token.
    base = stem.split("_")[0]
    candidate = os.path.join("input_pa2", base + ".block")
    return candidate


def _guess_nets(block_path: str) -> str:
    if block_path.endswith(".block"):
        return block_path[: -len(".block")] + ".nets"
    return block_path + ".nets"


def _render_one(
    rpt_path: str,
    block_override: str,
    nets_mode: str,
    out_path: str,
    show: bool,
    show_labels: bool,
    show_terminals: bool,
) -> bool:
    """Render a single report to an image. Returns True on success."""
    if not os.path.isfile(rpt_path):
        print("ERROR: report not found: %s" % rpt_path, file=sys.stderr)
        return False

    block_path = block_override or _guess_block(rpt_path)
    if not os.path.isfile(block_path):
        print(
            "ERROR: block file not found for %s (tried %s)" % (rpt_path, block_path),
            file=sys.stderr,
        )
        return False

    if nets_mode == "auto":
        nets_path = _guess_nets(block_path)
    elif nets_mode == "none" or not nets_mode:
        nets_path = ""
    else:
        nets_path = nets_mode

    nets: List[List[str]] = []
    if nets_path:
        if not os.path.isfile(nets_path):
            print(
                "WARNING: nets file not found, skipping: %s" % nets_path,
                file=sys.stderr,
            )
        else:
            nets = parse_nets(nets_path)

    try:
        report = parse_report(rpt_path)
    except Exception as exc:
        print("ERROR: failed to parse %s: %s" % (rpt_path, exc), file=sys.stderr)
        return False
    info = parse_block(block_path)

    plot_floorplan(
        report=report,
        info=info,
        nets=nets,
        out_path=out_path,
        show=show,
        show_labels=show_labels,
        show_terminals=show_terminals,
        title_suffix=os.path.basename(rpt_path),
    )
    return True


def main(argv: List[str]) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "report",
        nargs="?",
        help="Placement report file (.rpt). Omit when --all is used.",
    )
    ap.add_argument("--block", help="Matching .block file (auto-detected if omitted)")
    ap.add_argument(
        "--nets",
        help="Matching .nets file. Use 'auto' to try the sibling of --block, "
        "'none' to disable (default: none)",
        default="none",
    )
    ap.add_argument("--out", help="Output image file (default: <report>.png)")
    ap.add_argument(
        "--all",
        action="store_true",
        help="Batch mode: visualize every .rpt under --input-dir into --out-dir.",
    )
    ap.add_argument(
        "--input-dir",
        default="output",
        help="Directory to scan for .rpt files in --all mode (default: output)",
    )
    ap.add_argument(
        "--out-dir",
        default="output/viz",
        help="Directory to write PNGs into in --all mode (default: output/viz)",
    )
    ap.add_argument("--no-show", action="store_true", help="Do not open a window")
    ap.add_argument("--no-labels", action="store_true", help="Omit block name labels")
    ap.add_argument("--no-terminals", action="store_true", help="Hide terminal pins")
    args = ap.parse_args(argv)

    # ---------------- Batch mode ----------------
    if args.all:
        if args.report is not None:
            print(
                "ERROR: do not pass a positional report when using --all.",
                file=sys.stderr,
            )
            return 2
        if not os.path.isdir(args.input_dir):
            print(
                "ERROR: input directory not found: %s" % args.input_dir,
                file=sys.stderr,
            )
            return 2

        rpt_files = sorted(glob.glob(os.path.join(args.input_dir, "*.rpt")))
        if not rpt_files:
            print(
                "ERROR: no .rpt files found under %s" % args.input_dir,
                file=sys.stderr,
            )
            return 1

        os.makedirs(args.out_dir, exist_ok=True)
        print(
            "Visualizing %d report(s) from %s -> %s"
            % (len(rpt_files), args.input_dir, args.out_dir)
        )

        # In batch mode we never pop a window (ignore any lingering --show intent).
        show = False
        n_ok = 0
        n_fail = 0
        for rpt in rpt_files:
            stem = os.path.splitext(os.path.basename(rpt))[0]
            out_path = os.path.join(args.out_dir, stem + ".png")
            ok = _render_one(
                rpt_path=rpt,
                block_override=args.block or "",
                nets_mode=args.nets,
                out_path=out_path,
                show=show,
                show_labels=(not args.no_labels),
                show_terminals=(not args.no_terminals),
            )
            if ok:
                n_ok += 1
            else:
                n_fail += 1

        print("Done. success=%d fail=%d" % (n_ok, n_fail))
        return 0 if n_fail == 0 else 1

    # ---------------- Single-file mode ----------------
    if args.report is None:
        ap.error("report is required unless --all is specified")

    out = args.out or (os.path.splitext(args.report)[0] + ".png")
    ok = _render_one(
        rpt_path=args.report,
        block_override=args.block or "",
        nets_mode=args.nets,
        out_path=out,
        show=(not args.no_show),
        show_labels=(not args.no_labels),
        show_terminals=(not args.no_terminals),
    )
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
