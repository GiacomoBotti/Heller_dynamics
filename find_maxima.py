#!/usr/bin/env python3
"""Find the strongest local maxima in a two-column dataset.

The input is assumed to be whitespace-separated numeric data.
Column numbers are 1-based, i.e. --xcol 1 means the first column.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Find and print the strongest local maxima of a y column."
    )
    parser.add_argument("file", type=Path, help="Input data file")
    parser.add_argument(
        "--xcol", type=int, default=1,
        help="Column containing x values (1-based, default: 1)",
    )
    parser.add_argument(
        "--ycol", type=int, default=2,
        help="Column containing y values (1-based, default: 2)",
    )
    parser.add_argument(
        "-n", "--number", type=int, default=5,
        help="Number of maxima to print (default: 5)",
    )
    parser.add_argument(
        "--xmin", type=float, default=None,
        help="Minimum x value to consider (inclusive)",
    )
    parser.add_argument(
        "--xmax", type=float, default=None,
        help="Maximum x value to consider (inclusive)",
    )
    parser.add_argument(
        "--ascending", action="store_true",
        help="Print selected maxima in ascending x order instead of descending y order",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    if args.xcol < 1 or args.ycol < 1:
        print("Error: --xcol and --ycol must be >= 1.", file=sys.stderr)
        return 2
    if args.number < 1:
        print("Error: --number must be >= 1.", file=sys.stderr)
        return 2
    if args.xmin is not None and args.xmax is not None and args.xmin > args.xmax:
        print("Error: --xmin must be <= --xmax.", file=sys.stderr)
        return 2

    try:
        data = np.loadtxt(args.file)
    except OSError as exc:
        print(f"Error reading {args.file}: {exc}", file=sys.stderr)
        return 1
    except ValueError as exc:
        print(f"Error parsing {args.file}: {exc}", file=sys.stderr)
        return 1

    if data.ndim == 1:
        data = data.reshape(1, -1)

    ncols = data.shape[1]
    if args.xcol > ncols or args.ycol > ncols:
        print(
            f"Error: file has {ncols} columns, but --xcol {args.xcol} "
            f"and/or --ycol {args.ycol} was requested.",
            file=sys.stderr,
        )
        return 2

    x = data[:, args.xcol - 1]
    y = data[:, args.ycol - 1]

    # Apply x-range before searching for maxima.
    mask = np.ones_like(x, dtype=bool)
    if args.xmin is not None:
        mask &= x >= args.xmin
    if args.xmax is not None:
        mask &= x <= args.xmax

    x_sel = x[mask]
    y_sel = y[mask]

    if len(x_sel) < 3:
        print("Error: at least 3 points are needed to identify local maxima.", file=sys.stderr)
        return 1

    # A strict local maximum: y[i] > y[i-1] and y[i] > y[i+1].
    peak_mask = (y_sel[1:-1] > y_sel[:-2]) & (y_sel[1:-1] > y_sel[2:])
    peak_indices = np.flatnonzero(peak_mask) + 1

    if len(peak_indices) == 0:
        print("No local maxima found in the selected x-range.")
        return 0

    # Keep the strongest peaks by y value.
    peak_indices = peak_indices[np.argsort(y_sel[peak_indices])[::-1]]
    peak_indices = peak_indices[: args.number]

    if args.ascending:
        peak_indices = peak_indices[np.argsort(x_sel[peak_indices])]

    print(f"{'x':>18} {'y':>18}")
    print(f"{'-' * 18} {'-' * 18}")
    for i in peak_indices:
        print(f"{x_sel[i]:18.10g} {y_sel[i]:18.10g}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
