"""Summarise FASTA sequence lengths without NumPy."""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = "StORF_Reporter.Tools"

from ..io import read_fasta


def _percentile(values: list[int], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1 - weight) + ordered[upper] * weight


def seq_Lengths(URs: str) -> dict[str, float]:
    lengths = [len(sequence) for sequence in read_fasta(URs).values()]
    return {
        "count": float(len(lengths)),
        "longest": float(max(lengths)),
        "shortest": float(min(lengths)),
        "median": float(statistics.median(lengths)),
        "mean": float(statistics.mean(lengths)),
        "std": float(statistics.pstdev(lengths)),
        "p25": _percentile(lengths, 0.25),
        "p75": _percentile(lengths, 0.75),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Summarise UR FASTA lengths")
    parser.add_argument("-ur", "--URs", required=True)
    args = parser.parse_args(argv)
    for name, value in seq_Lengths(args.URs).items():
        print(f"{name}\t{value:.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
