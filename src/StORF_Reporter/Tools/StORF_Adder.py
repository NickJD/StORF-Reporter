"""Merge additional StORF features into a reference GFF3 without ORForise."""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = "StORF_Reporter.Tools"

from ..core import feature_matches_gene_ident, merge_intervals
from ..io import feature_line, read_fasta, read_gff, write_fasta_records
from ..workflow import atomic_output


def _order_key(feature, seqid_rank: dict[str, int]) -> tuple[int, int, int, str]:
    return (
        seqid_rank.get(feature.seqid, len(seqid_rank)),
        feature.interval.start,
        feature.interval.end,
        str(feature.attributes.get("ID", "")),
    )


def _write_merged_features(handle, reference, additions, sequences) -> None:
    """Insert additions by contig/coordinate while retaining a sorted reference verbatim."""

    seqid_rank = {seqid: number for number, seqid in enumerate(sequences)}
    additions = sorted(additions, key=lambda item: _order_key(item, seqid_rank))
    original_keys = [_order_key(feature, seqid_rank) for feature in reference.features]
    coordinate_sorted = all(left <= right for left, right in zip(original_keys, original_keys[1:]))
    header_written = False
    if not any(line.startswith("##gff-version") for line in reference.pre_fasta_lines):
        handle.write("##gff-version 3\n")
        header_written = True

    if not coordinate_sorted:
        for line in reference.pre_fasta_lines:
            if line.startswith("##gff-version"):
                if not header_written:
                    handle.write("##gff-version 3\n")
                    header_written = True
            elif line.startswith("#") or not line.strip():
                handle.write(line)
        if not header_written:
            handle.write("##gff-version 3\n")
        handle.write("# Original feature order normalized by sequence and coordinate.\n")
        combined = sorted(
            (*reference.features, *additions),
            key=lambda item: _order_key(item, seqid_rank),
        )
        for feature in combined:
            handle.write(feature_line(feature) + "\n")
        return

    feature_number = 0
    addition_number = 0
    for line in reference.pre_fasta_lines:
        if line.startswith("##gff-version"):
            if not header_written:
                handle.write("##gff-version 3\n")
                header_written = True
            continue
        if line.startswith("#") or not line.strip():
            handle.write(line)
            continue
        original = reference.features[feature_number]
        feature_number += 1
        while addition_number < len(additions) and _order_key(
            additions[addition_number], seqid_rank
        ) < _order_key(original, seqid_rank):
            handle.write(feature_line(additions[addition_number]) + "\n")
            addition_number += 1
        handle.write(line if line.endswith("\n") else line + "\n")
    if not header_written:
        handle.write("##gff-version 3\n")
    for feature in additions[addition_number:]:
        handle.write(feature_line(feature) + "\n")


def add_features(options: argparse.Namespace) -> tuple[int, int]:
    if options.overlap < 0:
        raise ValueError("maximum overlap must be non-negative")
    sequences = read_fasta(options.genome_DNA)
    reference = read_gff(options.reference_annotation)
    additional = read_gff(options.additional_annotation)
    occupied = defaultdict(list)
    for feature in reference.features:
        if feature_matches_gene_ident(feature, options.gene_ident):
            occupied[feature.seqid].append(feature.interval)
    occupied = defaultdict(
        list,
        {seqid: merge_intervals(intervals) for seqid, intervals in occupied.items()},
    )
    reference_ids = {
        str(feature.attributes["ID"]) for feature in reference.features if feature.attributes.get("ID")
    }
    addition_ids: set[str] = set()
    kept = []
    rejected = 0
    for feature in additional.features:
        if feature.seqid not in sequences or feature.interval.end > len(sequences[feature.seqid]):
            raise ValueError(f"additional feature is outside FASTA bounds: {feature.seqid}")
        identifier = str(feature.attributes.get("ID", ""))
        if identifier and (identifier in reference_ids or identifier in addition_ids):
            raise ValueError(f"duplicate feature ID while merging: {identifier}")
        if identifier:
            addition_ids.add(identifier)
        total_overlap = sum(
            feature.interval.intersection_length(interval) for interval in occupied[feature.seqid]
        )
        if total_overlap <= options.overlap:
            kept.append(feature)
        else:
            rejected += 1
    output = Path(options.output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with atomic_output(output, options.overwrite) as handle:
        _write_merged_features(handle, reference, kept, sequences)
        if reference.embedded_fasta:
            handle.write("##FASTA\n")
            write_fasta_records(handle, reference.embedded_fasta.items())
    return len(kept), rejected


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Add filtered StORF features to a reference GFF3")
    parser.add_argument("-dna", "--genome_DNA", required=True)
    parser.add_argument("-rt", "--reference_tool", help=argparse.SUPPRESS)
    parser.add_argument("-ref", "--reference_annotation", required=True)
    parser.add_argument(
        "-gi",
        "--gene_ident",
        default="CDS",
        help="exact feature types, Prokka/Bakta, or an attribute query such as ID=gene",
    )
    parser.add_argument("-at", "--additional_tool", default="StORF-Reporter", help=argparse.SUPPRESS)
    parser.add_argument("-add", "--additional_annotation", required=True)
    parser.add_argument("-olap", "--overlap", default=50, type=int)
    parser.add_argument("-o", "--output_file", required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    options = build_parser().parse_args(argv)
    kept, rejected = add_features(options)
    print(f"kept {kept} additional features; rejected {rejected} by overlap")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
