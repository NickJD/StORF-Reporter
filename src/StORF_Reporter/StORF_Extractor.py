"""Legacy ``StORF-Extractor`` command with strict GFF3/FASTA handling."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .genetic_code import translate
from .io import InputFormatError, feature_line, write_fasta_records
from .utils import parse_bool as _bool
from .workflow import (
    atomic_output,
    discover_pairs,
    extract_feature_sequence,
    load_annotated_genome,
    logical_stem,
)


def _is_storf(feature) -> bool:
    identifier = str(feature.attributes.get("ID", ""))
    return feature.source == "StORF-Reporter" or identifier.startswith("StORF_") or "_StORF_" in identifier


def _extract_one(
    options: argparse.Namespace, gff: Path, fasta: Path | None, number: int
) -> tuple[Path, Path | None, int]:
    document, sequences = load_annotated_genome(gff, fasta)
    selected = [feature for feature in document.features if _is_storf(feature)]
    if not selected:
        raise InputFormatError(f"no StORF-Reporter features found in {gff}")
    output_dir = Path(options.o_dir).resolve() if options.o_dir else gff.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    if options.o_name:
        stem = options.o_name if number == 0 else f"{options.o_name}_{number + 1}"
    else:
        stem = f"{logical_stem(gff)}_Extracted_StORFs"
    suffix = ".gz" if options.gz else ""
    fasta_out = output_dir / f"{stem}{'_aa' if options.translate else ''}.fasta{suffix}"
    gff_out = output_dir / f"{stem}.gff{suffix}" if options.gff_out else None
    records = []
    for feature in selected:
        region = options.sequence_region
        if region == "auto":
            region = "candidate"
        sequence = extract_feature_sequence(feature, sequences[feature.seqid], region)
        table = int(feature.attributes.get("GeneticCode", options.genetic_code))
        if options.translate:
            if region == "candidate" and not options.stop_inclusive:
                sequence = sequence[3:]
            sequence = translate(
                sequence,
                table,
                keep_terminal_stop=options.stop_ident,
                start_as_methionine=region == "coding",
            )
            if not options.stop_ident:
                sequence = sequence.replace("*", "")
        elif region == "candidate" and not options.stop_inclusive:
            sequence = sequence[3:]
        identifier = str(feature.attributes.get("ID"))
        records.append((identifier, sequence))
    with atomic_output(fasta_out, options.overwrite) as handle:
        write_fasta_records(handle, records, 60 if options.line_wrap else 0)
    if gff_out is not None:
        with atomic_output(gff_out, options.overwrite) as handle:
            handle.write("##gff-version 3\n")
            for seqid, sequence in sequences.items():
                handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
            for feature in selected:
                handle.write(feature_line(feature) + "\n")
    return fasta_out, gff_out, len(selected)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="StORF-Extractor",
        description=f"StORF-Reporter {StORF_Reporter_Version}: extract reported candidate sequences",
    )
    parser.add_argument("-storf_input", choices=("Combined", "Separate"), required=True)
    parser.add_argument("-p", dest="path", required=True)
    parser.add_argument("-gff_out", dest="gff_out", type=_bool, default=False)
    parser.add_argument("-oname", dest="o_name")
    parser.add_argument("-odir", dest="o_dir")
    parser.add_argument("-gz", dest="gz", type=_bool, default=False)
    parser.add_argument("-aa", dest="translate", type=_bool, default=False)
    parser.add_argument(
        "-spos",
        "--include-first-stop",
        dest="stop_inclusive",
        type=_bool,
        default=False,
        help="include the upstream delimiting stop in candidate-locus nucleotide or protein output",
    )
    parser.add_argument(
        "-sequence_region",
        "--sequence-region",
        dest="sequence_region",
        choices=("auto", "candidate", "coding", "feature"),
        default="auto",
        help=(
            "auto extracts the stop-derived candidate locus for nucleotide and protein output; "
            "coding must be selected explicitly to use CodingCoordinates; feature uses displayed "
            "GFF coordinates"
        ),
    )
    parser.add_argument(
        "-stop_ident",
        dest="stop_ident",
        type=_bool,
        default=False,
        help="retain '*' stop symbols in protein output (default: False)",
    )
    parser.add_argument("-lw", dest="line_wrap", type=_bool, default=True)
    parser.add_argument("--genetic-code", type=int, choices=(1, 4, 11), default=11)
    parser.add_argument("-verbose", dest="verbose", type=_bool, default=False)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("-v", "--version", action="version", version=StORF_Reporter_Version)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    try:
        path = Path(options.path).resolve()
        pairs = discover_pairs(path, recursive=path.is_dir())
        total = 0
        for number, (_, gff, fasta) in enumerate(pairs):
            fasta_out, gff_out, count = _extract_one(options, gff, fasta, number)
            total += count
            if options.verbose:
                print(f"{count} candidates: {fasta_out}" + (f"; {gff_out}" if gff_out else ""))
        if options.verbose:
            print(f"complete: {total} candidates extracted")
        return 0
    except (FileExistsError, FileNotFoundError, InputFormatError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
