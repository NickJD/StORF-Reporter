"""Legacy ``StORF-Finder`` command backed by the v2 prediction core."""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from contextlib import ExitStack
from pathlib import Path
from typing import TextIO

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .core import apply_overlap_filter, find_storf_candidates
from .io import feature_line, read_fasta, write_fasta_records
from .models import Feature, Interval, StorfCandidate, UnannotatedRegion
from .utils import parse_bool as _bool
from .workflow import (
    AnalysisConfig,
    atomic_output,
    candidate_feature,
    candidate_nucleotide_sequence,
    candidate_protein_sequence,
    logical_stem,
)


def _genetic_code(options: argparse.Namespace) -> int:
    explicit = getattr(options, "genetic_code", None)
    if explicit is not None:
        return int(explicit)
    stops = frozenset(item.strip().upper() for item in options.stop_codons.split(","))
    if stops == frozenset({"TAA", "TAG"}):
        return 4
    if stops == frozenset({"TAA", "TAG", "TGA"}):
        return 11
    raise ValueError(
        "custom stop sets are no longer guessed; use the standard set or TAA,TAG (genetic code 4)"
    )


def _region(sequence_id: str, sequence: str, unannotated: bool) -> UnannotatedRegion:
    sequence = sequence.upper()
    seqid = sequence_id
    origin = 0
    if unannotated and "_UR_" in sequence_id:
        seqid, coordinates = sequence_id.rsplit("_UR_", 1)
        try:
            start_text, end_text = coordinates.split("_", 1)
            start, end = int(start_text), int(end_text)
            if end - start + 1 == len(sequence):
                origin = start - 1
        except (ValueError, TypeError):
            origin = 0
    interval = Interval(origin, origin + len(sequence))
    return UnannotatedRegion(
        id=sequence_id,
        seqid=seqid,
        core=interval,
        extended=interval,
        sequence=sequence,
    )


def _find(options: argparse.Namespace, sequence_id: str, sequence: str) -> list[StorfCandidate]:
    if getattr(options, "partial_storf", False) or getattr(options, "whole_contig", False):
        raise ValueError(
            "partial/run-through candidates are disabled because they are not stop-to-stop ORFs; "
            "use the v2 evidence table to review edge candidates in a future dedicated model"
        )
    region = _region(sequence_id, sequence, bool(getattr(options, "unannotated", True)))
    short_mode = getattr(options, "short_storfs", False)
    include_short = short_mode not in {False, "False", None} or bool(
        getattr(options, "short_storfs_only", False)
    )
    candidates, _ = find_storf_candidates(
        region,
        genetic_code=_genetic_code(options),
        min_orf=int(options.min_orf),
        max_orf=int(options.max_orf),
        ambiguous_threshold=float(options.non_standard),
        include_short=include_short,
        short_min=30,
        short_max=120,
        require_start=bool(getattr(options, "start_filtering", False)),
        include_consecutive=bool(
            getattr(options, "con_storfs", False) or getattr(options, "con_only", False)
        ),
        consecutive_only=bool(getattr(options, "con_only", False)),
    )
    if getattr(options, "short_storfs_only", False):
        candidates = [item for item in candidates if item.candidate_type == "Short-StORF"]
    overlap_mode = {
        "both-strand": "all",
        "single-strand": "same-strand",
        "none": "none",
    }[options.olap_filtering]
    coordinate_basis = getattr(options, "overlap_coordinate_basis", "coding")
    if short_mode == "Olap" and not getattr(options, "short_storfs_only", False):
        ordinary = [item for item in candidates if item.candidate_type != "Short-StORF"]
        short = [item for item in candidates if item.candidate_type == "Short-StORF"]
        ordinary, _ = apply_overlap_filter(ordinary, int(options.overlap_nt), overlap_mode, coordinate_basis)
        short, _ = apply_overlap_filter(short, int(options.overlap_nt), overlap_mode, coordinate_basis)
        candidates = sorted(
            (*ordinary, *short),
            key=lambda item: (item.seqid, item.interval.start, item.interval.end, item.strand),
        )
    else:
        candidates, _ = apply_overlap_filter(
            candidates,
            int(options.overlap_nt),
            overlap_mode,
            coordinate_basis,
        )
    if getattr(options, "storf_order", "start_pos") == "strand":
        candidates.sort(
            key=lambda item: (item.seqid, item.strand == "-", item.interval.start, item.interval.end)
        )
    return candidates


def _legacy_dict(candidates: list[StorfCandidate]) -> OrderedDict[str, list[object]]:
    result: OrderedDict[str, list[object]] = OrderedDict()
    for index, candidate in enumerate(candidates):
        key = f"{candidate.interval.start},{candidate.interval.end}"
        result[key] = [
            candidate.sequence,
            str(candidate.frame),
            candidate.strand,
            candidate.length,
            candidate.candidate_type,
            index,
        ]
    return result


def _display_feature(candidate: StorfCandidate, options: argparse.Namespace) -> Feature:
    requested = getattr(options, "feature_type", "ORF")
    config = AnalysisConfig(
        feature_type="ORF" if requested == "CDS" else requested,
        gff_coordinate_mode=getattr(options, "gff_coordinate_mode", "stop-to-stop"),
        include_first_stop=bool(getattr(options, "stop_inclusive", False)),
        promote_tier="exploratory" if requested == "CDS" else None,
    )
    return candidate_feature(candidate, config)


def _write_candidates(
    candidates: list[StorfCandidate],
    options: argparse.Namespace,
    fasta_out: TextIO | None,
    aa_fasta_out: TextIO | None,
    gff_out: TextIO | None,
) -> None:
    nt_records = []
    aa_records = []
    for candidate in candidates:
        sequence = candidate_nucleotide_sequence(
            candidate,
            include_leading_stop=bool(getattr(options, "stop_inclusive", False)),
        )
        nt_records.append((candidate.id, sequence))
        aa = candidate_protein_sequence(
            candidate,
            coordinate_mode=getattr(options, "protein_coordinate_mode", "candidate"),
            include_leading_stop=bool(getattr(options, "stop_inclusive", False)),
            keep_terminal_stop=bool(getattr(options, "stop_ident", False)),
            remove_inframe_stops=bool(getattr(options, "remove_stop", False)),
        )
        aa_records.append((candidate.id, aa))
        if gff_out is not None:
            gff_out.write(feature_line(_display_feature(candidate, options)) + "\n")
    if fasta_out is not None:
        write_fasta_records(fasta_out, nt_records, 60 if getattr(options, "line_wrap", True) else 0)
    if aa_fasta_out is not None:
        write_fasta_records(aa_fasta_out, aa_records, 60 if getattr(options, "line_wrap", True) else 0)


def STORF_Finder(
    options: argparse.Namespace,
    sequence_info,
    sequence_id: str,
    fasta_out,
    aa_fasta_out,
    gff_out,
    split_index: int = 3,
):
    """Historical programmatic entry point, now routed through the v2 core.

    ``split_index`` remains accepted for source compatibility; v2 uses typed
    candidates rather than the positional records that required that offset.
    """

    _ = split_index
    sequence = sequence_info[1] if isinstance(sequence_info, (list, tuple)) else str(sequence_info)
    candidates = _find(options, sequence_id.lstrip(">"), sequence)
    if bool(getattr(options, "reporter", False)):
        return _legacy_dict(candidates)
    _write_candidates(candidates, options, fasta_out, aa_fasta_out, gff_out)
    return None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="StORF-Finder",
        description=f"StORF-Reporter {StORF_Reporter_Version}: find stop-to-stop ORF candidates",
    )
    parser.add_argument("-f", dest="fasta", required=True)
    parser.add_argument("-ua", dest="unannotated", type=_bool, default=True)
    parser.add_argument("-wc", dest="whole_contig", type=_bool, default=False)
    parser.add_argument("-ps", dest="partial_storf", type=_bool, default=False)
    parser.add_argument(
        "-olap_filt",
        dest="olap_filtering",
        choices=("none", "single-strand", "both-strand"),
        default="both-strand",
    )
    parser.add_argument(
        "--overlap-coordinate-basis",
        choices=("coding", "candidate"),
        default="coding",
        help="measure overlap using ranked coding or full stop-to-stop intervals",
    )
    parser.add_argument("-start_filt", dest="start_filtering", type=_bool, default=False)
    parser.add_argument("-con_storfs", dest="con_storfs", type=_bool, default=False)
    parser.add_argument("-con_only", dest="con_only", type=_bool, default=False)
    parser.add_argument(
        "-short_storfs", dest="short_storfs", choices=("False", "Nolap", "Olap"), default="False"
    )
    parser.add_argument("-short_storfs_only", dest="short_storfs_only", type=_bool, default=False)
    parser.add_argument(
        "-f_type",
        dest="feature_type",
        choices=("StORF", "CDS", "ORF", "open_reading_frame"),
        default="ORF",
    )
    parser.add_argument(
        "--gff-coordinate-mode",
        choices=("stop-to-stop", "start-to-stop"),
        default="stop-to-stop",
        help="GFF coordinates for unpromoted candidates (default: stop-to-stop)",
    )
    parser.add_argument("-minorf", dest="min_orf", type=int, default=99)
    parser.add_argument("-maxorf", dest="max_orf", type=int, default=60_000)
    parser.add_argument("-codons", dest="stop_codons", default="TAG,TGA,TAA")
    parser.add_argument("--genetic-code", type=int, choices=(1, 4, 11))
    parser.add_argument("-non_standard", dest="non_standard", type=float, default=0.20)
    parser.add_argument("-olap", dest="overlap_nt", type=int, default=50)
    parser.add_argument("-so", dest="storf_order", choices=("start_pos", "strand"), default="start_pos")
    parser.add_argument("-oname", dest="o_name")
    parser.add_argument("-odir", dest="o_dir")
    parser.add_argument("-gff", dest="gff", type=_bool, default=True)
    parser.add_argument("-aa", dest="translate", type=_bool, default=False)
    parser.add_argument("-aa_only", dest="aa_only", type=_bool, default=False)
    parser.add_argument("-lw", dest="line_wrap", type=_bool, default=True)
    parser.add_argument(
        "-spos",
        dest="stop_inclusive",
        type=_bool,
        default=False,
        help=("include the upstream delimiting stop in GFF3 and candidate-locus FASTA (default: False)"),
    )
    parser.add_argument("-stop_ident", dest="stop_ident", type=_bool, default=True)
    parser.add_argument(
        "-rs",
        dest="remove_stop",
        type=_bool,
        default=True,
        help="remove internal '*' stop symbols from candidate protein FASTA (default: True)",
    )
    parser.add_argument(
        "--protein-coordinate-mode",
        choices=("candidate", "coding"),
        default="candidate",
        help=(
            "translate the reported stop-derived candidate locus by default; use coding only to "
            "request ranked CodingCoordinates"
        ),
    )
    parser.add_argument("-gz", dest="gz", type=_bool, default=False)
    parser.add_argument("-verbose", dest="verbose", type=_bool, default=False)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("-v", "--version", action="version", version=StORF_Reporter_Version)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    try:
        sequences = read_fasta(options.fasta)
        input_path = Path(options.fasta).resolve()
        output_dir = Path(options.o_dir).resolve() if options.o_dir else input_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = options.o_name or f"{logical_stem(input_path)}_StORF-Finder"
        prefix = output_dir / stem
        suffix = ".gz" if options.gz else ""
        paths = {
            "gff": Path(f"{prefix}.gff{suffix}"),
            "nt": Path(f"{prefix}.fasta{suffix}"),
            "aa": Path(f"{prefix}_aa.fasta{suffix}"),
        }
        candidates_by_sequence = {
            seqid: _find(options, seqid, sequence) for seqid, sequence in sequences.items()
        }
        selected_paths = []
        if options.gff and not options.aa_only:
            selected_paths.append(paths["gff"])
        if not options.aa_only:
            selected_paths.append(paths["nt"])
        if options.translate or options.aa_only:
            selected_paths.append(paths["aa"])
        existing = [path for path in selected_paths if path.exists()]
        if existing and not options.overwrite:
            raise FileExistsError(f"output exists; use --overwrite: {existing[0]}")
        with ExitStack() as stack:
            gff_out = None
            nt_out = None
            aa_out = None
            if options.gff and not options.aa_only:
                gff_out = stack.enter_context(atomic_output(paths["gff"], options.overwrite))
                gff_out.write("##gff-version 3\n")
                for seqid, sequence in sequences.items():
                    gff_out.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
            if not options.aa_only:
                nt_out = stack.enter_context(atomic_output(paths["nt"], options.overwrite))
            if options.translate or options.aa_only:
                aa_out = stack.enter_context(atomic_output(paths["aa"], options.overwrite))
            for candidates in candidates_by_sequence.values():
                _write_candidates(candidates, options, nt_out, aa_out, gff_out)
        count = sum(map(len, candidates_by_sequence.values()))
        if options.verbose:
            print(f"{count} candidates found")
        return 0
    except (FileExistsError, FileNotFoundError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
