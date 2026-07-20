"""High-level analysis and reproducible output generation."""

from __future__ import annotations

import csv
import hashlib
import html
import json
import os
from collections import OrderedDict
from collections.abc import Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import TextIO

from .constants import StORF_Reporter_Version
from .core import (
    TIER_ORDER,
    apply_overlap_filter,
    extract_unannotated_regions,
    find_storf_candidates,
)
from .genetic_code import reverse_complement, translate
from .io import (
    GffDocument,
    InputFormatError,
    feature_line,
    open_text,
    read_fasta,
    read_gff,
    write_fasta_records,
)
from .models import Feature, Interval, RunSummary, StorfCandidate, UnannotatedRegion
from .validation import validate_candidate_cds, validate_con_storf_candidate


def logical_stem(path: str | Path) -> str:
    """Return a biological filename stem after compression/data suffixes."""

    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".gff3", ".gff", ".fasta", ".fna", ".fa"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def storf_query_id(value: str) -> str:
    """Return the GFF3 ID represented by a modern or historical FASTA query name.

    Historical StORF FASTA files placed ``;UR=...`` and other GFF-like fields
    directly in the identifier (without whitespace).  BLAST therefore retains
    the whole header in qseqid, whereas the corresponding GFF3 ``ID`` ends at
    the first semicolon.
    """

    return value.strip().lstrip(">").split(";", 1)[0]


def is_gff_path(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return name.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz"))


def is_fasta_path(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return name.endswith((".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz"))


def feature_attribute_interval(feature: Feature, attribute: str) -> Interval | None:
    """Parse one of the explicit v2 coordinate attributes on a feature.

    Candidate and coding coordinates are genomic, one-based, inclusive values,
    regardless of strand.  Returning ``None`` for the literal ``none`` keeps
    start-less candidates representable without guessing a coding interval.
    """

    text = str(feature.attributes.get(attribute, "")).strip()
    if not text or text.casefold() == "none":
        return None
    try:
        start_text, end_text = text.split("-", 1)
        return Interval.from_gff(int(start_text), int(end_text))
    except (TypeError, ValueError) as error:
        identifier = feature.attributes.get("ID", "StORF feature")
        raise InputFormatError(
            f"{identifier} has invalid {attribute}={text!r}; expected START-END"
        ) from error


def storf_sequence_interval(feature: Feature, region: str) -> Interval:
    """Resolve the candidate, coding, or displayed interval of a StORF feature."""

    if region == "feature":
        return feature.interval
    if region == "candidate":
        return feature_attribute_interval(feature, "CandidateCoordinates") or feature.interval
    if region == "coding":
        interval = feature_attribute_interval(feature, "CodingCoordinates")
        if interval is None and feature.type == "CDS":
            # Compatibility for pre-v2 CDS records without dual coordinates.
            interval = feature.interval
        if interval is None:
            identifier = feature.attributes.get("ID", "StORF feature")
            raise InputFormatError(f"{identifier} has no supported coding interval")
        return interval
    raise ValueError("sequence region must be candidate, coding, or feature")


def extract_feature_sequence(feature: Feature, genomic_sequence: str, region: str) -> str:
    """Extract an oriented StORF sequence using an explicit coordinate representation."""

    interval = storf_sequence_interval(feature, region)
    if interval.end > len(genomic_sequence):
        identifier = feature.attributes.get("ID", "StORF feature")
        raise InputFormatError(f"{identifier} {region} coordinates exceed sequence {feature.seqid!r}")
    sequence = genomic_sequence[interval.start : interval.end]
    return reverse_complement(sequence) if feature.strand == "-" else sequence


def candidate_nucleotide_sequence(
    candidate: StorfCandidate,
    include_leading_stop: bool = False,
) -> str:
    """Return the oriented candidate locus using the historical first-stop policy."""

    return candidate.sequence if include_leading_stop else candidate.sequence[3:]


def candidate_protein_sequence(
    candidate: StorfCandidate,
    coordinate_mode: str = "candidate",
    include_leading_stop: bool = False,
    keep_terminal_stop: bool = False,
    remove_inframe_stops: bool = False,
) -> str:
    """Translate an explicitly selected candidate or ranked coding representation."""

    candidate_sequence = candidate_nucleotide_sequence(
        candidate,
        include_leading_stop=include_leading_stop,
    )
    if coordinate_mode == "candidate":
        sequence = candidate_sequence
        start_as_methionine = False
    elif coordinate_mode == "coding":
        sequence = candidate.coding_sequence or candidate_sequence
        start_as_methionine = candidate.coding_sequence is not None
    else:
        raise ValueError("protein coordinate mode must be candidate or coding")
    protein = translate(
        sequence,
        candidate.genetic_code,
        keep_terminal_stop=keep_terminal_stop,
        start_as_methionine=start_as_methionine,
    )
    return protein.replace("*", "") if remove_inframe_stops else protein


@dataclass(frozen=True)
class AnalysisConfig:
    excluded_types: tuple[str, ...] = ("CDS",)
    gene_ident: str | None = None
    min_ur: int = 30
    max_ur: int = 100_000
    extension: int = 50
    min_orf: int = 99
    max_orf: int = 60_000
    genetic_code: int = 11
    ambiguous_threshold: float = 0.20
    include_short: bool = False
    short_only: bool = False
    short_overlap_mode: str = "exclude"
    short_min: int = 30
    short_max: int = 120
    require_start: bool = False
    include_consecutive: bool = False
    consecutive_only: bool = False
    overlap_mode: str = "all"
    overlap_coordinate_basis: str = "coding"
    maximum_overlap: int = 50
    maximum_annotation_overlap: int = 50
    feature_type: str = "ORF"
    gff_coordinate_mode: str = "stop-to-stop"
    include_first_stop: bool = False
    protein_coordinate_mode: str = "candidate"
    remove_inframe_stops: bool = False
    promote_tier: str | None = None
    primary_tier: str = "exploratory"


@dataclass
class AnalysisResult:
    sequences: OrderedDict[str, str]
    features: list[Feature]
    regions: list[UnannotatedRegion]
    candidates: list[StorfCandidate]
    summary: RunSummary
    document: GffDocument | None = None
    input_gff: Path | None = None
    input_fasta: Path | None = None
    sample: str = "sample"


def _normalise_types(types: Iterable[str]) -> tuple[str, ...]:
    result = tuple(dict.fromkeys(item.strip() for item in types if item.strip()))
    if not result:
        raise ValueError("at least one excluded feature type is required")
    return result


def read_external_evidence(path: str | Path | None) -> dict[str, dict[str, float]]:
    if path is None:
        return {}
    evidence: dict[str, dict[str, float]] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "candidate_id" not in reader.fieldnames:
            raise InputFormatError("evidence TSV requires a candidate_id column")
        for row_number, row in enumerate(reader, 2):
            candidate_id = (row.get("candidate_id") or "").strip()
            if not candidate_id:
                raise InputFormatError(f"missing candidate_id in evidence TSV row {row_number}")
            values: dict[str, float] = {}
            for name, value in row.items():
                if name == "candidate_id" or value in (None, ""):
                    continue
                try:
                    numeric = float(value)
                except ValueError as error:
                    raise InputFormatError(
                        f"evidence value {name}={value!r} in row {row_number} is not numeric"
                    ) from error
                if not 0 <= numeric <= 1:
                    raise InputFormatError(
                        f"evidence value {name} in row {row_number} must be between 0 and 1"
                    )
                values[name] = numeric
            evidence[candidate_id] = values
    return evidence


def load_annotated_genome(
    gff_path: str | Path,
    fasta_path: str | Path | None = None,
) -> tuple[GffDocument, OrderedDict[str, str]]:
    document = read_gff(gff_path)
    if fasta_path is not None:
        sequences = read_fasta(fasta_path)
    elif document.embedded_fasta is not None:
        sequences = document.embedded_fasta
    else:
        raise InputFormatError(f"{gff_path} has no embedded FASTA; provide a matching FASTA file")
    missing = sorted({feature.seqid for feature in document.features} - set(sequences))
    if missing:
        preview = ", ".join(missing[:5])
        raise InputFormatError(f"GFF3 sequence IDs absent from FASTA: {preview}")
    return document, sequences


def pyrodigal_features(
    sequences: Mapping[str, str],
    training: str = "longest",
) -> list[Feature]:
    try:
        import pyrodigal
    except ImportError as error:
        raise RuntimeError("Pyrodigal mode requires the pyrodigal package") from error
    if training not in {"longest", "individual", "meta"}:
        raise ValueError("Pyrodigal training must be longest, individual or meta")

    shared_finder = None
    if training == "meta":
        shared_finder = pyrodigal.GeneFinder(meta=True)
    elif training == "longest":
        longest = max(sequences.values(), key=len)
        if len(longest) < 20_000:
            shared_finder = pyrodigal.GeneFinder(meta=True)
        else:
            shared_finder = pyrodigal.GeneFinder()
            shared_finder.train(longest)

    features: list[Feature] = []
    for seqid, sequence in sequences.items():
        if training == "individual":
            if len(sequence) < 20_000:
                finder = pyrodigal.GeneFinder(meta=True)
            else:
                finder = pyrodigal.GeneFinder()
                finder.train(sequence)
        else:
            finder = shared_finder
        assert finder is not None
        genes = finder.find_genes(sequence)
        for number, gene in enumerate(genes, 1):
            confidence = gene.confidence() if callable(gene.confidence) else gene.confidence
            attributes = {
                "ID": f"pyrodigal_{seqid}_{number}",
                "Name": f"pyrodigal_{number}",
                "partial": f"{int(gene.partial_begin)}{int(gene.partial_end)}",
                "start_type": str(gene.start_type),
                "rbs_motif": str(gene.rbs_motif),
                "confidence": f"{float(confidence):.3f}",
                "transl_table": str(gene.translation_table),
            }
            features.append(
                Feature(
                    seqid=seqid,
                    source="Pyrodigal",
                    type="CDS",
                    interval=Interval.from_gff(int(gene.begin), int(gene.end)),
                    score=f"{float(gene.score):.3f}",
                    strand="+" if int(gene.strand) == 1 else "-",
                    phase="0",
                    attributes=attributes,
                )
            )
    return features


def analyse(
    sequences: OrderedDict[str, str],
    features: list[Feature],
    config: AnalysisConfig,
    external_evidence: Mapping[str, Mapping[str, float]] | None = None,
    document: GffDocument | None = None,
    input_gff: Path | None = None,
    input_fasta: Path | None = None,
    sample: str = "sample",
) -> AnalysisResult:
    excluded = _normalise_types(config.excluded_types) if config.gene_ident is None else config.excluded_types
    regions = extract_unannotated_regions(
        sequences,
        features,
        excluded_types=excluded,
        min_length=config.min_ur,
        max_length=config.max_ur,
        extension=config.extension,
        gene_ident=config.gene_ident,
    )
    raw_candidates: list[StorfCandidate] = []
    rejected_ambiguous = 0
    for region in regions:
        found, rejected = find_storf_candidates(
            region,
            genetic_code=config.genetic_code,
            min_orf=config.min_orf,
            max_orf=config.max_orf,
            ambiguous_threshold=config.ambiguous_threshold,
            include_short=config.include_short,
            short_min=config.short_min,
            short_max=config.short_max,
            require_start=config.require_start,
            include_consecutive=config.include_consecutive,
            consecutive_only=config.consecutive_only,
            external_evidence=external_evidence,
        )
        raw_candidates.extend(found)
        rejected_ambiguous += rejected
    if config.short_only:
        raw_candidates = [
            candidate for candidate in raw_candidates if candidate.candidate_type == "Short-StORF"
        ]
    annotation_eligible = [
        candidate
        for candidate in raw_candidates
        if candidate.evidence.overlap_annotated_nt <= config.maximum_annotation_overlap
    ]
    annotation_removed = len(raw_candidates) - len(annotation_eligible)
    if config.short_overlap_mode not in {"exclude", "allow"}:
        raise ValueError("short-StORF overlap mode must be exclude or allow")
    if config.short_overlap_mode == "allow" and not config.short_only:
        ordinary = [
            candidate for candidate in annotation_eligible if candidate.candidate_type != "Short-StORF"
        ]
        short = [candidate for candidate in annotation_eligible if candidate.candidate_type == "Short-StORF"]
        ordinary, ordinary_removed = apply_overlap_filter(
            ordinary,
            maximum_overlap=config.maximum_overlap,
            mode=config.overlap_mode,
            coordinate_basis=config.overlap_coordinate_basis,
        )
        short, short_removed = apply_overlap_filter(
            short,
            maximum_overlap=config.maximum_overlap,
            mode=config.overlap_mode,
            coordinate_basis=config.overlap_coordinate_basis,
        )
        candidates = sorted(
            (*ordinary, *short),
            key=lambda item: (item.seqid, item.interval.start, item.interval.end, item.strand),
        )
        removed = ordinary_removed + short_removed
    else:
        candidates, removed = apply_overlap_filter(
            annotation_eligible,
            maximum_overlap=config.maximum_overlap,
            mode=config.overlap_mode,
            coordinate_basis=config.overlap_coordinate_basis,
        )
    summary = RunSummary(
        sequences=len(sequences),
        annotations=len(features),
        unannotated_regions=len(regions),
        candidates_before_filtering=len(raw_candidates),
        candidates=len(candidates),
        filtered_overlap=removed,
        filtered_annotation_overlap=annotation_removed,
        rejected_ambiguous=rejected_ambiguous,
    )
    for candidate in candidates:
        summary.tiers[candidate.evidence.tier] += 1
    return AnalysisResult(
        sequences=sequences,
        features=features,
        regions=regions,
        candidates=candidates,
        summary=summary,
        document=document,
        input_gff=input_gff,
        input_fasta=input_fasta,
        sample=sample,
    )


def candidate_feature(candidate: StorfCandidate, config: AnalysisConfig) -> Feature:
    feature_type = config.feature_type
    if feature_type == "CDS":
        raise ValueError("CDS is reserved for validated promotion; set promote_tier instead")
    if config.gff_coordinate_mode not in {"stop-to-stop", "start-to-stop"}:
        raise ValueError("GFF coordinate mode must be stop-to-stop or start-to-stop")
    interval = candidate.interval
    phase = "."
    promoted = False
    minimum_coding_length = max(
        6,
        (config.short_min if candidate.candidate_type == "Short-StORF" else config.min_orf) - 3,
    )
    cds_validation = validate_candidate_cds(candidate, minimum_coding_length)
    if candidate.candidate_type == "Con-StORF":
        candidate_validation = validate_con_storf_candidate(candidate, minimum_coding_length)
        candidate_validation_policy = "con-storf-readthrough-v1"
    else:
        candidate_validation = cds_validation
        candidate_validation_policy = "complete-orf-v1"
    if config.promote_tier is not None:
        if config.promote_tier not in TIER_ORDER:
            raise ValueError("promotion tier must be high, medium or exploratory")
        if (
            candidate.candidate_type != "Con-StORF"
            and TIER_ORDER[candidate.evidence.tier] >= TIER_ORDER[config.promote_tier]
            and cds_validation.passed
        ):
            feature_type = "CDS"
            assert candidate.coding_interval is not None
            interval = candidate.coding_interval
            phase = "0"
            promoted = True
    reported_coordinate_mode = "stop-to-stop"
    if config.gff_coordinate_mode == "start-to-stop" and candidate.coding_interval is not None:
        interval = candidate.coding_interval
        reported_coordinate_mode = "start-to-stop"
    else:
        if config.include_first_stop:
            reported_coordinate_mode = "stop-to-stop-inclusive"
        elif candidate.strand == "+":
            interval = Interval(candidate.interval.start + 3, candidate.interval.end)
            reported_coordinate_mode = "stop-to-stop-first-stop-excluded"
        else:
            interval = Interval(candidate.interval.start, candidate.interval.end - 3)
            reported_coordinate_mode = "stop-to-stop-first-stop-excluded"

    middle_stop_coordinates: list[str] = []
    for offset in candidate.middle_stop_offsets:
        if candidate.strand == "+":
            middle = Interval(candidate.interval.start + offset, candidate.interval.start + offset + 3)
        else:
            middle = Interval(candidate.interval.end - offset - 3, candidate.interval.end - offset)
        start, end = middle.to_gff()
        middle_stop_coordinates.append(f"{start}-{end}")
    attributes: dict[str, object] = {
        "ID": candidate.id,
        "Name": candidate.id,
        "UR": candidate.region_id,
        "StORF_Type": candidate.candidate_type,
        "CandidateCoordinates": f"{candidate.interval.start + 1}-{candidate.interval.end}",
        "CodingCoordinates": (
            f"{candidate.coding_interval.start + 1}-{candidate.coding_interval.end}"
            if candidate.coding_interval is not None
            else "none"
        ),
        "CandidateLength": candidate.length,
        "CodingLength": candidate.coding_length if candidate.coding_length is not None else "none",
        "Frame": candidate.frame,
        "GeneticCode": candidate.genetic_code,
        "StartStop": candidate.start_stop,
        "EndStop": candidate.end_stop,
        "StartCodon": candidate.evidence.start_codon or "none",
        "RBS": candidate.evidence.rbs_motif or "none",
        "RBSDistance": (
            candidate.evidence.rbs_distance if candidate.evidence.rbs_distance is not None else "none"
        ),
        "StartSelectionPolicy": "ranked-rbs-start-codon-length-v1",
        "StartSiteScore": (
            f"{candidate.evidence.start_site_score:.3f}"
            if candidate.evidence.start_site_score is not None
            else "none"
        ),
        "StartCandidateCount": candidate.evidence.start_candidate_count,
        "MiddleStopCount": len(candidate.middle_stop_offsets),
        "MiddleStops": ",".join(candidate.middle_stop_codons) or "none",
        "MiddleStopCoordinates": ",".join(middle_stop_coordinates) or "none",
        "StartBeforeFirstMiddleStop": (
            str(candidate.start_before_first_middle_stop).lower()
            if candidate.candidate_type == "Con-StORF"
            else "not_applicable"
        ),
        "EvidenceTier": candidate.evidence.tier,
        "EvidenceScore": f"{candidate.evidence.score:.3f}",
        "ScoreModel": "heuristic-v1-not-a-probability",
        "CoreFraction": f"{candidate.evidence.core_fraction:.4f}",
        "AnnotationOverlap": candidate.evidence.overlap_annotated_nt,
        "AmbiguousFraction": f"{candidate.evidence.ambiguous_fraction:.4f}",
        "CandidateValidation": "passed" if candidate_validation.passed else "failed",
        "CandidateValidationPolicy": candidate_validation_policy,
        "CandidateValidationIssues": (
            "none" if candidate_validation.passed else " | ".join(candidate_validation.issues)
        ),
        "CDSValidation": "passed" if cds_validation.passed else "failed",
        "CDSValidationPolicy": "complete-orf-v1",
        "CDSValidationIssues": "none" if cds_validation.passed else " | ".join(cds_validation.issues),
        "PromotedToCDS": str(promoted).lower(),
        "GFFCoordinateMode": reported_coordinate_mode,
    }
    for name, value in sorted(candidate.evidence.external.items()):
        attributes[f"Evidence_{name}"] = f"{value:.3f}"
    return Feature(
        seqid=candidate.seqid,
        source="StORF-Reporter",
        type=feature_type,
        interval=interval,
        score=f"{candidate.evidence.score:.3f}",
        strand=candidate.strand,
        phase=phase,
        attributes=attributes,
    )


def partition_candidates(
    candidates: Sequence[StorfCandidate],
    minimum_tier: str,
) -> tuple[list[StorfCandidate], list[StorfCandidate]]:
    """Split candidates into primary and lower-tier review outputs."""

    if minimum_tier not in TIER_ORDER:
        raise ValueError("primary tier must be high, medium or exploratory")
    threshold = TIER_ORDER[minimum_tier]
    primary = [candidate for candidate in candidates if TIER_ORDER[candidate.evidence.tier] >= threshold]
    review = [candidate for candidate in candidates if TIER_ORDER[candidate.evidence.tier] < threshold]
    return primary, review


@contextmanager
def atomic_output(path: str | Path, overwrite: bool = False) -> Iterator[TextIO]:
    path = Path(path)
    if path.exists() and not overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.stem}.tmp{path.suffix}")
    try:
        with open_text(temporary, "wt") as handle:
            yield handle
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def ensure_outputs_available(paths: Iterable[str | Path], overwrite: bool = False) -> None:
    if overwrite:
        return
    existing = [Path(path) for path in paths if Path(path).exists()]
    if existing:
        raise FileExistsError(f"output exists; use --overwrite: {existing[0]}")


def _normalised_original_lines(document: GffDocument) -> list[str]:
    result: list[str] = []
    saw_version = False
    for line in document.pre_fasta_lines:
        if line.startswith("##gff-version"):
            if not saw_version:
                result.append("##gff-version 3\n")
                saw_version = True
            continue
        result.append(line if line.endswith("\n") else line + "\n")
    if not saw_version:
        result.insert(0, "##gff-version 3\n")
    return result


def _feature_order_key(feature: Feature, seqid_rank: Mapping[str, int]) -> tuple[int, int, int, str]:
    return (
        seqid_rank.get(feature.seqid, len(seqid_rank)),
        feature.interval.start,
        feature.interval.end,
        str(feature.attributes.get("ID", "")),
    )


def _write_storf_notice(handle: TextIO, config: AnalysisConfig) -> None:
    handle.write(f"# StORF-Reporter {StORF_Reporter_Version}\n")
    handle.write("# EvidenceScore is a prioritisation heuristic, not a probability.\n")
    if config.gene_ident is not None:
        gene_ident = config.gene_ident.replace("\r", " ").replace("\n", " ")
        handle.write(f"# GeneIdent={gene_ident}\n")
    handle.write(f"# PrimaryMinimumTier={config.primary_tier}\n")


def _write_candidate_feature_lines(handle: TextIO, features: Sequence[Feature]) -> None:
    if not features:
        return
    for feature in features:
        handle.write(feature_line(feature) + "\n")
    handle.write("###\n")


def _write_interleaved_features(
    handle: TextIO,
    result: AnalysisResult,
    config: AnalysisConfig,
    candidates: Sequence[StorfCandidate],
) -> None:
    """Preserve a sorted input document while inserting candidates by locus."""

    assert result.document is not None
    lines = _normalised_original_lines(result.document)
    seqid_rank = {seqid: index for index, seqid in enumerate(result.sequences)}
    original_keys = [_feature_order_key(feature, seqid_rank)[:2] for feature in result.document.features]
    coordinate_sorted = all(left <= right for left, right in zip(original_keys, original_keys[1:]))
    candidate_features = sorted(
        (candidate_feature(candidate, config) for candidate in candidates),
        key=lambda feature: _feature_order_key(feature, seqid_rank),
    )

    if not coordinate_sorted:
        feature_started = False
        notice_written = False
        for line in lines:
            is_feature = bool(line.strip()) and not line.startswith("#")
            if is_feature:
                feature_started = True
                continue
            if not feature_started:
                handle.write(line)
                if line.startswith("##gff-version") and not notice_written:
                    _write_storf_notice(handle, config)
                    handle.write("# Original feature order normalized by sequence and coordinate.\n")
                    notice_written = True
        combined = [*result.document.features, *candidate_features]
        combined.sort(key=lambda feature: _feature_order_key(feature, seqid_rank))
        for feature in combined:
            handle.write(feature_line(feature) + "\n")
        return

    feature_number = 0
    candidate_number = 0
    notice_written = False
    for line in lines:
        is_feature = bool(line.strip()) and not line.startswith("#")
        if not is_feature:
            handle.write(line)
            if line.startswith("##gff-version") and not notice_written:
                _write_storf_notice(handle, config)
                notice_written = True
            continue
        original = result.document.features[feature_number]
        feature_number += 1
        original_position = _feature_order_key(original, seqid_rank)[:2]
        pending: list[Feature] = []
        while candidate_number < len(candidate_features):
            candidate = candidate_features[candidate_number]
            if _feature_order_key(candidate, seqid_rank)[:2] >= original_position:
                break
            pending.append(candidate)
            candidate_number += 1
        _write_candidate_feature_lines(handle, pending)
        handle.write(line)
    _write_candidate_feature_lines(handle, candidate_features[candidate_number:])


def write_augmented_gff(
    result: AnalysisResult,
    config: AnalysisConfig,
    path: str | Path,
    overwrite: bool = False,
    embed_fasta: bool | None = None,
    candidates: Sequence[StorfCandidate] | None = None,
) -> Path:
    path = Path(path)
    selected = result.candidates if candidates is None else candidates
    with atomic_output(path, overwrite) as handle:
        if result.document is not None:
            _write_interleaved_features(handle, result, config, selected)
        else:
            handle.write("##gff-version 3\n")
            for seqid, sequence in result.sequences.items():
                handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
            _write_storf_notice(handle, config)
            seqid_rank = {seqid: index for index, seqid in enumerate(result.sequences)}
            features = [*result.features, *(candidate_feature(item, config) for item in selected)]
            for feature in sorted(features, key=lambda item: _feature_order_key(item, seqid_rank)):
                handle.write(feature_line(feature) + "\n")
        should_embed = (
            result.document is not None and result.document.embedded_fasta is not None
            if embed_fasta is None
            else embed_fasta
        )
        if should_embed:
            handle.write("##FASTA\n")
            write_fasta_records(handle, result.sequences.items())
    return path


def write_candidate_gff(
    result: AnalysisResult,
    config: AnalysisConfig,
    path: str | Path,
    candidates: Sequence[StorfCandidate],
    overwrite: bool = False,
) -> Path:
    """Write a candidate-only GFF3 used for the lower-tier review partition."""

    path = Path(path)
    seqid_rank = {seqid: index for index, seqid in enumerate(result.sequences)}
    features = sorted(
        (candidate_feature(candidate, config) for candidate in candidates),
        key=lambda feature: _feature_order_key(feature, seqid_rank),
    )
    with atomic_output(path, overwrite) as handle:
        handle.write("##gff-version 3\n")
        for seqid, sequence in result.sequences.items():
            handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
        _write_storf_notice(handle, config)
        for feature in features:
            handle.write(feature_line(feature) + "\n")
    return path


def write_region_outputs(
    result: AnalysisResult,
    prefix: str | Path,
    gzip_output: bool = False,
    overwrite: bool = False,
    gene_ident: str | None = None,
) -> tuple[Path, Path]:
    prefix = Path(prefix)
    suffix = ".gz" if gzip_output else ""
    fasta_path = Path(f"{prefix}_UR.fasta{suffix}")
    gff_path = Path(f"{prefix}_UR.gff{suffix}")
    ensure_outputs_available((fasta_path, gff_path), overwrite)
    with atomic_output(fasta_path, overwrite) as handle:
        write_fasta_records(handle, ((region.id, region.sequence) for region in result.regions))
    with atomic_output(gff_path, overwrite) as handle:
        handle.write("##gff-version 3\n")
        if gene_ident is not None:
            comment_value = gene_ident.replace("\r", " ").replace("\n", " ")
            handle.write(f"# GeneIdent={comment_value}\n")
        for seqid, sequence in result.sequences.items():
            handle.write(f"##sequence-region {seqid} 1 {len(sequence)}\n")
        for region in result.regions:
            start, end = region.extended.to_gff()
            core_start, core_end = region.core.to_gff()
            attributes = (
                f"ID={region.id};OriginalUR={core_start}-{core_end};"
                f"LeftExtension={region.left_extension};RightExtension={region.right_extension}"
            )
            handle.write(
                f"{region.seqid}\tUR-Extractor\tunannotated_region\t{start}\t{end}\t.\t.\t.\t{attributes}\n"
            )
    return fasta_path, gff_path


def write_candidate_outputs(
    result: AnalysisResult,
    prefix: str | Path,
    gzip_output: bool = False,
    overwrite: bool = False,
    candidates: Sequence[StorfCandidate] | None = None,
    name_suffix: str = "",
    line_width: int = 60,
    include_leading_stop: bool = False,
    keep_terminal_stop: bool = False,
    protein_coordinate_mode: str = "candidate",
    remove_inframe_stops: bool = False,
) -> dict[str, Path]:
    prefix = Path(prefix)
    suffix = ".gz" if gzip_output else ""
    selected = result.candidates if candidates is None else candidates
    stem = f"{prefix}_StORFs{name_suffix}"
    paths = {
        "nucleotide_fasta": Path(f"{stem}.fasta{suffix}"),
        "protein_fasta": Path(f"{stem}_aa.fasta{suffix}"),
        "evidence_tsv": Path(f"{stem}.tsv{suffix}"),
    }
    ensure_outputs_available(paths.values(), overwrite)
    with atomic_output(paths["nucleotide_fasta"], overwrite) as handle:
        write_fasta_records(
            handle,
            (
                (
                    candidate.id,
                    candidate_nucleotide_sequence(
                        candidate,
                        include_leading_stop=include_leading_stop,
                    ),
                )
                for candidate in selected
            ),
            line_width,
        )
    with atomic_output(paths["protein_fasta"], overwrite) as handle:
        records = (
            (
                candidate.id,
                candidate_protein_sequence(
                    candidate,
                    coordinate_mode=protein_coordinate_mode,
                    include_leading_stop=include_leading_stop,
                    keep_terminal_stop=keep_terminal_stop,
                    remove_inframe_stops=remove_inframe_stops,
                ),
            )
            for candidate in selected
        )
        write_fasta_records(handle, records, line_width)
    with atomic_output(paths["evidence_tsv"], overwrite) as handle:
        columns = [
            "candidate_id",
            "seqid",
            "start",
            "end",
            "strand",
            "candidate_length",
            "coding_length",
            "region_id",
            "candidate_type",
            "genetic_code",
            "start_codon",
            "start_site_score",
            "start_candidate_count",
            "rbs_motif",
            "rbs_distance",
            "middle_stop_count",
            "middle_stop_codons",
            "core_fraction",
            "annotation_overlap_nt",
            "ambiguous_fraction",
            "evidence_score",
            "evidence_tier",
            "score_components",
        ]
        handle.write("\t".join(columns) + "\n")
        for candidate in selected:
            values = [
                candidate.id,
                candidate.seqid,
                str(candidate.interval.start + 1),
                str(candidate.interval.end),
                candidate.strand,
                str(candidate.length),
                "" if candidate.coding_length is None else str(candidate.coding_length),
                candidate.region_id,
                candidate.candidate_type,
                str(candidate.genetic_code),
                candidate.evidence.start_codon or "",
                (
                    ""
                    if candidate.evidence.start_site_score is None
                    else f"{candidate.evidence.start_site_score:.3f}"
                ),
                str(candidate.evidence.start_candidate_count),
                candidate.evidence.rbs_motif or "",
                "" if candidate.evidence.rbs_distance is None else str(candidate.evidence.rbs_distance),
                str(len(candidate.middle_stop_offsets)),
                ",".join(candidate.middle_stop_codons),
                f"{candidate.evidence.core_fraction:.6f}",
                str(candidate.evidence.overlap_annotated_nt),
                f"{candidate.evidence.ambiguous_fraction:.6f}",
                f"{candidate.evidence.score:.3f}",
                candidate.evidence.tier,
                json.dumps(candidate.evidence.components, sort_keys=True, separators=(",", ":")),
            ]
            handle.write("\t".join(values) + "\n")
    return paths


def write_combined_cds_storf_fasta(
    result: AnalysisResult,
    path: str | Path,
    amino_acid: bool = False,
    genetic_code: int = 11,
    overwrite: bool = False,
    line_width: int = 60,
    include_leading_stop: bool = False,
    protein_coordinate_mode: str = "candidate",
    remove_inframe_stops: bool = False,
) -> Path:
    """Write original CDS records followed by prioritised StORF candidates."""

    path = Path(path)
    records: list[tuple[str, str]] = []
    used: set[str] = set()
    cds_number = 0
    for feature in result.features:
        if feature.type != "CDS" or feature.seqid not in result.sequences:
            continue
        cds_number += 1
        identifier = feature.attributes.get("ID") or f"original_CDS_{cds_number}"
        identifier = str(identifier)
        if identifier in used:
            identifier = f"{identifier}_{cds_number}"
        used.add(identifier)
        sequence = result.sequences[feature.seqid][feature.interval.start : feature.interval.end]
        if feature.strand == "-":
            sequence = reverse_complement(sequence)
        if amino_acid:
            table = int(feature.attributes.get("transl_table", genetic_code))
            sequence = translate(
                sequence,
                table,
                keep_terminal_stop=False,
                start_as_methionine=True,
            )
        records.append((identifier, sequence))
    for candidate in result.candidates:
        if amino_acid:
            sequence = candidate_protein_sequence(
                candidate,
                coordinate_mode=protein_coordinate_mode,
                include_leading_stop=include_leading_stop,
                remove_inframe_stops=remove_inframe_stops,
            )
        else:
            sequence = candidate_nucleotide_sequence(
                candidate,
                include_leading_stop=include_leading_stop,
            )
        records.append((candidate.id, sequence))
    with atomic_output(path, overwrite) as handle:
        write_fasta_records(handle, records, line_width)
    return path


def _sha256(path: Path | None) -> str | None:
    if path is None:
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_manifest(
    result: AnalysisResult,
    config: AnalysisConfig,
    path: str | Path,
    outputs: Mapping[str, Path],
    overwrite: bool = False,
) -> Path:
    path = Path(path)
    payload = {
        "schema_version": 1,
        "tool": "StORF-Reporter",
        "version": StORF_Reporter_Version,
        "sample": result.sample,
        "inputs": {
            "gff": str(result.input_gff) if result.input_gff else None,
            "gff_sha256": _sha256(result.input_gff),
            "fasta": str(result.input_fasta) if result.input_fasta else None,
            "fasta_sha256": _sha256(result.input_fasta),
        },
        "configuration": asdict(config),
        "summary": asdict(result.summary),
        "outputs": {name: str(output) for name, output in outputs.items()},
        "scientific_notice": (
            "EvidenceScore is an interpretable prioritisation heuristic and is not "
            "a calibrated probability of a protein-coding gene."
        ),
    }
    with atomic_output(path, overwrite) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return path


def write_html_report(
    result: AnalysisResult,
    path: str | Path,
    overwrite: bool = False,
    primary_tier: str = "exploratory",
) -> Path:
    if primary_tier not in TIER_ORDER:
        raise ValueError("primary tier must be high, medium or exploratory")
    threshold = TIER_ORDER[primary_tier]
    path = Path(path)
    with atomic_output(path, overwrite) as handle:
        handle.write("<!doctype html><html lang='en'><head><meta charset='utf-8'>")
        handle.write(f"<title>StORF-Reporter – {html.escape(result.sample)}</title>")
        handle.write(
            "<style>body{font:15px system-ui;margin:2rem;color:#17212b}"
            "table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd4dc;"
            "padding:.4rem;text-align:left}th{background:#edf3f6}.high{background:#d8f3dc}"
            ".medium{background:#fff3bf}.exploratory{background:#f1f3f5}code{font-size:.85em}</style>"
        )
        handle.write("</head><body>")
        handle.write(f"<h1>StORF-Reporter: {html.escape(result.sample)}</h1>")
        handle.write(
            "<p><strong>Scientific notice:</strong> evidence scores prioritise review; "
            "they are not calibrated probabilities.</p>"
        )
        summary = result.summary
        handle.write(
            f"<p>{summary.sequences} sequences · {summary.annotations} annotations · "
            f"{summary.unannotated_regions} URs · {summary.candidates} retained candidates "
            f"({summary.filtered_overlap} overlap-filtered) · primary minimum tier: "
            f"{html.escape(primary_tier)}</p>"
        )
        handle.write("<table><thead><tr><th>ID</th><th>Locus</th><th>Candidate/Coding length</th>")
        handle.write(
            "<th>Start/RBS</th><th>Score</th><th>Tier</th><th>Output</th><th>UR</th></tr></thead><tbody>"
        )
        for candidate in result.candidates:
            evidence = candidate.evidence
            locus = (
                f"{candidate.seqid}:{candidate.interval.start + 1}-"
                f"{candidate.interval.end}({candidate.strand})"
            )
            handle.write(f"<tr class='{evidence.tier}'>")
            handle.write(f"<td><code>{html.escape(candidate.id)}</code></td>")
            coding_length = candidate.coding_length if candidate.coding_length is not None else "none"
            handle.write(f"<td>{html.escape(locus)}</td><td>{candidate.length} / {coding_length}</td>")
            handle.write(
                f"<td>{html.escape(evidence.start_codon or 'none')} / "
                f"{html.escape(evidence.rbs_motif or 'none')}</td>"
            )
            handle.write(f"<td>{evidence.score:.3f}</td><td>{evidence.tier}</td>")
            destination = "primary" if TIER_ORDER[evidence.tier] >= threshold else "below-tier review"
            handle.write(f"<td>{destination}</td>")
            handle.write(f"<td><code>{html.escape(candidate.region_id)}</code></td></tr>")
        handle.write("</tbody></table></body></html>\n")
    return path


def discover_pairs(
    input_path: str | Path,
    explicit_fasta: str | Path | None = None,
    recursive: bool = False,
) -> list[tuple[str, Path, Path | None]]:
    path = Path(input_path).resolve()
    if path.is_file():
        gffs = [path]
    elif path.is_dir():
        iterator = path.rglob("*") if recursive else path.glob("*")
        derived_markers = (
            "StORF-Reporter",
            "_UR.gff",
            "_UR_StORF-Finder",
            "_StORF-Finder",
            "_Extracted_StORFs",
            "_NoHit_StORFs",
            "_StORFs_below-",
        )
        gffs = sorted(
            item
            for item in iterator
            if item.is_file()
            and is_gff_path(item)
            and not any(marker in item.name for marker in derived_markers)
        )
    else:
        raise FileNotFoundError(f"input does not exist: {path}")
    if not gffs:
        raise FileNotFoundError(f"no GFF/GFF3 inputs found under {path}")
    if explicit_fasta is not None and len(gffs) != 1:
        raise ValueError("an explicit FASTA can only be paired with one GFF3 input")

    pairs: list[tuple[str, Path, Path | None]] = []
    for gff in gffs:
        document = read_gff(gff)
        fasta: Path | None = Path(explicit_fasta).resolve() if explicit_fasta else None
        if fasta is None and document.embedded_fasta is None:
            stem = logical_stem(gff)
            candidates = [
                gff.parent / f"{stem}{extension}"
                for extension in (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
            ]
            candidates.extend(
                sorted(item for item in gff.parent.iterdir() if item.is_file() and is_fasta_path(item))
            )
            fasta = next((item for item in candidates if item.exists()), None)
            unique_fastas = {item.resolve() for item in candidates if item.exists()}
            same_stem = [item for item in unique_fastas if logical_stem(item) == logical_stem(gff)]
            if same_stem:
                fasta = sorted(same_stem)[0]
            elif len(unique_fastas) != 1:
                raise InputFormatError(
                    f"could not unambiguously pair {gff} with a FASTA; use --fasta or a manifest"
                )
        pairs.append((logical_stem(gff), gff, fasta))
    return pairs


def read_sample_manifest(path: str | Path) -> list[tuple[str, Path, Path | None]]:
    manifest = Path(path).resolve()
    pairs: list[tuple[str, Path, Path | None]] = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or not {"sample", "gff"}.issubset(reader.fieldnames):
            raise InputFormatError("sample manifest requires sample and gff columns")
        for row_number, row in enumerate(reader, 2):
            sample = (row.get("sample") or "").strip()
            gff_text = (row.get("gff") or "").strip()
            fasta_text = (row.get("fasta") or "").strip()
            if not sample or not gff_text:
                raise InputFormatError(f"manifest row {row_number} is missing sample or gff")
            gff = (
                (manifest.parent / gff_text).resolve() if not Path(gff_text).is_absolute() else Path(gff_text)
            )
            fasta = None
            if fasta_text:
                fasta = (
                    (manifest.parent / fasta_text).resolve()
                    if not Path(fasta_text).is_absolute()
                    else Path(fasta_text)
                )
            pairs.append((sample, gff, fasta))
    return pairs
