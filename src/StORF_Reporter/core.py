"""Correctness-first StORF discovery and evidence scoring.

The functions in this module are deliberately independent of command-line
parsing and file naming so they can be exercised with tiny synthetic tests.
"""

from __future__ import annotations

import hashlib
import re
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace

from .genetic_code import reverse_complement, start_codons, stop_codons
from .models import (
    CandidateEvidence,
    Feature,
    Interval,
    StorfCandidate,
    UnannotatedRegion,
)

DEFAULT_EXCLUDED_FEATURES = frozenset(
    {
        "CDS",
        "gene",
        "mRNA",
        "misc_RNA",
        "ncRNA",
        "oriC",
        "pseudo",
        "pseudogene",
        "regulatory_region",
        "repeat_region",
        "rRNA",
        "tmRNA",
        "tRNA",
        "CRISPR",
    }
)

# Historical ``-gene_ident Prokka`` policy.  Keep this list explicit rather
# than deriving it from DEFAULT_EXCLUDED_FEATURES: the legacy public option is
# a reproducibility contract and intentionally has its original semantics.
PROKKA_GENE_IDENT_TYPES = (
    "misc_RNA",
    "gene",
    "mRNA",
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "CRISPR",
    "ncRNA",
    "regulatory_region",
    "oriC",
    "pseudo",
)

RBS_MOTIFS = ("AGGAGG", "GGAGG", "AGGAG", "GAGG", "GGAG")
TIER_ORDER = {"exploratory": 0, "medium": 1, "high": 2}

_START_CODON_PRIOR = {"ATG": 20.0, "GTG": 15.0, "TTG": 12.0, "CTG": 8.0}
_RBS_MOTIF_PRIOR = {"AGGAGG": 12.0, "GGAGG": 10.0, "AGGAG": 10.0, "GAGG": 8.0, "GGAG": 8.0}


@dataclass(frozen=True)
class _StartSelection:
    offset: int | None
    codon: str | None
    rbs_motif: str | None
    rbs_distance: int | None
    score: float | None
    candidate_count: int


def resolve_gene_ident(gene_ident: str | Iterable[str]) -> tuple[str, tuple[str, ...]]:
    """Resolve the historical gene-ident expression into a selector.

    The returned mode is either ``"type"`` for exact GFF3 feature-type
    matching or ``"attribute"`` for a single ``KEY=PREFIX`` attribute query.
    ``ID=gene`` therefore matches Ensembl identifiers such as
    ``ID=gene:ER3413_1`` without pretending that every occupied record has the
    exact feature type ``gene``.
    """

    if isinstance(gene_ident, str):
        text = gene_ident.strip()
    else:
        text = ",".join(str(item).strip() for item in gene_ident if str(item).strip())
    if not text:
        raise ValueError("gene-ident must specify at least one feature type or attribute query")
    if text.casefold() in {"prokka", "bakta"}:
        return "type", PROKKA_GENE_IDENT_TYPES

    identifiers = tuple(item.strip() for item in text.split(",") if item.strip())
    if len(identifiers) == 1 and "=" in identifiers[0]:
        key, prefix = (item.strip() for item in identifiers[0].split("=", 1))
        if not key or not prefix:
            raise ValueError("attribute gene-ident must use KEY=PREFIX, for example ID=gene")
        return "attribute", (key, prefix)
    if any("=" in identifier for identifier in identifiers):
        raise ValueError(
            "attribute gene-ident cannot be mixed with feature types; use one query such as ID=gene"
        )
    return "type", identifiers


def feature_matches_gene_ident(feature: Feature, gene_ident: str | Iterable[str]) -> bool:
    """Return whether ``feature`` occupies sequence under ``gene_ident``."""

    mode, values = resolve_gene_ident(gene_ident)
    if mode == "type":
        return feature.type in frozenset(values)
    key, prefix = values
    value = feature.attributes.get(key)
    return value is not None and str(value).startswith(prefix)


def merge_intervals(intervals: Iterable[Interval]) -> list[Interval]:
    ordered = sorted(intervals)
    if not ordered:
        return []
    merged = [ordered[0]]
    for interval in ordered[1:]:
        previous = merged[-1]
        if interval.start <= previous.end:
            merged[-1] = Interval(previous.start, max(previous.end, interval.end))
        else:
            merged.append(interval)
    return merged


def extract_unannotated_regions(
    sequences: Mapping[str, str],
    features: Sequence[Feature],
    excluded_types: Iterable[str] = ("CDS",),
    min_length: int = 30,
    max_length: int = 100_000,
    extension: int = 50,
    gene_ident: str | Iterable[str] | None = None,
) -> list[UnannotatedRegion]:
    """Return complements of selected annotation intervals.

    UR length thresholds apply to the unextended core. Extensions are clipped
    to sequence bounds and retained separately so overlap with annotations is
    measurable rather than hidden.
    """

    if min_length < 1:
        raise ValueError("minimum UR length must be at least 1")
    if max_length < min_length:
        raise ValueError("maximum UR length must be >= minimum UR length")
    if extension < 0:
        raise ValueError("UR extension must be >= 0")

    excluded = frozenset(excluded_types)
    selector = resolve_gene_ident(gene_ident) if gene_ident is not None else None
    selector_types = frozenset(selector[1]) if selector is not None and selector[0] == "type" else frozenset()
    by_seqid: dict[str, list[Interval]] = defaultdict(list)
    for feature in features:
        if feature.seqid not in sequences:
            raise ValueError(f"GFF3 sequence {feature.seqid!r} is absent from the FASTA input")
        seq_length = len(sequences[feature.seqid])
        if feature.interval.end > seq_length:
            raise ValueError(
                f"feature {feature.seqid}:{feature.interval.start + 1}-"
                f"{feature.interval.end} exceeds sequence length {seq_length}"
            )
        if selector is None:
            occupies_sequence = feature.type in excluded
        elif selector[0] == "type":
            occupies_sequence = feature.type in selector_types
        else:
            key, prefix = selector[1]
            value = feature.attributes.get(key)
            occupies_sequence = value is not None and str(value).startswith(prefix)
        if occupies_sequence:
            by_seqid[feature.seqid].append(feature.interval)

    regions: list[UnannotatedRegion] = []
    for seqid, raw_sequence in sequences.items():
        sequence = raw_sequence.upper()
        seq_length = len(sequence)
        annotations = merge_intervals(by_seqid.get(seqid, ()))
        gaps: list[Interval] = []
        cursor = 0
        for annotation in annotations:
            if annotation.start > cursor:
                gaps.append(Interval(cursor, annotation.start))
            cursor = max(cursor, annotation.end)
        if cursor < seq_length:
            gaps.append(Interval(cursor, seq_length))

        for core in gaps:
            if not min_length <= core.length <= max_length:
                continue
            extended = Interval(
                max(0, core.start - extension),
                min(seq_length, core.end + extension),
            )
            gff_start, gff_end = core.to_gff()
            region_id = f"{seqid}_UR_{gff_start}_{gff_end}"
            regions.append(
                UnannotatedRegion(
                    id=region_id,
                    seqid=seqid,
                    core=core,
                    extended=extended,
                    sequence=sequence[extended.start : extended.end],
                )
            )
    return regions


def _stable_candidate_id(
    seqid: str,
    interval: Interval,
    strand: str,
    code: int,
    candidate_type: str,
) -> str:
    identity = f"{seqid}|{interval.start}|{interval.end}|{strand}|{code}|{candidate_type}"
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:14]
    return f"StORF_{digest}"


def _rbs(
    sequence: str,
    start_offset: int,
    upstream_context: str = "",
) -> tuple[str | None, int | None, float]:
    """Return the strongest Shine-Dalgarno-like motif for one proposed start.

    ``upstream_context`` permits motif detection immediately before the
    delimiting stop codon rather than artificially truncating the search at
    the stop-to-stop candidate boundary.
    """

    prefix = upstream_context[-20:]
    contextual = prefix + sequence
    contextual_start = len(prefix) + start_offset
    window_start = max(0, contextual_start - 20)
    window_end = max(window_start, contextual_start - 3)
    upstream = contextual[window_start:window_end]
    best: tuple[str, int, float] | None = None
    for motif in RBS_MOTIFS:
        for match in re.finditer(motif, upstream):
            distance = contextual_start - (window_start + match.end())
            spacing = max(0.0, 6.0 - abs(distance - 7.0))
            quality = _RBS_MOTIF_PRIOR[motif] + spacing
            proposal = (motif, distance, quality)
            if best is None or (quality, -abs(distance - 7), len(motif)) > (
                best[2],
                -abs(best[1] - 7),
                len(best[0]),
            ):
                best = proposal
    return best if best is not None else (None, None, 0.0)


def _best_start(
    sequence: str,
    code: int,
    upstream_context: str = "",
    search_end: int | None = None,
    minimum_coding_length: int = 6,
) -> _StartSelection:
    """Rank every permitted in-frame start using local sequence evidence.

    RBS support, RBS spacing, start-codon prior and retained coding length are
    considered.  Earlier coordinates are used only as the final deterministic
    tie-break.  For Con-StORFs, ``search_end`` is the first middle stop so a
    downstream ordinary-ORF start cannot masquerade as a Con-StORF start.
    """

    allowed = set(start_codons(code))
    upper = min(len(sequence) - 3, search_end if search_end is not None else len(sequence) - 3)
    proposals: list[tuple[tuple[float, int, float, float, int, int], _StartSelection]] = []
    scoring_end = search_end + 3 if search_end is not None else len(sequence)
    maximum_coding_length = max(1, scoring_end - 3)
    for offset in range(3, upper, 3):
        codon = sequence[offset : offset + 3]
        if codon not in allowed:
            continue
        coding_length = scoring_end - offset
        if coding_length < minimum_coding_length:
            continue
        motif, distance, rbs_quality = _rbs(sequence, offset, upstream_context)
        codon_prior = _START_CODON_PRIOR.get(codon, 6.0)
        coverage = 40.0 * coding_length / maximum_coding_length
        rbs_support = 15.0 + min(5.0, rbs_quality / 4.0) if motif else 0.0
        score = codon_prior + coverage + rbs_support
        rank = (
            score,
            int(motif is not None),
            rbs_quality,
            codon_prior,
            coding_length,
            -offset,
        )
        proposals.append(
            (
                rank,
                _StartSelection(offset, codon, motif, distance, round(score, 3), 0),
            )
        )
    if not proposals:
        return _StartSelection(None, None, None, None, None, 0)
    selected = max(proposals, key=lambda item: item[0])[1]
    return replace(selected, candidate_count=len(proposals))


def _score_evidence(
    sequence: str,
    start_offset: int | None,
    start_codon: str | None,
    rbs_motif: str | None,
    ambiguous_fraction: float,
    core_fraction: float,
    external: Mapping[str, float] | None,
) -> tuple[float, str, dict[str, float], tuple[str, ...]]:
    coding_length = len(sequence) - (start_offset or 3)
    components: dict[str, float] = {
        "length": min(25.0, max(0.0, (coding_length - 30) / 12.0)),
        "core_coverage": round(20.0 * core_fraction, 3),
        "ambiguity": round(-50.0 * ambiguous_fraction, 3),
    }
    if start_codon is not None:
        components["start_codon"] = {"ATG": 20.0, "GTG": 15.0, "TTG": 12.0}.get(start_codon, 8.0)
    else:
        components["start_codon"] = 0.0
    components["rbs"] = 15.0 if rbs_motif else 0.0

    clean_external: dict[str, float] = {}
    for name, value in (external or {}).items():
        numeric = max(0.0, min(1.0, float(value)))
        clean_external[name] = numeric
        components[f"external_{name}"] = 20.0 * numeric

    score = round(max(0.0, min(100.0, sum(components.values()))), 3)
    tier = "high" if score >= 65 else "medium" if score >= 40 else "exploratory"
    notes: list[str] = ["EvidenceScore is a transparent prioritisation heuristic, not a probability."]
    if start_codon is None:
        notes.append("No supported in-frame start codon was found.")
    if core_fraction < 1.0:
        notes.append("Candidate overlaps the annotation-extension flank.")
    if ambiguous_fraction:
        notes.append("Candidate contains non-ACGT nucleotides.")
    return score, tier, components, tuple(notes)


def _build_candidate(
    region: UnannotatedRegion,
    oriented_sequence: str,
    local_start: int,
    local_end: int,
    strand: str,
    oriented_frame: int,
    candidate_type: str,
    genetic_code: int,
    minimum_coding_length: int,
    external_evidence: Mapping[str, Mapping[str, float]] | None,
) -> StorfCandidate:
    sequence = oriented_sequence[local_start:local_end]
    middle_stop_offsets = tuple(
        offset
        for offset in range(3, len(sequence) - 3, 3)
        if sequence[offset : offset + 3] in stop_codons(genetic_code)
    )
    if strand == "+":
        interval = Interval(
            region.extended.start + local_start,
            region.extended.start + local_end,
        )
        frame = oriented_frame + 1
    else:
        region_length = len(oriented_sequence)
        interval = Interval(
            region.extended.start + region_length - local_end,
            region.extended.start + region_length - local_start,
        )
        frame = oriented_frame + 4

    candidate_id = _stable_candidate_id(region.seqid, interval, strand, genetic_code, candidate_type)
    upstream_context = oriented_sequence[max(0, local_start - 20) : local_start]
    search_end = middle_stop_offsets[0] if candidate_type == "Con-StORF" else None
    start = _best_start(
        sequence,
        genetic_code,
        upstream_context,
        search_end,
        minimum_coding_length,
    )
    ambiguous = sum(base not in "ACGT" for base in sequence) / len(sequence)
    core_overlap = interval.intersection_length(region.core)
    core_fraction = core_overlap / interval.length
    external = (external_evidence or {}).get(candidate_id, {})
    score, tier, components, notes = _score_evidence(
        sequence,
        start.offset,
        start.codon,
        start.rbs_motif,
        ambiguous,
        core_fraction,
        external,
    )
    evidence = CandidateEvidence(
        start_codon=start.codon,
        start_offset=start.offset,
        rbs_motif=start.rbs_motif,
        rbs_distance=start.rbs_distance,
        start_site_score=start.score,
        start_candidate_count=start.candidate_count,
        ambiguous_fraction=ambiguous,
        core_fraction=core_fraction,
        overlap_annotated_nt=interval.length - core_overlap,
        score=score,
        tier=tier,
        components=components,
        external=dict(external),
        notes=notes,
    )
    return StorfCandidate(
        id=candidate_id,
        seqid=region.seqid,
        interval=interval,
        strand=strand,
        frame=frame,
        sequence=sequence,
        region_id=region.id,
        candidate_type=candidate_type,
        start_stop=sequence[:3],
        end_stop=sequence[-3:],
        genetic_code=genetic_code,
        evidence=evidence,
        middle_stop_offsets=middle_stop_offsets,
    )


def find_storf_candidates(
    region: UnannotatedRegion,
    genetic_code: int = 11,
    min_orf: int = 99,
    max_orf: int = 60_000,
    ambiguous_threshold: float = 0.20,
    include_short: bool = False,
    short_min: int = 30,
    short_max: int = 120,
    require_start: bool = False,
    include_consecutive: bool = False,
    consecutive_only: bool = False,
    external_evidence: Mapping[str, Mapping[str, float]] | None = None,
) -> tuple[list[StorfCandidate], int]:
    if min_orf < 6 or min_orf % 3:
        raise ValueError("minimum StORF length must be a multiple of 3 and at least 6")
    if max_orf < min_orf:
        raise ValueError("maximum StORF length must be >= minimum StORF length")
    if not 0 <= ambiguous_threshold <= 1:
        raise ValueError("ambiguous threshold must be between 0 and 1")
    if include_short and (short_min < 6 or short_max < short_min):
        raise ValueError("invalid short-StORF length range")

    stops = stop_codons(genetic_code)
    candidates: list[StorfCandidate] = []
    rejected_ambiguous = 0
    for strand, oriented in (("+", region.sequence.upper()), ("-", reverse_complement(region.sequence))):
        for frame in range(3):
            positions = [
                offset
                for offset in range(frame, len(oriented) - 2, 3)
                if oriented[offset : offset + 3] in stops
            ]
            candidate_spans: list[tuple[int, int, str]] = []
            adjacent_lengths = [second + 3 - first for first, second in zip(positions, positions[1:])]
            if not consecutive_only:
                for first, second, length in zip(positions, positions[1:], adjacent_lengths):
                    if min_orf <= length <= max_orf:
                        candidate_spans.append((first, second, "StORF"))
                    elif include_short and short_min <= length <= short_max:
                        candidate_spans.append((first, second, "Short-StORF"))
            if include_consecutive or consecutive_only:
                # A Con-StORF is the maximal run of at least two adjacent,
                # qualifying ordinary StORFs.  Such a run can contain one or
                # many middle stops, matching the original public definition.
                run_start: int | None = None
                for index in range(len(adjacent_lengths) + 1):
                    qualifying = (
                        index < len(adjacent_lengths) and min_orf <= adjacent_lengths[index] <= max_orf
                    )
                    if qualifying and run_start is None:
                        run_start = index
                    if not qualifying and run_start is not None:
                        edge_count = index - run_start
                        if edge_count >= 2:
                            first = positions[run_start]
                            second = positions[index]
                            if second + 3 - first <= max_orf:
                                candidate_spans.append((first, second, "Con-StORF"))
                        run_start = None
            for first, second, candidate_type in candidate_spans:
                local_start, local_end = first, second + 3
                candidate = _build_candidate(
                    region,
                    oriented,
                    local_start,
                    local_end,
                    strand,
                    frame,
                    candidate_type,
                    genetic_code,
                    max(
                        6,
                        (short_min if candidate_type == "Short-StORF" else min_orf) - 3,
                    ),
                    external_evidence,
                )
                # Flanks help find boundary-crossing candidates, but a candidate
                # entirely inside an extension is not an unannotated-region hit.
                if not candidate.interval.intersects(region.core):
                    continue
                if candidate.evidence.ambiguous_fraction > ambiguous_threshold:
                    rejected_ambiguous += 1
                    continue
                if require_start and candidate.evidence.start_codon is None:
                    continue
                candidates.append(candidate)
    candidates.sort(key=lambda item: (item.seqid, item.interval.start, item.interval.end, item.strand))
    return candidates, rejected_ambiguous


def apply_overlap_filter(
    candidates: Sequence[StorfCandidate],
    maximum_overlap: int = 50,
    mode: str = "all",
    coordinate_basis: str = "coding",
) -> tuple[list[StorfCandidate], int]:
    """Greedily retain the strongest candidate in conflicts.

    Ties are deterministic: evidence score, then length, genomic coordinate,
    strand and stable identifier. Conflicts use ranked coding intervals when
    available unless ``coordinate_basis="candidate"`` is requested.
    """

    if maximum_overlap < 0:
        raise ValueError("maximum overlap must be >= 0")
    if mode not in {"none", "same-strand", "all"}:
        raise ValueError("overlap mode must be none, same-strand or all")
    if coordinate_basis not in {"coding", "candidate"}:
        raise ValueError("overlap coordinate basis must be coding or candidate")
    if mode == "none":
        return list(candidates), 0

    ranked = sorted(
        candidates,
        key=lambda item: (
            -item.evidence.score,
            -item.length,
            item.seqid,
            item.interval.start,
            item.interval.end,
            item.strand,
            item.id,
        ),
    )
    kept: list[StorfCandidate] = []
    removed = 0
    by_seqid: dict[str, list[StorfCandidate]] = defaultdict(list)
    for candidate in ranked:
        conflicts = False
        candidate_interval = (
            candidate.coding_interval
            if coordinate_basis == "coding" and candidate.coding_interval is not None
            else candidate.interval
        )
        for existing in by_seqid[candidate.seqid]:
            if mode == "same-strand" and candidate.strand != existing.strand:
                continue
            existing_interval = (
                existing.coding_interval
                if coordinate_basis == "coding" and existing.coding_interval is not None
                else existing.interval
            )
            if candidate_interval.intersection_length(existing_interval) > maximum_overlap:
                conflicts = True
                break
        if conflicts:
            removed += 1
        else:
            kept.append(candidate)
            by_seqid[candidate.seqid].append(candidate)
    kept.sort(key=lambda item: (item.seqid, item.interval.start, item.interval.end, item.strand))
    return kept, removed
