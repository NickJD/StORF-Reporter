"""Typed records used by the StORF-Reporter v2 analysis core.

Coordinates are zero-based, half-open everywhere in memory.  Conversion to
one-based inclusive coordinates happens only when GFF3 is read or written.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field


@dataclass(frozen=True, order=True)
class Interval:
    """A zero-based, half-open interval ``[start, end)``."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0:
            raise ValueError(f"interval start must be >= 0, got {self.start}")
        if self.end < self.start:
            raise ValueError(f"interval end {self.end} is before start {self.start}")

    @property
    def length(self) -> int:
        return self.end - self.start

    def intersects(self, other: "Interval") -> bool:
        return self.start < other.end and other.start < self.end

    def intersection_length(self, other: "Interval") -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))

    def clip(self, lower: int, upper: int) -> "Interval":
        return Interval(max(lower, self.start), min(upper, self.end))

    def to_gff(self) -> tuple[int, int]:
        if self.length == 0:
            raise ValueError("empty intervals cannot be represented in GFF3")
        return self.start + 1, self.end

    @classmethod
    def from_gff(cls, start: int, end: int) -> "Interval":
        if start < 1 or end < start:
            raise ValueError(f"invalid GFF3 coordinates: {start}-{end}")
        return cls(start - 1, end)


@dataclass(frozen=True)
class Feature:
    seqid: str
    source: str
    type: str
    interval: Interval
    score: str = "."
    strand: str = "."
    phase: str = "."
    attributes: Mapping[str, str] = field(default_factory=dict)
    raw_line: str | None = None


@dataclass(frozen=True)
class UnannotatedRegion:
    id: str
    seqid: str
    core: Interval
    extended: Interval
    sequence: str

    @property
    def left_extension(self) -> int:
        return self.core.start - self.extended.start

    @property
    def right_extension(self) -> int:
        return self.extended.end - self.core.end


@dataclass(frozen=True)
class CandidateEvidence:
    start_codon: str | None
    start_offset: int | None
    rbs_motif: str | None
    rbs_distance: int | None
    start_site_score: float | None
    start_candidate_count: int
    ambiguous_fraction: float
    core_fraction: float
    overlap_annotated_nt: int
    score: float
    tier: str
    components: Mapping[str, float] = field(default_factory=dict)
    external: Mapping[str, float] = field(default_factory=dict)
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class StorfCandidate:
    id: str
    seqid: str
    interval: Interval
    strand: str
    frame: int
    sequence: str
    region_id: str
    candidate_type: str
    start_stop: str
    end_stop: str
    genetic_code: int
    evidence: CandidateEvidence
    middle_stop_offsets: tuple[int, ...] = ()

    @property
    def length(self) -> int:
        return self.interval.length

    @property
    def coding_interval(self) -> Interval | None:
        offset = self.evidence.start_offset
        if offset is None:
            return None
        if self.strand == "+":
            return Interval(self.interval.start + offset, self.interval.end)
        return Interval(self.interval.start, self.interval.end - offset)

    @property
    def coding_sequence(self) -> str | None:
        offset = self.evidence.start_offset
        return None if offset is None else self.sequence[offset:]

    @property
    def coding_length(self) -> int | None:
        sequence = self.coding_sequence
        return None if sequence is None else len(sequence)

    @property
    def middle_stop_codons(self) -> tuple[str, ...]:
        return tuple(self.sequence[offset : offset + 3] for offset in self.middle_stop_offsets)

    @property
    def start_before_first_middle_stop(self) -> bool:
        return bool(
            self.middle_stop_offsets
            and self.evidence.start_offset is not None
            and self.evidence.start_offset < self.middle_stop_offsets[0]
        )


@dataclass
class RunSummary:
    sequences: int = 0
    annotations: int = 0
    unannotated_regions: int = 0
    candidates_before_filtering: int = 0
    candidates: int = 0
    filtered_overlap: int = 0
    filtered_annotation_overlap: int = 0
    rejected_ambiguous: int = 0
    tiers: dict[str, int] = field(default_factory=lambda: {"high": 0, "medium": 0, "exploratory": 0})
    warnings: list[str] = field(default_factory=list)
