"""Structural validation used by ``storf augment`` and ``storf validate``.

Generic GFF3 validation deliberately does not require every pre-existing CDS
to have a start and stop codon: real annotations can contain partial and
multi-part CDS features.  The stricter translation checks are applied only to
CDS features produced by StORF-Reporter.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from .genetic_code import reverse_complement, start_codons, stop_codons
from .io import InputFormatError, read_fasta, read_gff
from .models import Interval, StorfCandidate


@dataclass(frozen=True)
class ValidationIssue:
    severity: str
    code: str
    message: str

    def as_dict(self) -> dict[str, str]:
        return asdict(self)


@dataclass(frozen=True)
class CodingValidation:
    """Result of validating one proposed complete coding sequence."""

    passed: bool
    issues: tuple[str, ...]


def validate_coding_sequence(
    sequence: str,
    genetic_code: int,
    minimum_length: int = 6,
) -> CodingValidation:
    """Check the structural requirements for StORF-to-CDS promotion.

    ``sequence`` must be oriented 5' to 3' and include the terminal stop.
    These checks establish a structurally complete ORF, not that it is a
    translated or functional gene.
    """

    sequence = sequence.upper()
    problems: list[str] = []
    if len(sequence) < minimum_length:
        problems.append(f"coding sequence is shorter than {minimum_length} nt")
    if len(sequence) % 3:
        problems.append("coding-sequence length is not divisible by three")
    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        problems.append(f"coding sequence contains ambiguous bases: {''.join(invalid)}")
    if len(sequence) >= 3 and sequence[:3] not in start_codons(genetic_code):
        problems.append(f"{sequence[:3]} is not a start codon under genetic code {genetic_code}")
    if len(sequence) >= 3 and sequence[-3:] not in stop_codons(genetic_code):
        problems.append(f"{sequence[-3:]} is not a stop codon under genetic code {genetic_code}")
    if len(sequence) >= 6 and not len(sequence) % 3:
        internal = [
            sequence[offset : offset + 3]
            for offset in range(3, len(sequence) - 3, 3)
            if sequence[offset : offset + 3] in stop_codons(genetic_code)
        ]
        if internal:
            problems.append(f"coding sequence contains {len(internal)} internal in-frame stop codon(s)")
    return CodingValidation(not problems, tuple(problems))


def validate_candidate_cds(
    candidate: StorfCandidate,
    minimum_length: int = 6,
) -> CodingValidation:
    """Validate whether a discovered candidate can be represented as CDS."""

    if candidate.coding_sequence is None:
        return CodingValidation(False, ("no supported in-frame start codon was found",))
    return validate_coding_sequence(candidate.coding_sequence, candidate.genetic_code, minimum_length)


def validate_con_storf_candidate(
    candidate: StorfCandidate,
    minimum_first_segment_length: int = 6,
) -> CodingValidation:
    """Validate the structure of a potential multi-stop readthrough candidate.

    A Con-StORF is not a conventional complete CDS.  Its expected middle stop
    codons are retained, and the selected start must occur before the first of
    them in translation orientation.  More than one middle stop is permitted.
    """

    sequence = candidate.sequence.upper()
    problems: list[str] = []
    if candidate.candidate_type != "Con-StORF":
        problems.append("candidate is not labelled Con-StORF")
    if len(sequence) % 3:
        problems.append("candidate length is not divisible by three")
    invalid = sorted(set(sequence) - set("ACGT"))
    if invalid:
        problems.append(f"candidate contains ambiguous bases: {''.join(invalid)}")
    stops = stop_codons(candidate.genetic_code)
    if len(sequence) >= 3 and sequence[:3] not in stops:
        problems.append("candidate does not begin at a supported delimiting stop")
    if len(sequence) >= 3 and sequence[-3:] not in stops:
        problems.append("candidate does not end at a supported delimiting stop")
    if not candidate.middle_stop_offsets:
        problems.append("Con-StORF has no middle in-frame stop codon")
    if candidate.evidence.start_offset is None:
        problems.append("no supported in-frame start codon was found before the first middle stop")
    elif candidate.middle_stop_offsets:
        first_middle = candidate.middle_stop_offsets[0]
        if candidate.evidence.start_offset >= first_middle:
            problems.append("selected start codon is not before the first middle stop")
        first_segment_length = first_middle + 3 - candidate.evidence.start_offset
        if first_segment_length < minimum_first_segment_length:
            problems.append(
                f"start-to-first-middle-stop segment is shorter than {minimum_first_segment_length} nt"
            )
    return CodingValidation(not problems, tuple(problems))


def _storf_cds_issues(feature, sequences) -> list[ValidationIssue]:
    """Apply complete-ORF checks to the ranked coding locus of a promoted feature."""

    if feature.seqid not in sequences:
        return []
    identifier = str(feature.attributes.get("ID", "unnamed StORF-Reporter CDS"))
    preliminary: list[ValidationIssue] = []
    if feature.strand not in {"+", "-"}:
        preliminary.append(
            ValidationIssue(
                "error",
                "storf_cds_strand",
                f"{identifier} must have '+' or '-' strand for sequence validation",
            )
        )
    if feature.phase != "0":
        preliminary.append(
            ValidationIssue(
                "error",
                "storf_cds_phase",
                f"{identifier} must have phase 0 for its ranked reading frame",
            )
        )
    if preliminary:
        return preliminary
    try:
        genetic_code = int(feature.attributes.get("GeneticCode", 11))
        # Resolve the code before extracting sequence so unsupported values get
        # one clear validation error rather than an uncaught exception.
        start_codons(genetic_code)
    except (TypeError, ValueError):
        return [
            ValidationIssue(
                "error",
                "storf_cds_genetic_code",
                f"{identifier} has an unsupported GeneticCode attribute",
            )
        ]

    if "CodingCoordinates" in feature.attributes:
        coding, coding_issue = _coordinate_attribute(feature, "CodingCoordinates")
        if coding_issue is not None:
            return [coding_issue]
        if coding is None:
            return [
                ValidationIssue(
                    "error",
                    "storf_cds_coordinates",
                    f"{identifier} CDS has no ranked CodingCoordinates",
                )
            ]
    else:
        # Pre-v2 StORF CDS records used the displayed interval as their only
        # sequence locus; retain validation support for those files.
        coding = feature.interval
    sequence = sequences[feature.seqid][coding.start : coding.end]
    if feature.strand == "-":
        sequence = reverse_complement(sequence)
    result = validate_coding_sequence(sequence, genetic_code)
    issues: list[ValidationIssue] = []
    code_by_prefix = {
        "coding sequence is shorter": "storf_cds_length",
        "coding-sequence length": "storf_cds_frame",
        "coding sequence contains ambiguous": "storf_cds_ambiguous",
        "coding sequence contains": "storf_cds_internal_stop",
    }
    for problem in result.issues:
        code = next((value for prefix, value in code_by_prefix.items() if problem.startswith(prefix)), None)
        if code is None:
            code = "storf_cds_start" if "start codon" in problem else "storf_cds_stop"
        issues.append(ValidationIssue("error", code, f"{identifier}: {problem}"))
    return issues


def _coordinate_attribute(feature, name: str) -> tuple[Interval | None, ValidationIssue | None]:
    text = str(feature.attributes.get(name, "")).strip()
    if not text or text.casefold() == "none":
        return None, None
    try:
        start_text, end_text = text.split("-", 1)
        return Interval.from_gff(int(start_text), int(end_text)), None
    except (TypeError, ValueError):
        identifier = str(feature.attributes.get("ID", "unnamed StORF feature"))
        return None, ValidationIssue(
            "error",
            "storf_coordinate_attribute",
            f"{identifier} has invalid {name}={text!r}; expected START-END",
        )


def _integer_attribute(feature, name: str) -> tuple[int | None, ValidationIssue | None]:
    text = str(feature.attributes.get(name, "")).strip()
    if not text or text.casefold() == "none":
        return None, None
    try:
        value = int(text)
    except ValueError:
        identifier = str(feature.attributes.get("ID", "unnamed StORF feature"))
        return None, ValidationIssue(
            "error", "storf_length_attribute", f"{identifier} has invalid {name}={text!r}"
        )
    return value, None


def _v2_storf_candidate_issues(feature, sequences) -> list[ValidationIssue]:
    """Validate the dual candidate/coding representation used by v2 outputs."""

    if "CandidateCoordinates" not in feature.attributes:
        return []
    identifier = str(feature.attributes.get("ID", "unnamed StORF feature"))
    issues: list[ValidationIssue] = []
    candidate, candidate_issue = _coordinate_attribute(feature, "CandidateCoordinates")
    coding, coding_issue = _coordinate_attribute(feature, "CodingCoordinates")
    candidate_length, candidate_length_issue = _integer_attribute(feature, "CandidateLength")
    coding_length, coding_length_issue = _integer_attribute(feature, "CodingLength")
    for issue in (candidate_issue, coding_issue, candidate_length_issue, coding_length_issue):
        if issue is not None:
            issues.append(issue)
    if candidate is None:
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_coordinates",
                f"{identifier} has no CandidateCoordinates",
            )
        )
        return issues
    if candidate_length != candidate.length:
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_length",
                f"{identifier} CandidateLength does not match CandidateCoordinates",
            )
        )
    if coding is not None:
        if candidate.intersection_length(coding) != coding.length:
            issues.append(
                ValidationIssue(
                    "error",
                    "storf_coding_containment",
                    f"{identifier} CodingCoordinates are outside CandidateCoordinates",
                )
            )
        if coding_length != coding.length:
            issues.append(
                ValidationIssue(
                    "error",
                    "storf_coding_length",
                    f"{identifier} CodingLength does not match CodingCoordinates",
                )
            )
    elif coding_length is not None:
        issues.append(
            ValidationIssue(
                "error",
                "storf_coding_length",
                f"{identifier} has CodingLength but no CodingCoordinates",
            )
        )

    coordinate_mode = str(feature.attributes.get("GFFCoordinateMode", ""))
    expected_display = candidate
    if coordinate_mode in {"coding", "start-to-stop"} and coding is not None:
        expected_display = coding
    elif coordinate_mode == "stop-to-stop-first-stop-excluded":
        expected_display = (
            Interval(candidate.start + 3, candidate.end)
            if feature.strand == "+"
            else Interval(candidate.start, candidate.end - 3)
        )
    if (
        coordinate_mode
        in {
            "stop-to-stop",
            "stop-to-stop-inclusive",
            "stop-to-stop-first-stop-excluded",
            "coding",
            "start-to-stop",
        }
        and feature.interval != expected_display
    ):
        issues.append(
            ValidationIssue(
                "error",
                "storf_display_coordinates",
                f"{identifier} GFF columns disagree with GFFCoordinateMode={coordinate_mode}",
            )
        )
    if feature.type == "CDS" and coding is None:
        issues.append(
            ValidationIssue(
                "error",
                "storf_cds_coordinates",
                f"{identifier} CDS has no ranked CodingCoordinates",
            )
        )

    if feature.seqid not in sequences or candidate.end > len(sequences[feature.seqid]):
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_bounds",
                f"{identifier} CandidateCoordinates exceed the matching sequence",
            )
        )
        return issues
    if feature.strand not in {"+", "-"}:
        issues.append(
            ValidationIssue("error", "storf_candidate_strand", f"{identifier} must have '+' or '-' strand")
        )
        return issues
    try:
        genetic_code = int(feature.attributes.get("GeneticCode", 11))
        stops = stop_codons(genetic_code)
    except (TypeError, ValueError):
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_genetic_code",
                f"{identifier} has an unsupported GeneticCode attribute",
            )
        )
        return issues

    candidate_sequence = sequences[feature.seqid][candidate.start : candidate.end]
    if feature.strand == "-":
        candidate_sequence = reverse_complement(candidate_sequence)
    if candidate.length % 3:
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_frame",
                f"{identifier} CandidateLength is not divisible by three",
            )
        )
    if candidate_sequence[:3] not in stops or candidate_sequence[-3:] not in stops:
        issues.append(
            ValidationIssue(
                "error",
                "storf_candidate_stops",
                f"{identifier} CandidateCoordinates do not delimit in-frame stop codons",
            )
        )
    internal_stops = [
        candidate_sequence[offset : offset + 3]
        for offset in range(3, len(candidate_sequence) - 3, 3)
        if candidate_sequence[offset : offset + 3] in stops
    ]
    candidate_type = str(feature.attributes.get("StORF_Type", "StORF"))
    middle_count, middle_count_issue = _integer_attribute(feature, "MiddleStopCount")
    if middle_count_issue is not None:
        issues.append(middle_count_issue)
    if candidate_type == "Con-StORF":
        if middle_count != len(internal_stops) or not internal_stops:
            issues.append(
                ValidationIssue(
                    "error",
                    "con_storf_middle_stops",
                    f"{identifier} MiddleStopCount does not match its internal in-frame stops",
                )
            )
    elif internal_stops:
        issues.append(
            ValidationIssue(
                "error",
                "storf_internal_stop",
                f"{identifier} ordinary StORF contains an internal in-frame stop",
            )
        )

    if coding is not None and coding.end <= len(sequences[feature.seqid]):
        coding_sequence = sequences[feature.seqid][coding.start : coding.end]
        if feature.strand == "-":
            coding_sequence = reverse_complement(coding_sequence)
        declared_start = str(feature.attributes.get("StartCodon", "none"))
        if declared_start != coding_sequence[:3]:
            issues.append(
                ValidationIssue(
                    "error",
                    "storf_start_coordinate",
                    f"{identifier} StartCodon does not match CodingCoordinates",
                )
            )
        structural = validate_coding_sequence(coding_sequence, genetic_code)
        if feature.attributes.get("CDSValidation") == "passed" and not structural.passed:
            issues.append(
                ValidationIssue(
                    "error",
                    "storf_saved_cds_validation",
                    f"{identifier} is marked CDSValidation=passed but its coding sequence fails",
                )
            )
    return issues


def validate_inputs(
    gff_path: str | Path,
    fasta_path: str | Path | None = None,
) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    try:
        document = read_gff(gff_path)
    except (OSError, InputFormatError) as error:
        return [ValidationIssue("error", "gff_parse", str(error))]

    try:
        sequences = read_fasta(fasta_path) if fasta_path else document.embedded_fasta
    except (OSError, InputFormatError) as error:
        return [ValidationIssue("error", "fasta_parse", str(error))]
    if sequences is None:
        issues.append(
            ValidationIssue("error", "missing_fasta", "GFF3 has no embedded FASTA and no FASTA was supplied")
        )
        return issues

    first_nonempty = next((line.strip() for line in document.pre_fasta_lines if line.strip()), "")
    if first_nonempty != "##gff-version 3":
        issues.append(
            ValidationIssue(
                "warning",
                "gff_version",
                "first non-empty line should be exactly '##gff-version 3'",
            )
        )

    seen_ids: set[str] = set()
    for feature in document.features:
        if feature.seqid not in sequences:
            issues.append(
                ValidationIssue(
                    "error", "missing_seqid", f"{feature.seqid!r} is present in GFF3 but absent from FASTA"
                )
            )
            continue
        if feature.interval.end > len(sequences[feature.seqid]):
            issues.append(
                ValidationIssue(
                    "error",
                    "coordinate_bounds",
                    (
                        f"{feature.seqid}:{feature.interval.start + 1}-"
                        f"{feature.interval.end} exceeds sequence length"
                    ),
                )
            )
        if feature.type == "CDS" and feature.phase not in {"0", "1", "2"}:
            issues.append(
                ValidationIssue(
                    "error",
                    "cds_phase",
                    (
                        f"CDS {feature.seqid}:{feature.interval.start + 1}-"
                        f"{feature.interval.end} has phase {feature.phase!r}"
                    ),
                )
            )
        if (
            feature.type == "CDS"
            and feature.interval.end <= len(sequences[feature.seqid])
            and (
                feature.source == "StORF-Reporter"
                or str(feature.attributes.get("ID", "")).startswith("StORF_")
            )
        ):
            issues.extend(_storf_cds_issues(feature, sequences))
        if feature.source == "StORF-Reporter" or str(feature.attributes.get("ID", "")).startswith("StORF_"):
            issues.extend(_v2_storf_candidate_issues(feature, sequences))
        feature_id = feature.attributes.get("ID")
        if feature_id:
            if feature_id in seen_ids:
                issues.append(ValidationIssue("error", "duplicate_id", f"duplicate GFF3 ID {feature_id!r}"))
            seen_ids.add(feature_id)
        elif feature.source == "StORF-Reporter":
            issues.append(
                ValidationIssue("error", "missing_id", "StORF-Reporter feature has no ID attribute")
            )

    for seqid, sequence in sequences.items():
        invalid = sorted(set(sequence.upper()) - set("ACGTRYSWKMBDHVN.-"))
        if invalid:
            issues.append(
                ValidationIssue(
                    "warning",
                    "non_iupac",
                    f"sequence {seqid!r} contains non-IUPAC characters: {''.join(invalid)}",
                )
            )
    return issues
