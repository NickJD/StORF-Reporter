from __future__ import annotations

from dataclasses import replace

import pytest

from StORF_Reporter.core import (
    apply_overlap_filter,
    extract_unannotated_regions,
    feature_matches_gene_ident,
    find_storf_candidates,
    resolve_gene_ident,
)
from StORF_Reporter.genetic_code import reverse_complement, translate
from StORF_Reporter.models import Feature, Interval, UnannotatedRegion
from StORF_Reporter.workflow import AnalysisConfig, candidate_feature


def feature(
    start: int,
    end: int,
    kind: str = "CDS",
    seqid: str = "ctg",
    attributes: dict[str, str] | None = None,
) -> Feature:
    return Feature(
        seqid,
        "test",
        kind,
        Interval.from_gff(start, end),
        phase="0" if kind == "CDS" else ".",
        attributes=attributes or {},
    )


def region(sequence: str, core: Interval | None = None) -> UnannotatedRegion:
    interval = Interval(0, len(sequence))
    return UnannotatedRegion("ctg_UR_1_x", "ctg", core or interval, interval, sequence)


def storf_sequence(coding_codons: int = 31, start: str = "ATG") -> str:
    return "TAA" + start + ("AAA" * (coding_codons - 1)) + "TAG"


def test_interval_round_trip() -> None:
    interval = Interval.from_gff(1, 100)
    assert interval == Interval(0, 100)
    assert interval.to_gff() == (1, 100)


def test_no_annotations_returns_whole_contig() -> None:
    regions = extract_unannotated_regions({"ctg": "A" * 100}, [], min_length=1)
    assert [(item.core.start, item.core.end) for item in regions] == [(0, 100)]


def test_whole_contig_annotation_returns_no_ur() -> None:
    regions = extract_unannotated_regions({"ctg": "A" * 100}, [feature(1, 100)], min_length=1)
    assert regions == []


def test_unsorted_overlapping_features_are_merged_before_complement() -> None:
    features = [feature(701, 800), feature(101, 250), feature(201, 300)]
    regions = extract_unannotated_regions({"ctg": "A" * 1000}, features, min_length=1, extension=0)
    assert [item.core.to_gff() for item in regions] == [(1, 100), (301, 700), (801, 1000)]


def test_terminal_ur_uses_last_feature_end() -> None:
    regions = extract_unannotated_regions({"ctg": "A" * 1000}, [feature(101, 200)], min_length=1, extension=0)
    assert [item.core.to_gff() for item in regions] == [(1, 100), (201, 1000)]


def test_extension_is_bounded_and_core_is_retained() -> None:
    sequence = "".join(str(index % 10) for index in range(100))
    regions = extract_unannotated_regions({"ctg": sequence}, [feature(21, 80)], min_length=1, extension=10)
    assert regions[0].core == Interval(0, 20)
    assert regions[0].extended == Interval(0, 30)
    assert regions[0].sequence == sequence[:30]
    assert regions[1].core == Interval(80, 100)
    assert regions[1].extended == Interval(70, 100)


def test_gene_ident_attribute_query_matches_attribute_prefix_across_types() -> None:
    features = [
        feature(11, 20, "gene", attributes={"ID": "gene:protein"}),
        feature(31, 40, "ncRNA_gene", attributes={"ID": "gene:ncrna"}),
        feature(51, 60, "pseudogene", attributes={"ID": "gene:pseudo"}),
        feature(71, 80, "gene", attributes={"ID": "not-a-gene-prefix"}),
    ]
    regions = extract_unannotated_regions(
        {"ctg": "A" * 100},
        features,
        min_length=1,
        extension=0,
        gene_ident="ID=gene",
    )
    assert [item.core.to_gff() for item in regions] == [
        (1, 10),
        (21, 30),
        (41, 50),
        (61, 100),
    ]


def test_gene_ident_type_lists_are_exact_and_prokka_alias_is_stable() -> None:
    cds = feature(11, 20, "CDS")
    rrna = feature(31, 40, "rRNA")
    repeat = feature(51, 60, "repeat_region")
    assert feature_matches_gene_ident(cds, "CDS,rRNA")
    assert feature_matches_gene_ident(rrna, "CDS,rRNA")
    assert not feature_matches_gene_ident(repeat, "CDS,rRNA")
    assert feature_matches_gene_ident(rrna, "Prokka")
    assert not feature_matches_gene_ident(repeat, "Prokka")
    assert resolve_gene_ident("Bakta") == resolve_gene_ident("Prokka")


def test_gene_ident_rejects_mixed_attribute_and_type_selectors() -> None:
    with pytest.raises(ValueError, match="cannot be mixed"):
        resolve_gene_ident("ID=gene,CDS")


def test_out_of_bounds_annotation_fails() -> None:
    with pytest.raises(ValueError, match="exceeds sequence length"):
        extract_unannotated_regions({"ctg": "A" * 100}, [feature(1, 101)])


def test_lowercase_and_uppercase_candidates_are_identical() -> None:
    sequence = storf_sequence()
    upper, _ = find_storf_candidates(region(sequence), min_orf=99)
    lower, _ = find_storf_candidates(region(sequence.lower()), min_orf=99)
    assert [(item.interval, item.strand, item.sequence) for item in upper] == [
        (item.interval, item.strand, item.sequence) for item in lower
    ]
    assert any(item.strand == "+" for item in upper)


def test_reverse_strand_coordinates_and_sequence() -> None:
    oriented = storf_sequence()
    genomic = "CCC" + reverse_complement(oriented) + "CCC"
    candidates, _ = find_storf_candidates(region(genomic), min_orf=99)
    negative = [item for item in candidates if item.strand == "-"]
    assert len(negative) == 1
    assert negative[0].interval == Interval(3, 3 + len(oriented))
    assert negative[0].sequence == oriented


def test_ambiguous_threshold_is_enforced() -> None:
    sequence = "TAA" + "ATG" + ("NNN" * 30) + "TAG"
    candidates, rejected = find_storf_candidates(region(sequence), min_orf=99, ambiguous_threshold=0.0)
    assert candidates == []
    assert rejected == 1


def test_short_storf_has_both_bounds() -> None:
    exactly_120 = "TAA" + "AAA" * 38 + "TAG"
    candidates, _ = find_storf_candidates(region(exactly_120), min_orf=300, max_orf=600, include_short=True)
    assert [item.length for item in candidates if item.strand == "+"] == [120]
    too_long = "TAA" + "AAA" * 39 + "TAG"
    candidates, _ = find_storf_candidates(region(too_long), min_orf=300, max_orf=600, include_short=True)
    assert not [item for item in candidates if item.strand == "+"]


def test_code_four_does_not_treat_tga_as_stop() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 10 + "TGA" + "AAA" * 20 + "TAG"
    code_11, _ = find_storf_candidates(region(sequence), genetic_code=11, min_orf=6)
    code_4, _ = find_storf_candidates(region(sequence), genetic_code=4, min_orf=6)
    plus_11 = [item.length for item in code_11 if item.strand == "+"]
    plus_4 = [item.length for item in code_4 if item.strand == "+"]
    assert len(plus_11) == 2
    assert plus_4 == [len(sequence)]
    assert translate("TGA", 4) == "W"


def test_alternative_initiation_codon_translates_as_methionine_when_requested() -> None:
    sequence = "GTGAAATAG"
    assert translate(sequence, 11, keep_terminal_stop=False) == "VK"
    assert translate(sequence, 11, keep_terminal_stop=False, start_as_methionine=True) == "MK"


def test_consecutive_storf_is_explicitly_labelled() -> None:
    sequence = "TAA" + "AAA" * 10 + "TGA" + "AAA" * 10 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6, include_consecutive=True)
    assert sum(item.candidate_type == "Con-StORF" for item in candidates) == 1


def test_consecutive_storf_can_span_multiple_middle_stops() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 9 + "TGA" + "AAA" * 10 + "TAG" + "AAA" * 10 + "TAA"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6, include_consecutive=True)
    consecutive = [item for item in candidates if item.candidate_type == "Con-StORF"]
    assert len(consecutive) == 1
    assert consecutive[0].middle_stop_codons == ("TGA", "TAG")
    assert consecutive[0].start_before_first_middle_stop


def test_consecutive_storf_does_not_select_a_start_after_first_middle_stop() -> None:
    sequence = "TAA" + "AAA" * 10 + "TGA" + "ATG" + "AAA" * 9 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6, include_consecutive=True)
    consecutive = next(item for item in candidates if item.candidate_type == "Con-StORF")
    feature = candidate_feature(consecutive, AnalysisConfig(min_orf=6))
    assert consecutive.evidence.start_codon is None
    assert feature.attributes["CandidateValidation"] == "failed"
    assert "before the first middle stop" in feature.attributes["CandidateValidationIssues"]


def test_start_selection_ranks_rbs_and_codon_context_instead_of_taking_first_start() -> None:
    sequence = "TAA" + "TTG" + "AAA" + "AGGAGG" + "AAA" + "ATG" + "AAA" * 30 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6)
    candidate = next(item for item in candidates if item.strand == "+")
    assert candidate.evidence.start_candidate_count == 2
    assert candidate.evidence.start_codon == "ATG"
    assert candidate.evidence.start_offset == 18
    assert candidate.evidence.rbs_motif == "AGGAGG"
    promoted = candidate_feature(
        candidate,
        AnalysisConfig(min_orf=6, promote_tier="exploratory"),
    )
    assert promoted.type == "CDS"
    assert promoted.interval == Interval(candidate.interval.start + 3, candidate.interval.end)
    assert promoted.interval != candidate.coding_interval
    assert promoted.attributes["CodingCoordinates"] == (
        f"{candidate.coding_interval.start + 1}-{candidate.coding_interval.end}"
    )


def test_rbs_search_uses_context_before_stop_to_stop_candidate() -> None:
    sequence = "AGGAGG" + "TAA" + "ATG" + "AAA" * 30 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6)
    candidate = next(item for item in candidates if item.strand == "+")
    assert candidate.interval.start == 6
    assert candidate.evidence.rbs_motif == "AGGAGG"
    assert candidate.evidence.rbs_distance == 3


def test_start_selection_rejects_too_short_rbs_supported_internal_start() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 27 + "AGGAGG" + "AAA" + "ATG" + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=99)
    candidate = next(item for item in candidates if item.strand == "+")
    assert candidate.evidence.start_offset == 3
    assert candidate.evidence.start_codon == "ATG"
    assert candidate.coding_length == 99


def test_evidence_is_transparent_and_cds_promotion_requires_start() -> None:
    with_start, _ = find_storf_candidates(region(storf_sequence()), min_orf=99)
    candidate = next(item for item in with_start if item.strand == "+")
    assert candidate.evidence.start_codon == "ATG"
    assert "length" in candidate.evidence.components
    promoted = candidate_feature(candidate, AnalysisConfig(feature_type="ORF", promote_tier="exploratory"))
    assert promoted.type == "CDS"
    assert promoted.phase == "0"
    assert promoted.interval.start == candidate.interval.start + 3

    no_start_sequence = "TAA" + "AAA" * 31 + "TAG"
    no_start, _ = find_storf_candidates(region(no_start_sequence), min_orf=99)
    candidate = next(item for item in no_start if item.strand == "+")
    unpromoted = candidate_feature(candidate, AnalysisConfig(feature_type="ORF", promote_tier="exploratory"))
    assert unpromoted.type == "ORF"
    assert unpromoted.phase == "."


def test_cds_promotion_rejects_consecutive_storf_with_internal_stop() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 10 + "TGA" + "AAA" * 10 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6, include_consecutive=True)
    consecutive = next(item for item in candidates if item.candidate_type == "Con-StORF")
    feature = candidate_feature(
        consecutive,
        AnalysisConfig(min_orf=6, feature_type="ORF", promote_tier="exploratory"),
    )
    assert feature.type == "ORF"
    assert feature.attributes["CandidateValidation"] == "passed"
    assert feature.attributes["CandidateValidationPolicy"] == "con-storf-readthrough-v1"
    assert feature.attributes["CDSValidation"] == "failed"
    assert "internal in-frame stop" in feature.attributes["CDSValidationIssues"]
    assert feature.attributes["MiddleStopCount"] == 1
    assert feature.attributes["StartBeforeFirstMiddleStop"] == "true"


def test_gff_coordinate_mode_and_explicit_length_attributes() -> None:
    candidates, _ = find_storf_candidates(region(storf_sequence()), min_orf=99)
    candidate = next(item for item in candidates if item.strand == "+")
    stop_to_stop = candidate_feature(candidate, AnalysisConfig(min_orf=99))
    stop_inclusive = candidate_feature(
        candidate,
        AnalysisConfig(min_orf=99, include_first_stop=True),
    )
    start_to_stop = candidate_feature(
        candidate,
        AnalysisConfig(min_orf=99, gff_coordinate_mode="start-to-stop"),
    )
    promoted = candidate_feature(
        candidate,
        AnalysisConfig(
            min_orf=99,
            gff_coordinate_mode="stop-to-stop",
            promote_tier="exploratory",
        ),
    )

    assert stop_to_stop.interval == Interval(candidate.interval.start + 3, candidate.interval.end)
    assert stop_to_stop.attributes["GFFCoordinateMode"] == "stop-to-stop-first-stop-excluded"
    assert stop_inclusive.interval == candidate.interval
    assert stop_inclusive.attributes["GFFCoordinateMode"] == "stop-to-stop-inclusive"
    assert start_to_stop.interval == candidate.coding_interval
    assert start_to_stop.attributes["GFFCoordinateMode"] == "start-to-stop"
    assert promoted.interval == Interval(candidate.interval.start + 3, candidate.interval.end)
    assert promoted.attributes["GFFCoordinateMode"] == "stop-to-stop-first-stop-excluded"
    assert "Length" not in stop_to_stop.attributes
    assert stop_to_stop.attributes["CandidateLength"] == candidate.length
    assert stop_to_stop.attributes["CodingLength"] == candidate.coding_length


def test_consecutive_storf_can_replace_overlapping_storfs_but_never_become_cds() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 30 + "TGA" + "AAA" * 31 + "TAG"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=99, include_consecutive=True)
    filtered, _ = apply_overlap_filter(candidates, maximum_overlap=0, mode="all")
    assert len(filtered) == 1
    assert filtered[0].candidate_type == "Con-StORF"
    feature = candidate_feature(filtered[0], AnalysisConfig(min_orf=99, promote_tier="exploratory"))
    assert feature.attributes["CandidateValidation"] == "passed"
    assert feature.type == "ORF"
    assert feature.attributes["PromotedToCDS"] == "false"


def test_overlap_filter_prefers_evidence_then_is_deterministic() -> None:
    sequence = "TAA" + "ATG" + "AAA" * 10 + "TAG" + "AAA" * 20 + "TGA"
    candidates, _ = find_storf_candidates(region(sequence), min_orf=6, include_consecutive=True)
    filtered_a, removed_a = apply_overlap_filter(candidates, maximum_overlap=0, mode="all")
    filtered_b, removed_b = apply_overlap_filter(list(reversed(candidates)), maximum_overlap=0, mode="all")
    assert [item.id for item in filtered_a] == [item.id for item in filtered_b]
    assert removed_a == removed_b


def test_overlap_filter_can_measure_coding_instead_of_stop_to_stop_intervals() -> None:
    candidates, _ = find_storf_candidates(region(storf_sequence()), min_orf=99)
    first = next(item for item in candidates if item.strand == "+")
    second = replace(
        first,
        id=f"{first.id}_shifted",
        interval=Interval(first.interval.end - 9, first.interval.end - 9 + first.length),
    )
    coding, _ = apply_overlap_filter(
        [first, second], maximum_overlap=6, mode="all", coordinate_basis="coding"
    )
    candidate, _ = apply_overlap_filter(
        [first, second], maximum_overlap=6, mode="all", coordinate_basis="candidate"
    )
    assert len(coding) == 2
    assert len(candidate) == 1
