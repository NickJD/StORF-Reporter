from __future__ import annotations

import gzip
from collections import OrderedDict
from pathlib import Path

import pytest

from StORF_Reporter.genetic_code import translate
from StORF_Reporter.io import InputFormatError, read_fasta, read_gff
from StORF_Reporter.models import Interval
from StORF_Reporter.validation import validate_inputs
from StORF_Reporter.workflow import (
    AnalysisConfig,
    analyse,
    discover_pairs,
    extract_feature_sequence,
    load_annotated_genome,
    pyrodigal_features,
    write_augmented_gff,
    write_candidate_outputs,
    write_combined_cds_storf_fasta,
    write_manifest,
)


def write_fixture(tmp_path: Path, lowercase: bool = False) -> tuple[Path, Path]:
    candidate = "TAA" + "ATG" + "AAA" * 30 + "TAG"
    sequence = "A" * 100 + candidate + "A" * 50
    if lowercase:
        sequence = sequence.lower()
    fasta = tmp_path / "genome.fna"
    fasta.write_text(f">ctg\n{sequence}\n", encoding="utf-8")
    gff = tmp_path / "genome.gff"
    gff.write_text(
        "##gff-version 3\n"
        f"##sequence-region ctg 1 {len(sequence)}\n"
        "ctg\ttest\tCDS\t1\t100\t.\t+\t0\tID=gene1\n",
        encoding="utf-8",
    )
    return gff, fasta


def test_parser_accepts_feature_as_second_line(tmp_path: Path) -> None:
    gff = tmp_path / "second-line.gff"
    gff.write_text(
        "##gff-version 3\nctg\ttest\tCDS\t1\t3\t.\t+\t0\tID=x\n",
        encoding="utf-8",
    )
    assert len(read_gff(gff).features) == 1


def test_duplicate_fasta_ids_fail(tmp_path: Path) -> None:
    fasta = tmp_path / "duplicate.fa"
    fasta.write_text(">x\nAAA\n>x\nCCC\n", encoding="utf-8")
    with pytest.raises(InputFormatError, match="duplicate"):
        read_fasta(fasta)


def test_non_iupac_fasta_fails_instead_of_producing_empty_output(tmp_path: Path) -> None:
    fasta = tmp_path / "invalid.fa"
    fasta.write_text(">x\nACGT7\n", encoding="utf-8")
    with pytest.raises(InputFormatError, match="non-IUPAC"):
        read_fasta(fasta)


def test_gzip_is_detected_by_magic_not_name(tmp_path: Path) -> None:
    fasta = tmp_path / "misnamed.fa"
    with gzip.open(fasta, "wt") as handle:
        handle.write(">x\naaa\n")
    assert read_fasta(fasta) == OrderedDict([("x", "AAA")])


def test_final_unannotated_region_is_analysed_and_written(tmp_path: Path) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, overlap_mode="none")
    result = analyse(
        sequences,
        document.features,
        config,
        document=document,
        input_gff=gff,
        input_fasta=fasta,
        sample="fixture",
    )
    assert any(candidate.interval.start >= 100 for candidate in result.candidates)
    output = tmp_path / "augmented.gff"
    write_augmented_gff(result, config, output)
    parsed = read_gff(output)
    assert any(feature.source == "StORF-Reporter" for feature in parsed.features)
    assert not [issue for issue in validate_inputs(output, fasta) if issue.severity == "error"]


def test_augmented_gff_interleaves_candidates_by_contig_coordinate(tmp_path: Path) -> None:
    candidate = "TAA" + "ATG" + "AAA" * 30 + "TAG"
    sequence = "A" * 100 + candidate + "A" * 30
    fasta = tmp_path / "ordered.fna"
    fasta.write_text(f">ctg\n{sequence}\n")
    gff = tmp_path / "ordered.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "ctg\ttest\tCDS\t1\t100\t.\t+\t0\tID=left\n"
        f"ctg\ttest\tCDS\t220\t{len(sequence)}\t.\t+\t0\tID=right\n"
    )
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, extension=0, overlap_mode="none")
    result = analyse(sequences, document.features, config, document=document, sample="ordered")
    output = tmp_path / "ordered-augmented.gff3"
    write_augmented_gff(result, config, output)
    features = read_gff(output).features
    starts = [feature.interval.start for feature in features]
    assert starts == sorted(starts)
    storf_index = next(index for index, feature in enumerate(features) if feature.source == "StORF-Reporter")
    assert features[storf_index - 1].attributes.get("ID") == "left"
    assert features[storf_index + 1].attributes.get("ID") == "right"


def test_all_outputs_are_real_gzip(tmp_path: Path) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, overlap_mode="none")
    result = analyse(sequences, document.features, config, document=document, sample="fixture")
    augmented = tmp_path / "augmented.gff.gz"
    write_augmented_gff(result, config, augmented)
    outputs = write_candidate_outputs(result, tmp_path / "fixture", gzip_output=True)
    for path in [augmented, *outputs.values()]:
        assert path.read_bytes()[:2] == b"\x1f\x8b"
        with gzip.open(path, "rt") as handle:
            assert handle.read(1)


def test_default_candidate_outputs_keep_stop_derived_locus_but_omit_first_stop(
    tmp_path: Path,
) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, overlap_mode="none")
    result = analyse(sequences, document.features, config, document=document, sample="fixture")
    candidate = next(item for item in result.candidates if item.strand == "+")

    augmented = tmp_path / "default.gff"
    write_augmented_gff(result, config, augmented)
    feature = next(item for item in read_gff(augmented).features if item.attributes.get("ID") == candidate.id)
    assert feature.type in {"ORF", "CDS"}
    assert feature.interval == Interval(candidate.interval.start + 3, candidate.interval.end)
    assert feature.attributes["CandidateCoordinates"] == (
        f"{candidate.interval.start + 1}-{candidate.interval.end}"
    )

    outputs = write_candidate_outputs(result, tmp_path / "default")
    assert read_fasta(outputs["nucleotide_fasta"])[candidate.id] == candidate.sequence[3:]

    inclusive_config = AnalysisConfig(min_orf=99, overlap_mode="none", include_first_stop=True)
    inclusive = tmp_path / "inclusive.gff"
    write_augmented_gff(result, inclusive_config, inclusive)
    inclusive_feature = next(
        item for item in read_gff(inclusive).features if item.attributes.get("ID") == candidate.id
    )
    assert inclusive_feature.interval == candidate.interval
    inclusive_outputs = write_candidate_outputs(
        result,
        tmp_path / "inclusive",
        include_leading_stop=True,
    )
    assert read_fasta(inclusive_outputs["nucleotide_fasta"])[candidate.id] == candidate.sequence


def test_promoted_cds_type_keeps_stop_derived_gff_locus_and_validates(tmp_path: Path) -> None:
    sequence = "TAA" + "TTG" + "AAA" + "AGGAGG" + "AAA" + "ATG" + "AAA" * 30 + "TAG"
    fasta = tmp_path / "ranked-start.fna"
    fasta.write_text(f">ctg\n{sequence}\n")
    config = AnalysisConfig(
        min_ur=1,
        max_ur=len(sequence),
        extension=0,
        min_orf=6,
        overlap_mode="none",
        promote_tier="exploratory",
    )
    result = analyse(OrderedDict([("ctg", sequence)]), [], config, sample="ranked-start")
    candidate = next(item for item in result.candidates if item.strand == "+")
    assert candidate.coding_interval is not None
    assert candidate.coding_interval.start > candidate.interval.start + 3

    output = tmp_path / "ranked-start.gff3"
    write_augmented_gff(result, config, output)
    promoted = next(item for item in read_gff(output).features if item.attributes.get("ID") == candidate.id)
    assert promoted.type == "CDS"
    assert promoted.interval == Interval(candidate.interval.start + 3, candidate.interval.end)
    assert promoted.interval != candidate.coding_interval
    assert extract_feature_sequence(promoted, sequence, "coding") == candidate.coding_sequence
    assert not [issue for issue in validate_inputs(output, fasta) if issue.severity == "error"]

    default_outputs = write_candidate_outputs(result, tmp_path / "candidate-protein")
    coding_outputs = write_candidate_outputs(
        result,
        tmp_path / "coding-protein",
        protein_coordinate_mode="coding",
    )

    def proteins(path: Path) -> dict[str, str]:
        records: dict[str, str] = {}
        identifier: str | None = None
        for line in path.read_text().splitlines():
            if line.startswith(">"):
                identifier = line[1:].split()[0]
                records[identifier] = ""
            elif identifier is not None:
                records[identifier] += line
        return records

    candidate_protein = proteins(default_outputs["protein_fasta"])[candidate.id]
    coding_protein = proteins(coding_outputs["protein_fasta"])[candidate.id]
    assert candidate_protein == translate(
        candidate.sequence[3:],
        candidate.genetic_code,
        keep_terminal_stop=False,
        start_as_methionine=False,
    )
    assert coding_protein == translate(
        candidate.coding_sequence,
        candidate.genetic_code,
        keep_terminal_stop=False,
        start_as_methionine=True,
    )
    assert candidate_protein != coding_protein


def test_combined_nucleotide_fasta_never_substitutes_ranked_coding_sequence(
    tmp_path: Path,
) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    result = analyse(
        sequences,
        document.features,
        AnalysisConfig(min_orf=99, overlap_mode="none"),
        document=document,
        sample="fixture",
    )
    candidate = next(item for item in result.candidates if item.coding_sequence is not None)
    path = write_combined_cds_storf_fasta(result, tmp_path / "combined.fasta")
    assert read_fasta(path)[candidate.id] == candidate.sequence[3:]


def test_internal_stop_removal_is_opt_in_and_protein_only(tmp_path: Path) -> None:
    sequence = "TAA" + "ATG" + "AAA" * 30 + "TGA" + "AAA" * 31 + "TAG"
    result = analyse(
        OrderedDict([("ctg", sequence)]),
        [],
        AnalysisConfig(
            min_ur=1,
            max_ur=len(sequence),
            extension=0,
            min_orf=99,
            overlap_mode="none",
            include_consecutive=True,
            consecutive_only=True,
        ),
        sample="con",
    )
    candidate = result.candidates[0]
    preserved = write_candidate_outputs(result, tmp_path / "preserved")
    removed = write_candidate_outputs(
        result,
        tmp_path / "removed",
        remove_inframe_stops=True,
    )

    def protein(path: Path) -> str:
        return "".join(line for line in path.read_text().splitlines() if not line.startswith(">"))

    assert "*" in protein(preserved["protein_fasta"])
    assert "*" not in protein(removed["protein_fasta"])
    assert read_fasta(preserved["nucleotide_fasta"])[candidate.id] == candidate.sequence[3:]
    assert read_fasta(removed["nucleotide_fasta"])[candidate.id] == candidate.sequence[3:]


def test_manifest_has_hashes_configuration_and_scientific_notice(tmp_path: Path) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, overlap_mode="none")
    result = analyse(
        sequences,
        document.features,
        config,
        document=document,
        input_gff=gff,
        input_fasta=fasta,
        sample="fixture",
    )
    path = tmp_path / "manifest.json"
    write_manifest(result, config, path, {"gff": tmp_path / "out.gff"})
    payload = __import__("json").loads(path.read_text())
    assert len(payload["inputs"]["gff_sha256"]) == 64
    assert payload["configuration"]["genetic_code"] == 11
    assert "not" in payload["scientific_notice"]


def test_pyrodigal_individual_mode_works_on_short_contig() -> None:
    pytest.importorskip("pyrodigal")
    sequences = OrderedDict([("ctg", "ATG" + "AAA" * 100 + "TAA")])
    features = pyrodigal_features(sequences, training="individual")
    assert isinstance(features, list)


def test_directory_discovery_ignores_previous_tool_outputs(tmp_path: Path) -> None:
    gff, _ = write_fixture(tmp_path)
    (tmp_path / "genome_UR.gff").write_text(gff.read_text())
    (tmp_path / "genome_StORFs_below-high.gff").write_text(gff.read_text())
    (tmp_path / "genome_UR.fasta").write_text(">ctg_UR_1_10\nAAAAAAAAAA\n")
    pairs = discover_pairs(tmp_path)
    assert [(sample, path.name) for sample, path, _ in pairs] == [("genome", "genome.gff")]


def test_directory_discovery_pairs_compressed_gff_and_fasta(tmp_path: Path) -> None:
    gff, fasta = write_fixture(tmp_path)
    compressed_gff = tmp_path / "compressed.gff3.gz"
    compressed_fasta = tmp_path / "compressed.fna.gz"
    with gzip.open(compressed_gff, "wt") as handle:
        handle.write(gff.read_text())
    with gzip.open(compressed_fasta, "wt") as handle:
        handle.write(fasta.read_text())
    gff.unlink()
    fasta.unlink()
    pairs = discover_pairs(tmp_path)
    assert pairs == [("compressed", compressed_gff.resolve(), compressed_fasta.resolve())]


def test_validator_rejects_structurally_invalid_storf_cds(tmp_path: Path) -> None:
    fasta = tmp_path / "invalid-storf.fna"
    fasta.write_text(">ctg\nAAA" + "AAA" * 30 + "TAG\n")
    gff = tmp_path / "invalid-storf.gff3"
    gff.write_text(
        "##gff-version 3\n"
        f"ctg\tStORF-Reporter\tCDS\t1\t{len('AAA' + 'AAA' * 30 + 'TAG')}\t.\t+\t0\t"
        "ID=StORF_invalid;GeneticCode=11\n"
    )
    issues = validate_inputs(gff, fasta)
    assert any(issue.code == "storf_cds_start" for issue in issues)


def test_validator_checks_v2_candidate_coordinate_and_length_attributes(tmp_path: Path) -> None:
    gff, fasta = write_fixture(tmp_path)
    document, sequences = load_annotated_genome(gff, fasta)
    config = AnalysisConfig(min_orf=99, overlap_mode="none", promote_tier=None)
    result = analyse(sequences, document.features, config, document=document, sample="fixture")
    output = tmp_path / "augmented.gff"
    write_augmented_gff(result, config, output)
    reported = next(feature for feature in read_gff(output).features if feature.source == "StORF-Reporter")
    length = reported.attributes["CandidateLength"]
    output.write_text(output.read_text().replace(f"CandidateLength={length}", "CandidateLength=1", 1))
    issues = validate_inputs(output, fasta)
    assert any(issue.code == "storf_candidate_length" for issue in issues)
