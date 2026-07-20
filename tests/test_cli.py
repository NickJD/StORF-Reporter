from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path

import pytest

from StORF_Reporter.cli import build_parser
from StORF_Reporter.io import read_fasta, read_gff

ROOT = Path(__file__).parents[1]
ENV = {**os.environ, "PYTHONPATH": str(ROOT / "src")}
DIRECT_ENV = os.environ.copy()
DIRECT_ENV.pop("PYTHONPATH", None)


def run_module(module: str, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-m", module, *map(str, args)],
        cwd=cwd or ROOT,
        env=ENV,
        text=True,
        capture_output=True,
        check=False,
    )


def run_script(script: str, *args: str) -> subprocess.CompletedProcess[str]:
    package_dir = ROOT / "src" / "StORF_Reporter"
    return subprocess.run(
        [sys.executable, script, *map(str, args)],
        cwd=package_dir,
        env=DIRECT_ENV,
        text=True,
        capture_output=True,
        check=False,
    )


def read_sequence_records(path: Path) -> dict[str, str]:
    """Read nucleotide or protein FASTA without imposing an alphabet."""

    records: dict[str, str] = {}
    identifier: str | None = None
    for raw in path.read_text().splitlines():
        if raw.startswith(">"):
            identifier = raw[1:].split()[0]
            records[identifier] = ""
        elif raw and identifier is not None:
            records[identifier] += raw
    return records


@pytest.mark.parametrize(
    "module",
    [
        "StORF_Reporter.cli",
        "StORF_Reporter.StORF_Reporter",
        "StORF_Reporter.StORF_Finder",
        "StORF_Reporter.UR_Extractor",
        "StORF_Reporter.StORF_Extractor",
        "StORF_Reporter.StORF_Remover",
    ],
)
def test_versions_exit_zero_without_required_inputs(module: str) -> None:
    flag = "--version"
    result = run_module(module, flag)
    assert result.returncode == 0, result.stderr
    assert "2.0.0" in result.stdout


def test_cds_promotion_uses_default_or_user_selected_tier() -> None:
    parser = build_parser()
    default_args = parser.parse_args(["augment", "--input", "sample.gff"])
    promoted_args = parser.parse_args(["augment", "--input", "sample.gff", "--promote-cds", "high"])
    disabled_args = parser.parse_args(["augment", "--input", "sample.gff", "--no-promote-cds"])
    assert default_args.promote_cds == "high"
    assert promoted_args.promote_cds == "high"
    assert disabled_args.promote_cds is None
    assert default_args.gff_coordinate_mode == "stop-to-stop"
    assert default_args.include_first_stop is False
    assert default_args.protein_coordinate_mode == "candidate"
    assert default_args.remove_inframe_stops is False

    start_coordinates = parser.parse_args(
        [
            "augment",
            "--input",
            "sample.gff",
            "--gff-coordinate-mode",
            "start-to-stop",
        ]
    )
    assert start_coordinates.gff_coordinate_mode == "start-to-stop"

    inclusive = parser.parse_args(
        [
            "augment",
            "--input",
            "sample.gff",
            "--include-first-stop",
            "--protein-coordinate-mode",
            "coding",
            "--remove-inframe-stops",
        ]
    )
    assert inclusive.include_first_stop is True
    assert inclusive.protein_coordinate_mode == "coding"
    assert inclusive.remove_inframe_stops is True


@pytest.mark.parametrize(
    ("script", "flag", "expected"),
    [
        ("cli.py", "--version", "2.0.0"),
        ("StORF_Reporter.py", "--version", "2.0.0"),
        ("StORF_Finder.py", "--version", "2.0.0"),
        ("UR_Extractor.py", "--version", "2.0.0"),
        ("StORF_Extractor.py", "--version", "2.0.0"),
        ("StORF_Remover.py", "--version", "2.0.0"),
        ("Tools/StORF_Adder.py", "--help", "usage:"),
        ("Tools/UR_Lenghts.py", "--help", "usage:"),
        ("Tools/Un_StORFed.py", "--help", "usage:"),
    ],
)
def test_source_scripts_work_without_an_installed_package(script: str, flag: str, expected: str) -> None:
    result = run_script(script, flag)
    assert result.returncode == 0, result.stderr
    assert expected in result.stdout


def fixture(tmp_path: Path) -> tuple[Path, Path]:
    candidate = "TAA" + "ATG" + "AAA" * 30 + "TAG"
    sequence = "A" * 100 + candidate + "A" * 30
    fasta = tmp_path / "sample.fa"
    fasta.write_text(f">ctg\n{sequence}\n")
    gff = tmp_path / "sample.gff"
    gff.write_text("##gff-version 3\nctg\ttest\tCDS\t1\t100\t.\t+\t0\tID=gene1\n")
    return gff, fasta


def gene_ident_fixture(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "selector.fa"
    fasta.write_text(">ctg\n" + "A" * 100 + "\n")
    gff = tmp_path / "selector.gff"
    gff.write_text(
        "##gff-version 3\n"
        "ctg\ttest\tncRNA_gene\t31\t40\t.\t+\t.\tID=gene:ncrna\n"
        "ctg\ttest\tgene\t61\t70\t.\t+\t.\tID=not_gene_prefix\n"
    )
    return gff, fasta


def test_reporter_default_output_directory_no_longer_crashes(tmp_path: Path) -> None:
    gff, _ = fixture(tmp_path)
    result = run_module(
        "StORF_Reporter.StORF_Reporter",
        "-anno",
        "Feature_Types",
        "Single_Genome",
        "-p",
        str(gff),
        "-minorf",
        "99",
        "-olap_filt",
        "none",
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "sample_StORF-Reporter_Extended.gff").exists()


def test_standalone_finder_supports_standard_fasta_gzip_and_aa_only(tmp_path: Path) -> None:
    _, fasta = fixture(tmp_path)
    result = run_module(
        "StORF_Reporter.StORF_Finder",
        "-f",
        str(fasta),
        "-ua",
        "False",
        "-gz",
        "True",
        "-aa_only",
        "True",
        "-minorf",
        "99",
        "-olap_filt",
        "none",
        "-odir",
        str(tmp_path),
        "-oname",
        "finder",
    )
    assert result.returncode == 0, result.stderr
    output = tmp_path / "finder_aa.fasta.gz"
    assert output.read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(output, "rt") as handle:
        assert ">StORF_" in handle.read()


def test_standalone_finder_honours_gff_false(tmp_path: Path) -> None:
    _, fasta = fixture(tmp_path)
    result = run_module(
        "StORF_Reporter.StORF_Finder",
        "-f",
        str(fasta),
        "-ua",
        "False",
        "-gff",
        "False",
        "-minorf",
        "99",
        "-olap_filt",
        "none",
        "-odir",
        str(tmp_path),
        "-oname",
        "nogff",
    )
    assert result.returncode == 0, result.stderr
    assert not (tmp_path / "nogff.gff").exists()
    assert (tmp_path / "nogff.fasta").exists()


def test_unified_augment_dry_run_and_execution(tmp_path: Path) -> None:
    gff, fasta = fixture(tmp_path)
    dry = run_module("StORF_Reporter.cli", "augment", "--input", str(gff), "--fasta", str(fasta), "--dry-run")
    assert dry.returncode == 0, dry.stderr
    assert str(gff.resolve()) in dry.stdout
    run = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(tmp_path / "results"),
        "--overlap-mode",
        "none",
        "--promote-cds",
        "exploratory",
    )
    assert run.returncode == 0, run.stderr
    sample_dir = tmp_path / "results" / "sample"
    assert (sample_dir / "sample_StORF-Reporter.gff").exists()
    assert (sample_dir / "sample_manifest.json").exists()
    assert (sample_dir / "sample_report.html").exists()
    assert (sample_dir / "sample_UR.gff").exists()
    assert (sample_dir / "sample_UR.fasta").exists()
    assert (sample_dir / "sample_validation.json").exists()
    manifest = __import__("json").loads((sample_dir / "sample_manifest.json").read_text())
    assert manifest["configuration"]["gene_ident"] == "CDS"
    document = read_gff(sample_dir / "sample_StORF-Reporter.gff")
    reported = [feature for feature in document.features if feature.source == "StORF-Reporter"]
    assert reported and reported[0].type == "CDS"
    assert reported[0].attributes["CDSValidation"] == "passed"
    unsafe = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--feature-type",
        "CDS",
    )
    assert unsafe.returncode == 2
    assert "invalid choice" in unsafe.stderr


def test_gene_ident_attribute_mode_is_shared_by_unified_and_legacy_ur_commands(
    tmp_path: Path,
) -> None:
    gff, fasta = gene_ident_fixture(tmp_path)
    unified_dir = tmp_path / "unified"
    unified = run_module(
        "StORF_Reporter.cli",
        "extract-urs",
        "--gff",
        str(gff),
        "--fasta",
        str(fasta),
        "--gene-ident",
        "ID=gene",
        "--min-ur",
        "1",
        "--extension",
        "0",
        "--outdir",
        str(unified_dir),
    )
    assert unified.returncode == 0, unified.stderr

    profile_dir = tmp_path / "profile"
    profile = run_module(
        "StORF_Reporter.cli",
        "extract-urs",
        "--gff",
        str(gff),
        "--fasta",
        str(fasta),
        "--profile",
        "ensembl",
        "--min-ur",
        "1",
        "--extension",
        "0",
        "--outdir",
        str(profile_dir),
    )
    assert profile.returncode == 0, profile.stderr

    legacy_dir = tmp_path / "legacy"
    legacy = run_module(
        "StORF_Reporter.UR_Extractor",
        "-gff",
        str(gff),
        "-f",
        str(fasta),
        "-gene_ident",
        "ID=gene",
        "-min_len",
        "1",
        "-ex_len",
        "0",
        "-odir",
        str(legacy_dir),
    )
    assert legacy.returncode == 0, legacy.stderr

    reporter_dir = tmp_path / "reporter"
    reporter = run_module(
        "StORF_Reporter.StORF_Reporter",
        "-anno",
        "Ensembl",
        "Single_Genome",
        "-p",
        str(gff),
        "-odir",
        str(reporter_dir),
        "--overwrite",
    )
    assert reporter.returncode == 0, reporter.stderr

    expected = [(1, 30), (41, 100)]
    unified_regions = read_gff(unified_dir / "selector_UR.gff").features
    profile_regions = read_gff(profile_dir / "selector_UR.gff").features
    legacy_regions = read_gff(legacy_dir / "selector_UR.gff").features
    assert [feature.interval.to_gff() for feature in unified_regions] == expected
    assert [feature.interval.to_gff() for feature in profile_regions] == expected
    assert [feature.interval.to_gff() for feature in legacy_regions] == expected
    reporter_manifest = __import__("json").loads(
        (reporter_dir / "selector_StORF-Reporter_Extended_manifest.json").read_text()
    )
    assert reporter_manifest["configuration"]["gene_ident"] == "ID=gene"


def test_augment_rejects_invalid_input_before_analysis(tmp_path: Path) -> None:
    gff, fasta = fixture(tmp_path)
    gff.write_text(gff.read_text().replace("\t+\t0\t", "\t+\t.\t"))
    result = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(tmp_path / "invalid-results"),
    )
    assert result.returncode == 2
    assert "input validation failed" in result.stderr
    assert not (tmp_path / "invalid-results").exists()


def test_augment_partitions_candidates_below_primary_tier(tmp_path: Path) -> None:
    gff, fasta = fixture(tmp_path)
    outdir = tmp_path / "tiered-results"
    result = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(outdir),
        "--overlap-mode",
        "none",
        "--primary-tier",
        "high",
    )
    assert result.returncode == 0, result.stderr
    sample_dir = outdir / "sample"
    primary = read_gff(sample_dir / "sample_StORF-Reporter.gff")
    assert not [feature for feature in primary.features if feature.source == "StORF-Reporter"]
    review = read_gff(sample_dir / "sample_StORFs_below-high.gff")
    assert [feature for feature in review.features if feature.source == "StORF-Reporter"]
    assert (sample_dir / "sample_StORFs.fasta").read_text() == ""
    assert ">StORF_" in (sample_dir / "sample_StORFs_below-high.fasta").read_text()
    validation = __import__("json").loads((sample_dir / "sample_validation.json").read_text())
    partition = validation["candidate_cds_validation"]
    assert partition["primary_output"] == 0
    assert partition["below_primary_tier_output"] >= 1


def test_filter_can_promote_homology_supported_validated_candidate(tmp_path: Path) -> None:
    storf = "TAA" + "ATG" + "AAA" * 30 + "TAG"
    sequence = "A" * 100 + storf + "A" * 30 + storf + "A" * 30
    fasta = tmp_path / "sample.fa"
    fasta.write_text(f">ctg\n{sequence}\n")
    gff = tmp_path / "sample.gff"
    gff.write_text("##gff-version 3\nctg\ttest\tCDS\t1\t100\t.\t+\t0\tID=gene1\n")
    augmented_dir = tmp_path / "augmented"
    augment = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(augmented_dir),
        "--overlap-mode",
        "none",
    )
    assert augment.returncode == 0, augment.stderr
    augmented_gff = augmented_dir / "sample" / "sample_StORF-Reporter.gff"
    candidates = [
        feature for feature in read_gff(augmented_gff).features if feature.source == "StORF-Reporter"
    ]
    assert len(candidates) >= 2
    candidate = candidates[0]
    alignments = tmp_path / "hits.tsv"
    alignments.write_text(
        f"{candidate.attributes['ID']};UR=historical_header\tknown_protein\t95.0\t30\t0\t0\t"
        "1\t30\t1\t30\t1e-30\t80.0\n"
    )
    filtered = tmp_path / "homology-supported.gff"
    result = run_module(
        "StORF_Reporter.cli",
        "filter",
        "--gff",
        str(augmented_gff),
        "--fasta",
        str(fasta),
        "--alignments",
        str(alignments),
        "--promote-cds",
        "--keep-unmatched",
        "--out",
        str(filtered),
    )
    assert result.returncode == 0, result.stderr
    reported = [feature for feature in read_gff(filtered).features if feature.source == "StORF-Reporter"]
    promoted = next(feature for feature in reported if feature.attributes["ID"] == candidate.attributes["ID"])
    assert promoted.type == "CDS"
    assert promoted.attributes["PromotionEvidence"] == "homology-format6"
    assert promoted.attributes["HomologySubject"] == "known_protein"
    unmatched = [feature for feature in reported if feature.attributes["ID"] != candidate.attributes["ID"]]
    assert unmatched and all(feature.type == "ORF" for feature in unmatched)


def test_filter_does_not_promote_con_storf_from_homology_alone(tmp_path: Path) -> None:
    con_storf = "TAA" + "ATG" + "AAA" * 30 + "TGA" + "AAA" * 31 + "TAG"
    sequence = "A" * 100 + con_storf + "A" * 30
    fasta = tmp_path / "con.fa"
    fasta.write_text(f">ctg\n{sequence}\n")
    gff = tmp_path / "con.gff"
    gff.write_text("##gff-version 3\nctg\ttest\tCDS\t1\t100\t.\t+\t0\tID=gene1\n")
    augmented_dir = tmp_path / "con-augmented"
    augment = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(augmented_dir),
        "--overlap-mode",
        "none",
        "--include-consecutive",
        "--no-promote-cds",
    )
    assert augment.returncode == 0, augment.stderr
    augmented_gff = augmented_dir / "con" / "con_StORF-Reporter.gff"
    consecutive = next(
        feature
        for feature in read_gff(augmented_gff).features
        if feature.attributes.get("StORF_Type") == "Con-StORF"
    )
    alignments = tmp_path / "con-hits.tsv"
    alignments.write_text(
        f"{consecutive.attributes['ID']}\tknown_protein\t95.0\t60\t0\t0\t1\t60\t1\t60\t1e-30\t80.0\n"
    )
    filtered = tmp_path / "con-filtered.gff"
    result = run_module(
        "StORF_Reporter.cli",
        "filter",
        "--gff",
        str(augmented_gff),
        "--fasta",
        str(fasta),
        "--alignments",
        str(alignments),
        "--promote-cds",
        "--out",
        str(filtered),
    )
    assert result.returncode == 0, result.stderr
    retained = next(
        feature
        for feature in read_gff(filtered).features
        if feature.attributes.get("ID") == consecutive.attributes["ID"]
    )
    assert retained.type == "ORF"
    assert retained.attributes["HomologyCDSPromotion"] == "failed"
    assert "readthrough" in retained.attributes["HomologyCDSPromotionIssues"]


def test_extract_auto_round_trips_candidate_and_coding_sequences(tmp_path: Path) -> None:
    gff, fasta = fixture(tmp_path)
    results = tmp_path / "results"
    augment = run_module(
        "StORF_Reporter.cli",
        "augment",
        "--input",
        str(gff),
        "--fasta",
        str(fasta),
        "--outdir",
        str(results),
        "--overlap-mode",
        "none",
        "--no-promote-cds",
        "--embed-fasta",
    )
    assert augment.returncode == 0, augment.stderr
    sample = results / "sample"
    augmented = sample / "sample_StORF-Reporter.gff"

    extracted_nt = tmp_path / "extracted.fasta"
    extracted_aa = tmp_path / "extracted_aa.fasta"
    for output, extra in ((extracted_nt, ()), (extracted_aa, ("--amino-acid",))):
        result = run_module(
            "StORF_Reporter.cli",
            "extract",
            "--gff",
            str(augmented),
            "--out",
            str(output),
            *extra,
        )
        assert result.returncode == 0, result.stderr

    assert read_fasta(extracted_nt) == read_fasta(sample / "sample_StORFs.fasta")
    assert read_sequence_records(extracted_aa) == read_sequence_records(sample / "sample_StORFs_aa.fasta")

    legacy_aa = run_module(
        "StORF_Reporter.StORF_Extractor",
        "-storf_input",
        "Combined",
        "-p",
        str(augmented),
        "-aa",
        "True",
        "-stop_ident",
        "False",
        "-odir",
        str(tmp_path),
        "-oname",
        "legacy",
    )
    assert legacy_aa.returncode == 0, legacy_aa.stderr
    assert read_sequence_records(tmp_path / "legacy_aa.fasta") == read_sequence_records(
        sample / "sample_StORFs_aa.fasta"
    )


def test_legacy_short_only_and_overlap_mode_are_not_silently_ignored(tmp_path: Path) -> None:
    sequence = "TAA" + "AAA" * 10 + "TAG" + "AAA" * 40 + "TGA"
    fasta = tmp_path / "short.fa"
    fasta.write_text(f">ctg\n{sequence}\n")
    result = run_module(
        "StORF_Reporter.StORF_Finder",
        "-f",
        str(fasta),
        "-ua",
        "False",
        "-short_storfs",
        "Nolap",
        "-short_storfs_only",
        "True",
        "-olap_filt",
        "none",
        "-odir",
        str(tmp_path),
        "-oname",
        "short-only",
    )
    assert result.returncode == 0, result.stderr
    reported = read_gff(tmp_path / "short-only.gff").features
    assert reported
    assert {feature.attributes["StORF_Type"] for feature in reported} == {"Short-StORF"}


def test_storf_adder_uses_total_overlap_and_coordinate_order(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">ctg\n" + "A" * 200 + "\n")
    reference = tmp_path / "reference.gff"
    reference.write_text(
        "##gff-version 3\n"
        "ctg\ttest\tCDS\t1\t30\t.\t+\t0\tID=left1\n"
        "ctg\ttest\tCDS\t41\t70\t.\t+\t0\tID=left2\n"
        "ctg\ttest\tCDS\t120\t130\t.\t+\t0\tID=right\n"
    )
    additional = tmp_path / "additional.gff"
    additional.write_text(
        "##gff-version 3\n"
        "ctg\tStORF-Reporter\tORF\t20\t51\t.\t+\t.\tID=StORF_rejected\n"
        "ctg\tStORF-Reporter\tORF\t80\t90\t.\t+\t.\tID=StORF_kept\n"
    )
    output = tmp_path / "merged.gff"
    result = run_module(
        "StORF_Reporter.Tools.StORF_Adder",
        "--genome_DNA",
        str(fasta),
        "--reference_annotation",
        str(reference),
        "--additional_annotation",
        str(additional),
        "--overlap",
        "15",
        "--output_file",
        str(output),
    )
    assert result.returncode == 0, result.stderr
    features = read_gff(output).features
    ids = [feature.attributes["ID"] for feature in features]
    assert "StORF_rejected" not in ids
    assert ids == ["left1", "left2", "StORF_kept", "right"]


def test_legacy_reporter_numbers_custom_names_for_multiple_inputs(tmp_path: Path) -> None:
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    for sample in ("one", "two"):
        sequence = "A" * 100 + "TAA" + "ATG" + "AAA" * 30 + "TAG"
        (inputs / f"{sample}.fa").write_text(f">ctg_{sample}\n{sequence}\n")
        (inputs / f"{sample}.gff").write_text(
            f"##gff-version 3\nctg_{sample}\ttest\tCDS\t1\t100\t.\t+\t0\tID={sample}_gene\n"
        )
    outputs = tmp_path / "outputs"
    result = run_module(
        "StORF_Reporter.StORF_Reporter",
        "-anno",
        "Feature_Types",
        "Multiple_Genomes",
        "-p",
        str(inputs),
        "-odir",
        str(outputs),
        "-oname",
        "batch",
        "-olap_filt",
        "none",
    )
    assert result.returncode == 0, result.stderr
    assert (outputs / "batch_1.gff").exists()
    assert (outputs / "batch_2.gff").exists()
