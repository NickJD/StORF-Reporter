"""Unified StORF-Reporter v2 command-line interface."""

from __future__ import annotations

import argparse
import json
import math
import platform
import shutil
import sys
from pathlib import Path

if __package__ in {None, ""}:
    # Support ``python cli.py`` from a source checkout while retaining normal
    # package-relative imports for wheels, conda packages, and ``python -m``.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .core import extract_unannotated_regions, resolve_gene_ident
from .genetic_code import reverse_complement, supported_codes, translate
from .io import InputFormatError, feature_line, open_text, read_fasta, read_gff, write_fasta_records
from .models import Feature, Interval, RunSummary
from .validation import validate_coding_sequence, validate_inputs
from .workflow import (
    AnalysisConfig,
    AnalysisResult,
    analyse,
    atomic_output,
    candidate_feature,
    discover_pairs,
    ensure_outputs_available,
    extract_feature_sequence,
    load_annotated_genome,
    logical_stem,
    partition_candidates,
    read_external_evidence,
    read_sample_manifest,
    storf_query_id,
    write_augmented_gff,
    write_candidate_gff,
    write_candidate_outputs,
    write_html_report,
    write_manifest,
    write_region_outputs,
)

CITATION = (
    "Dimonaco NJ, Aubrey W, Kenobi K, Clare A. StORF-Reporter: finding genes "
    "missed by gene finders. Nucleic Acids Research. 2023;51(21):11504-11517. "
    "doi:10.1093/nar/gkad814"
)

GENE_IDENT_PROFILES = {
    "default": "CDS",
    "prokka": "Prokka",
    "bakta": "Bakta",
    "ensembl": "ID=gene",
    "cds-only": "CDS",
}


def _csv_types(value: str) -> tuple[str, ...]:
    types = tuple(item.strip() for item in value.split(",") if item.strip())
    if not types:
        raise argparse.ArgumentTypeError("provide at least one feature type")
    if any("=" in feature_type for feature_type in types):
        raise argparse.ArgumentTypeError("attribute queries must use --gene-ident, for example ID=gene")
    return types


def _gene_ident_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--profile",
        choices=("default", "prokka", "bakta", "ensembl", "cds-only"),
        default="default",
        help="annotation selector shorthand; --gene-ident or --exclude-types overrides it",
    )
    selector = parser.add_mutually_exclusive_group()
    selector.add_argument(
        "--gene-ident",
        "--gene_ident",
        help=(
            "annotated-region selector: an exact comma-separated feature-type list, Prokka/Bakta, "
            "or one attribute prefix query such as ID=gene (default profile: CDS)"
        ),
    )
    selector.add_argument(
        "--exclude-types",
        type=_csv_types,
        default=None,
        help="compatibility alias for a comma-separated exact feature-type --gene-ident value",
    )


def _analysis_arguments(parser: argparse.ArgumentParser) -> None:
    _gene_ident_arguments(parser)
    parser.add_argument("--min-ur", type=int, default=30)
    parser.add_argument("--max-ur", type=int, default=100_000)
    parser.add_argument("--extension", type=int, default=50)
    parser.add_argument("--min-orf", type=int, default=99)
    parser.add_argument("--max-orf", type=int, default=60_000)
    parser.add_argument(
        "--genetic-code",
        type=int,
        choices=supported_codes(),
        default=11,
        help="NCBI translation table; supported: 1, 4, 11",
    )
    parser.add_argument("--ambiguous-threshold", type=float, default=0.20)
    parser.add_argument("--include-short", action="store_true")
    parser.add_argument(
        "--short-only",
        action="store_true",
        help="report only 30–120 nt Short-StORFs; implies --include-short",
    )
    parser.add_argument(
        "--short-overlap",
        choices=("exclude", "allow"),
        default="exclude",
        help=(
            "exclude conflicts between short and ordinary StORFs (historical Nolap), or allow "
            "cross-class overlap while filtering each class independently (historical Olap)"
        ),
    )
    parser.add_argument("--short-min", type=int, default=30)
    parser.add_argument("--short-max", type=int, default=120)
    parser.add_argument("--require-start", action="store_true")
    parser.add_argument(
        "--include-consecutive",
        action="store_true",
        help=(
            "enable maximal Con-StORFs spanning one or more middle in-frame stops; Con-StORFs "
            "enter the same overlap ranking and can replace overlapping ordinary StORFs"
        ),
    )
    parser.add_argument(
        "--consecutive-only",
        action="store_true",
        help="only report maximal Con-StORFs spanning one or more middle in-frame stops",
    )
    parser.add_argument("--overlap-mode", choices=("none", "same-strand", "all"), default="all")
    parser.add_argument(
        "--overlap-coordinate-basis",
        choices=("coding", "candidate"),
        default="coding",
        help=(
            "measure candidate-candidate conflicts using ranked coding intervals by default; "
            "use candidate for historical full stop-to-stop overlap filtering"
        ),
    )
    parser.add_argument("--maximum-overlap", type=int, default=50)
    parser.add_argument(
        "--maximum-annotation-overlap",
        type=int,
        default=50,
        help="maximum candidate nucleotides allowed inside annotated extension flanks",
    )
    parser.add_argument(
        "--feature-type",
        choices=("ORF", "StORF"),
        default="ORF",
        help="GFF3 type for unpromoted candidates (default: ORF)",
    )
    parser.add_argument(
        "--gff-coordinate-mode",
        choices=("stop-to-stop", "start-to-stop"),
        default="stop-to-stop",
        help=(
            "coordinates for reported candidates: stop-to-stop-derived loci "
            "with the upstream delimiting stop omitted by default, or the ranked start-to-terminal-"
            "stop interval when a start exists; promotion changes feature type, not these coordinates"
        ),
    )
    parser.add_argument(
        "--include-first-stop",
        "--stop-inclusive",
        action="store_true",
        help=(
            "include the upstream delimiting stop in GFF3 and candidate-locus FASTA output; "
            "by default it remains in CandidateCoordinates but is omitted from displayed/output sequence"
        ),
    )
    parser.add_argument(
        "--remove-inframe-stops",
        action="store_true",
        help=(
            "remove '*' characters produced by internal in-frame stops from candidate protein FASTA; "
            "this does not alter genomic coordinates or nucleotide FASTA"
        ),
    )
    parser.add_argument(
        "--protein-coordinate-mode",
        choices=("candidate", "coding"),
        default="candidate",
        help=(
            "sequence represented in candidate protein FASTA: translate the reported stop-derived "
            "candidate locus by default; coding explicitly uses ranked CodingCoordinates when available"
        ),
    )
    parser.add_argument(
        "--primary-tier",
        choices=("high", "medium", "exploratory"),
        default="exploratory",
        help=(
            "lowest evidence tier included in the primary GFF/FASTA/TSV; lower tiers are written "
            "to separate review files (default: exploratory, include all)"
        ),
    )
    promotion = parser.add_mutually_exclusive_group()
    promotion.add_argument(
        "--promote-cds",
        choices=("high", "medium", "exploratory"),
        default="high",
        metavar="MIN_TIER",
        help="write candidates at/above this tier as CDS only after complete-ORF validation",
    )
    promotion.add_argument(
        "--no-promote-cds",
        dest="promote_cds",
        action="store_const",
        const=None,
        default=argparse.SUPPRESS,
        help="retain every candidate as ORF/StORF even when structural CDS validation passes",
    )
    parser.add_argument(
        "--evidence-tsv",
        help="optional TSV with candidate_id and 0–1 external evidence columns",
    )


def _selected_gene_ident(args: argparse.Namespace) -> str:
    if args.gene_ident is not None:
        return args.gene_ident
    if args.exclude_types:
        return ",".join(args.exclude_types)
    return GENE_IDENT_PROFILES[args.profile]


def _output_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--outdir", type=Path, default=Path("storf-results"))
    parser.add_argument("--gzip", action="store_true")
    parser.add_argument("--overwrite", action="store_true")


def _config(args: argparse.Namespace) -> AnalysisConfig:
    gene_ident = _selected_gene_ident(args)
    selector_mode, selector_values = resolve_gene_ident(gene_ident)
    return AnalysisConfig(
        excluded_types=selector_values if selector_mode == "type" else (),
        gene_ident=gene_ident,
        min_ur=args.min_ur,
        max_ur=args.max_ur,
        extension=args.extension,
        min_orf=args.min_orf,
        max_orf=args.max_orf,
        genetic_code=args.genetic_code,
        ambiguous_threshold=args.ambiguous_threshold,
        include_short=args.include_short or args.short_only,
        short_only=args.short_only,
        short_overlap_mode=args.short_overlap,
        short_min=args.short_min,
        short_max=args.short_max,
        require_start=args.require_start,
        include_consecutive=args.include_consecutive or args.consecutive_only,
        consecutive_only=args.consecutive_only,
        overlap_mode=args.overlap_mode,
        overlap_coordinate_basis=args.overlap_coordinate_basis,
        maximum_overlap=args.maximum_overlap,
        maximum_annotation_overlap=args.maximum_annotation_overlap,
        feature_type=args.feature_type,
        gff_coordinate_mode=args.gff_coordinate_mode,
        include_first_stop=args.include_first_stop,
        protein_coordinate_mode=args.protein_coordinate_mode,
        remove_inframe_stops=args.remove_inframe_stops,
        promote_tier=args.promote_cds,
        primary_tier=args.primary_tier,
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="storf",
        description="Correctness-first discovery and prioritisation of stop-to-stop ORFs",
        epilog=(
            "Workflow guide:\n"
            "  augment      Full annotated-genome workflow. Uses GFF3 plus FASTA or embedded FASTA,\n"
            "               validates each input, extracts unannotated regions (URs), discovers and\n"
            "               validates StORFs, then writes UR files, augmented GFF3, FASTA/TSV/report,\n"
            "               validation JSON, and manifest outputs. --promote-cds adds structurally\n"
            "               valid candidates at the selected evidence tier as CDS; --primary-tier can\n"
            "               keep lower-evidence calls in separate review GFF/FASTA/TSV files.\n"
            "               --include-consecutive enables maximal multi-stop Con-StORFs that compete\n"
            "               with and may replace ordinary StORFs. ORF/StORF GFF coordinates are full\n"
            "               stop-to-stop-derived by default with the upstream stop omitted;\n"
            "               --include-first-stop restores the literal full interval and\n"
            "               --gff-coordinate-mode selects ranked start-to-stop coordinates. Protein\n"
            "               FASTA also uses the candidate locus unless --protein-coordinate-mode coding\n"
            "               is explicitly selected; CDS promotion changes neither FASTA representation.\n"
            "  find         FASTA-only discovery. Treats each input sequence as searchable sequence\n"
            "               because no annotation context is available; use this for raw contigs or\n"
            "               pre-extracted UR FASTA, not for annotation-aware augmentation.\n"
            "  extract-urs  UR extraction only. Writes the unannotated regions implied by an\n"
            "               existing GFF3/FASTA pair without calling StORF discovery.\n"
            "  extract      Sequence extraction only. Pulls already-reported StORF-Reporter features\n"
            "               from a GFF3 into nucleotide or amino-acid FASTA.\n"
            "  filter       Independent homology-evidence step. Reads BLAST/DIAMOND format-6 output,\n"
            "               applies score/e-value/identity/length thresholds, and can promote supported,\n"
            "               structurally valid candidates to CDS while retaining other calls as ORFs.\n"
            "  validate     Input/output validation. Checks GFF3/FASTA agreement, coordinate bounds,\n"
            "               duplicate IDs, StORF IDs, CDS phase, sequence alphabet, and start/stop/\n"
            "               internal-stop structure for StORF-Reporter CDS features.\n"
            "  install-check\n"
            "               Installation/runtime health check. Reports the active Python, the Pyrodigal\n"
            "               dependency, and installed StORF-Reporter command wrappers.\n"
            "\n"
            "In short: use 'augment' for the usual GFF3+FASTA analysis, 'find' for FASTA-only\n"
            "candidate discovery, 'extract-urs' for URs only, 'extract' for sequences from\n"
            "an existing StORF-Reporter GFF3, and 'validate' before trusting or sharing results.\n"
            "The installed 'storf' command is just this unified CLI wrapper; it is not a\n"
            "biological StORF record."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=StORF_Reporter_Version)
    subparsers = parser.add_subparsers(dest="command", required=True)

    augment = subparsers.add_parser(
        "augment",
        help="full GFF3+FASTA workflow: extract URs, find StORFs, and augment annotations",
        description=(
            "Run the full annotation-aware StORF-Reporter workflow. This command loads an existing "
            "GFF3 annotation and its separate or embedded FASTA sequence, validates the pair, identifies "
            "and writes unannotated regions, searches those regions for StORF candidates, validates "
            "candidate coding structure, and validates the final augmented GFF3. A directory processes "
            "multiple unambiguous GFF3/FASTA pairs. Use --promote-cds to produce a mixed final annotation "
            "in which qualifying complete candidates are CDS and all other candidates remain ORF "
            "features. Use --primary-tier to restrict the primary GFF and matching "
            "candidate FASTA/TSV by evidence tier; lower-tier candidates are preserved in separate "
            "review files rather than discarded. Candidate features are inserted by contig and coordinate."
        ),
    )
    inputs = augment.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--input", type=Path, help="GFF3 file or directory")
    inputs.add_argument("--manifest", type=Path, help="TSV with sample, gff and optional fasta")
    augment.add_argument("--fasta", type=Path, help="FASTA paired with a single --input GFF3")
    augment.add_argument("--recursive", action="store_true")
    augment.add_argument("--dry-run", action="store_true")
    augment.add_argument("--embed-fasta", action=argparse.BooleanOptionalAction, default=None)
    _analysis_arguments(augment)
    _output_arguments(augment)

    find = subparsers.add_parser(
        "find",
        help="FASTA-only StORF discovery with no annotation/UR context",
        description=(
            "Search ordinary FASTA sequence for StORF candidates. Because no GFF3 is supplied, each "
            "input sequence is treated as searchable sequence rather than deriving URs from annotation. "
            "Use 'augment' when you want annotation-aware UR extraction and candidate insertion."
        ),
    )
    find.add_argument("--fasta", type=Path, required=True)
    find.add_argument("--sample")
    _analysis_arguments(find)
    _output_arguments(find)

    urs = subparsers.add_parser(
        "extract-urs",
        help="extract unannotated regions only; do not discover StORFs",
        description=(
            "Extract the unannotated regions (URs) implied by a GFF3/FASTA pair. This writes UR GFF3 "
            "and FASTA outputs only; it does not run StORF candidate discovery."
        ),
    )
    urs.add_argument("--gff", type=Path, required=True)
    urs.add_argument("--fasta", type=Path)
    _gene_ident_arguments(urs)
    urs.add_argument("--min-ur", type=int, default=30)
    urs.add_argument("--max-ur", type=int, default=100_000)
    urs.add_argument("--extension", type=int, default=50)
    _output_arguments(urs)

    extract = subparsers.add_parser(
        "extract",
        help="extract existing StORF-Reporter feature sequences from GFF3",
        description=(
            "Extract sequences for StORF-Reporter features that are already present in a GFF3 file. "
            "This is a sequence-retrieval step, not a discovery step."
        ),
    )
    extract.add_argument("--gff", type=Path, required=True)
    extract.add_argument("--fasta", type=Path)
    extract.add_argument("--out", type=Path, required=True)
    extract.add_argument("--amino-acid", action="store_true")
    extract.add_argument(
        "--sequence-region",
        choices=("auto", "candidate", "coding", "feature"),
        default="auto",
        help=(
            "coordinates to extract: auto uses the stop-derived candidate locus for nucleotide and "
            "protein output; coding must be selected explicitly to use CodingCoordinates; feature "
            "uses displayed GFF columns"
        ),
    )
    extract.add_argument(
        "--include-first-stop",
        action="store_true",
        help="retain the upstream delimiting stop in candidate-locus nucleotide or protein output",
    )
    extract.add_argument(
        "--remove-inframe-stops",
        action="store_true",
        help=(
            "remove '*' characters produced by internal in-frame stops from protein output; "
            "nucleotide sequence is never edited"
        ),
    )
    extract.add_argument("--genetic-code", type=int, choices=supported_codes(), default=11)
    extract.add_argument("--overwrite", action="store_true")

    filt = subparsers.add_parser(
        "filter",
        help="filter or promote StORFs using BLAST/DIAMOND format-6 evidence",
        description=(
            "Filter already-reported StORF features using BLAST or DIAMOND tabular format-6 evidence. "
            "Hits must pass the configured bit-score, e-value, identity, and alignment-length thresholds. "
            "By default only supported StORFs are retained. --keep-unmatched preserves unsupported calls "
            "at their existing type, and --promote-cds converts supported calls to CDS only after "
            "sequence-level complete-ORF validation."
        ),
    )
    filt.add_argument("--gff", type=Path, required=True)
    filt.add_argument(
        "--fasta",
        type=Path,
        help="matching FASTA; required for CDS promotion unless the GFF3 embeds FASTA",
    )
    filt.add_argument("--alignments", type=Path, required=True)
    filt.add_argument("--min-bitscore", type=float, default=30.0)
    filt.add_argument("--max-evalue", type=float, default=1e-5)
    filt.add_argument("--min-identity", type=float, default=0.0, help="minimum percent identity")
    filt.add_argument("--min-alignment-length", type=int, default=0, help="minimum aligned residues")
    filt.add_argument(
        "--keep-unmatched",
        action="store_true",
        help="keep StORFs without qualifying hits at their existing feature type",
    )
    filt.add_argument(
        "--promote-cds",
        action="store_true",
        help="promote qualifying, structurally validated StORFs to CDS",
    )
    filt.add_argument("--out", type=Path, required=True)
    filt.add_argument("--overwrite", action="store_true")

    validate = subparsers.add_parser(
        "validate",
        help="validate GFF3/FASTA agreement and StORF CDS structure",
        description=(
            "Validate parsing, sequence-ID agreement, coordinate bounds, duplicate IDs, StORF IDs, "
            "CDS phase, FASTA sequence alphabet, and complete-ORF structure for StORF-Reporter CDS "
            "features. Generic pre-existing CDS features may legitimately be partial or multi-part and "
            "are therefore not required to contain a start and terminal stop."
        ),
    )
    validate.add_argument("--gff", type=Path, required=True)
    validate.add_argument("--fasta", type=Path)
    validate.add_argument("--json", action="store_true")

    subparsers.add_parser(
        "install-check",
        help="report runtime, Pyrodigal dependency, and installed-command health",
        description=(
            "Report the active StORF-Reporter version, Python executable, Pyrodigal dependency, "
            "and installed command wrappers."
        ),
    )
    subparsers.add_parser("citation", help="print the preferred citation")
    return parser


def _pairs(args: argparse.Namespace) -> list[tuple[str, Path, Path | None]]:
    if args.manifest:
        return read_sample_manifest(args.manifest)
    return discover_pairs(args.input, args.fasta, recursive=args.recursive)


def _write_complete_result(
    result: AnalysisResult,
    config: AnalysisConfig,
    outdir: Path,
    gzip_output: bool,
    overwrite: bool,
    embed_fasta: bool | None,
    input_validation,
    include_region_outputs: bool = True,
) -> dict[str, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    prefix = outdir / result.sample
    primary_candidates, review_candidates = partition_candidates(result.candidates, config.primary_tier)
    suffix = ".gff.gz" if gzip_output else ".gff"
    augmented = Path(f"{prefix}_StORF-Reporter{suffix}")
    data_suffix = ".gz" if gzip_output else ""
    report_path = Path(f"{prefix}_report.html")
    manifest_path = Path(f"{prefix}_manifest.json")
    validation_path = Path(f"{prefix}_validation.json")
    expected = [
        augmented,
        Path(f"{prefix}_StORFs.fasta{data_suffix}"),
        Path(f"{prefix}_StORFs_aa.fasta{data_suffix}"),
        Path(f"{prefix}_StORFs.tsv{data_suffix}"),
        report_path,
        validation_path,
        manifest_path,
    ]
    review_gff = Path(f"{prefix}_StORFs_below-{config.primary_tier}{suffix}")
    review_name_suffix = f"_below-{config.primary_tier}"
    if review_candidates:
        expected.extend(
            (
                review_gff,
                Path(f"{prefix}_StORFs{review_name_suffix}.fasta{data_suffix}"),
                Path(f"{prefix}_StORFs{review_name_suffix}_aa.fasta{data_suffix}"),
                Path(f"{prefix}_StORFs{review_name_suffix}.tsv{data_suffix}"),
            )
        )
    if include_region_outputs:
        expected.extend((Path(f"{prefix}_UR.fasta{data_suffix}"), Path(f"{prefix}_UR.gff{data_suffix}")))
    ensure_outputs_available(expected, overwrite)
    write_augmented_gff(
        result,
        config,
        augmented,
        overwrite,
        embed_fasta,
        candidates=primary_candidates,
    )
    outputs = {"augmented_gff": augmented}
    if include_region_outputs:
        region_fasta, region_gff = write_region_outputs(
            result,
            prefix,
            gzip_output,
            overwrite,
            gene_ident=config.gene_ident,
        )
        outputs.update({"ur_fasta": region_fasta, "ur_gff": region_gff})
    outputs.update(
        write_candidate_outputs(
            result,
            prefix,
            gzip_output,
            overwrite,
            candidates=primary_candidates,
            include_leading_stop=config.include_first_stop,
            protein_coordinate_mode=config.protein_coordinate_mode,
            remove_inframe_stops=config.remove_inframe_stops,
        )
    )
    if review_candidates:
        outputs["below_tier_gff"] = write_candidate_gff(
            result,
            config,
            review_gff,
            review_candidates,
            overwrite,
        )
        review_outputs = write_candidate_outputs(
            result,
            prefix,
            gzip_output,
            overwrite,
            candidates=review_candidates,
            name_suffix=review_name_suffix,
            include_leading_stop=config.include_first_stop,
            protein_coordinate_mode=config.protein_coordinate_mode,
            remove_inframe_stops=config.remove_inframe_stops,
        )
        outputs.update({f"below_tier_{name}": path for name, path in review_outputs.items()})
    outputs["html_report"] = write_html_report(
        result,
        report_path,
        overwrite,
        primary_tier=config.primary_tier,
    )
    validation_fasta = result.input_fasta
    if (
        validation_fasta is None
        and result.document is not None
        and result.document.embedded_fasta is not None
    ):
        validation_fasta = result.input_gff
    output_validation = validate_inputs(augmented, validation_fasta)
    review_validation = validate_inputs(review_gff, validation_fasta) if review_candidates else []
    reported_features = [candidate_feature(candidate, config) for candidate in result.candidates]
    candidate_validation = {
        "total": len(reported_features),
        "primary_output": len(primary_candidates),
        "below_primary_tier_output": len(review_candidates),
        "primary_minimum_tier": config.primary_tier,
        "structurally_valid_for_cds": sum(
            feature.attributes.get("CDSValidation") == "passed" for feature in reported_features
        ),
        "not_structurally_valid_for_cds": sum(
            feature.attributes.get("CDSValidation") == "failed" for feature in reported_features
        ),
        "promoted_to_cds": sum(feature.type == "CDS" for feature in reported_features),
        "candidate_validation_passed": sum(
            feature.attributes.get("CandidateValidation") == "passed" for feature in reported_features
        ),
        "candidate_validation_failed": sum(
            feature.attributes.get("CandidateValidation") == "failed" for feature in reported_features
        ),
        "con_storf": {
            "total": sum(
                feature.attributes.get("StORF_Type") == "Con-StORF" for feature in reported_features
            ),
            "valid_readthrough_candidate": sum(
                feature.attributes.get("StORF_Type") == "Con-StORF"
                and feature.attributes.get("CandidateValidation") == "passed"
                for feature in reported_features
            ),
            "invalid_readthrough_candidate": sum(
                feature.attributes.get("StORF_Type") == "Con-StORF"
                and feature.attributes.get("CandidateValidation") == "failed"
                for feature in reported_features
            ),
            "promoted_to_cds": sum(
                feature.attributes.get("StORF_Type") == "Con-StORF" and feature.type == "CDS"
                for feature in reported_features
            ),
        },
    }
    validation_payload = {
        "status": (
            "failed"
            if any(
                issue.severity == "error"
                for issue in (*input_validation, *output_validation, *review_validation)
            )
            else "passed"
        ),
        "scientific_scope": (
            "Structural validation checks format, coordinates, sequence agreement, complete-ORF "
            "properties for CDS candidates, and Con-StORF start-before-middle-stop structure; "
            "it does not prove translation, readthrough, or function."
        ),
        "candidate_cds_validation": candidate_validation,
        "input": [issue.as_dict() for issue in input_validation],
        "output": [issue.as_dict() for issue in output_validation],
        "primary_output": [issue.as_dict() for issue in output_validation],
        "below_primary_tier_output": [issue.as_dict() for issue in review_validation],
    }
    with atomic_output(validation_path, overwrite) as handle:
        json.dump(validation_payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    outputs["validation"] = validation_path
    output_errors = [issue for issue in (*output_validation, *review_validation) if issue.severity == "error"]
    if output_errors:
        preview = "; ".join(f"[{issue.code}] {issue.message}" for issue in output_errors[:3])
        raise InputFormatError(f"generated GFF3 failed post-run validation: {preview}")
    outputs["manifest"] = manifest_path
    write_manifest(result, config, manifest_path, outputs, overwrite)
    return outputs


def _validate_augment_input(sample: str, gff: Path, fasta: Path | None):
    issues = validate_inputs(gff, fasta)
    warnings = [issue for issue in issues if issue.severity == "warning"]
    for issue in warnings:
        print(f"{sample}: WARNING [{issue.code}] {issue.message}", file=sys.stderr)
    errors = [issue for issue in issues if issue.severity == "error"]
    if errors:
        preview = "; ".join(f"[{issue.code}] {issue.message}" for issue in errors[:3])
        raise InputFormatError(f"{sample}: input validation failed: {preview}")
    return issues


def _print_summary(result: AnalysisResult, outputs: dict[str, Path] | None = None) -> None:
    summary = result.summary
    print(
        f"{result.sample}: {summary.sequences} sequences, {summary.unannotated_regions} URs, "
        f"{summary.candidates} candidates "
        f"(high={summary.tiers['high']}, medium={summary.tiers['medium']}, "
        f"exploratory={summary.tiers['exploratory']})"
    )
    if outputs:
        for name, path in outputs.items():
            print(f"  {name}: {path}")


def _augment(args: argparse.Namespace) -> int:
    pairs = _pairs(args)
    if args.dry_run:
        for sample, gff, fasta in pairs:
            print(f"{sample}\t{gff}\t{fasta or '[embedded FASTA]'}")
        return 0
    config = _config(args)
    evidence = read_external_evidence(args.evidence_tsv)
    for sample, gff, fasta in pairs:
        input_validation = _validate_augment_input(sample, gff, fasta)
        document, sequences = load_annotated_genome(gff, fasta)
        result = analyse(
            sequences,
            document.features,
            config,
            external_evidence=evidence,
            document=document,
            input_gff=gff,
            input_fasta=fasta,
            sample=sample,
        )
        outputs = _write_complete_result(
            result,
            config,
            args.outdir / sample,
            args.gzip,
            args.overwrite,
            args.embed_fasta,
            input_validation,
        )
        _print_summary(result, outputs)
    return 0


def _find(args: argparse.Namespace) -> int:
    sequences = read_fasta(args.fasta)
    config = _config(args)
    # Whole sequences are URs when there is no annotation.
    config = AnalysisConfig(
        **{**config.__dict__, "min_ur": 1, "max_ur": max(map(len, sequences.values())), "extension": 0}
    )
    result = analyse(
        sequences,
        [],
        config,
        external_evidence=read_external_evidence(args.evidence_tsv),
        input_fasta=args.fasta.resolve(),
        sample=args.sample or logical_stem(args.fasta),
    )
    outputs = _write_complete_result(
        result,
        config,
        args.outdir,
        args.gzip,
        args.overwrite,
        False,
        [],
        include_region_outputs=False,
    )
    _print_summary(result, outputs)
    return 0


def _extract_urs(args: argparse.Namespace) -> int:
    document, sequences = load_annotated_genome(args.gff, args.fasta)
    gene_ident = _selected_gene_ident(args)
    regions = extract_unannotated_regions(
        sequences,
        document.features,
        min_length=args.min_ur,
        max_length=args.max_ur,
        extension=args.extension,
        gene_ident=gene_ident,
    )
    result = AnalysisResult(
        sequences=sequences,
        features=document.features,
        regions=regions,
        candidates=[],
        summary=RunSummary(
            sequences=len(sequences), annotations=len(document.features), unannotated_regions=len(regions)
        ),
        document=document,
        input_gff=args.gff.resolve(),
        input_fasta=args.fasta.resolve() if args.fasta else None,
        sample=logical_stem(args.gff),
    )
    args.outdir.mkdir(parents=True, exist_ok=True)
    paths = write_region_outputs(
        result,
        args.outdir / logical_stem(args.gff),
        args.gzip,
        args.overwrite,
        gene_ident=gene_ident,
    )
    print(f"{len(regions)} URs written to {paths[0]} and {paths[1]}")
    return 0


def _extract(args: argparse.Namespace) -> int:
    document, sequences = load_annotated_genome(args.gff, args.fasta)
    if args.out.exists() and not args.overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {args.out}")
    records: list[tuple[str, str]] = []
    for feature in document.features:
        if feature.source != "StORF-Reporter" and not str(feature.attributes.get("ID", "")).startswith(
            "StORF_"
        ):
            continue
        region = args.sequence_region
        if region == "auto":
            region = "candidate"
        sequence = extract_feature_sequence(feature, sequences[feature.seqid], region)
        if args.amino_acid:
            if region == "candidate" and not args.include_first_stop:
                sequence = sequence[3:]
            genetic_code = int(feature.attributes.get("GeneticCode", args.genetic_code))
            sequence = translate(
                sequence,
                genetic_code,
                keep_terminal_stop=False,
                start_as_methionine=region == "coding",
            )
            if args.remove_inframe_stops:
                sequence = sequence.replace("*", "")
        elif region == "candidate" and not args.include_first_stop:
            sequence = sequence[3:]
        identifier = feature.attributes.get(
            "ID", f"{feature.seqid}_{feature.interval.start + 1}_{feature.interval.end}"
        )
        records.append((identifier, sequence))
    if not records:
        raise InputFormatError("no StORF-Reporter features found in the GFF3")
    with atomic_output(args.out, args.overwrite) as handle:
        write_fasta_records(handle, records)
    print(f"{len(records)} sequences written to {args.out}")
    return 0


def _coding_interval_from_feature(feature: Feature) -> Interval:
    text = str(feature.attributes.get("CodingCoordinates", ""))
    if (not text or text.casefold() == "none") and feature.type == "CDS":
        return feature.interval
    try:
        start_text, end_text = text.split("-", 1)
        interval = Interval.from_gff(int(start_text), int(end_text))
    except (TypeError, ValueError) as error:
        raise InputFormatError(
            f"{feature.attributes.get('ID', 'StORF feature')} has no valid CodingCoordinates; "
            "rerun augment with v2 or omit --promote-cds"
        ) from error
    if interval.start < feature.interval.start or interval.end > feature.interval.end:
        raise InputFormatError(
            f"{feature.attributes.get('ID', 'StORF feature')} has CodingCoordinates outside its feature"
        )
    return interval


def _homology_annotated_feature(
    feature: Feature,
    hit: dict[str, object],
    promote: bool,
    sequences,
) -> tuple[Feature, bool]:
    attributes: dict[str, object] = dict(feature.attributes)
    attributes.update(
        {
            "HomologyEvidence": "passed",
            "HomologySubject": hit["subject"],
            "HomologyIdentity": f"{hit['identity']:.3f}",
            "HomologyAlignmentLength": hit["alignment_length"],
            "HomologyEvalue": f"{hit['evalue']:.3g}",
            "HomologyBitscore": f"{hit['bitscore']:.3f}",
        }
    )
    if not promote:
        return Feature(
            feature.seqid,
            feature.source,
            feature.type,
            feature.interval,
            feature.score,
            feature.strand,
            feature.phase,
            attributes,
        ), False

    failure: str | None = None
    coding_interval = feature.interval
    try:
        if str(feature.attributes.get("StORF_Type", "")) == "Con-StORF":
            raise InputFormatError(
                "Con-StORFs are multi-stop readthrough candidates and cannot be promoted to a "
                "conventional CDS from homology alone"
            )
        coding_interval = _coding_interval_from_feature(feature)
        genetic_code = int(feature.attributes.get("GeneticCode", 11))
        sequence = sequences[feature.seqid][coding_interval.start : coding_interval.end]
        if feature.strand == "-":
            sequence = reverse_complement(sequence)
        structural = validate_coding_sequence(sequence, genetic_code)
        if not structural.passed:
            failure = " | ".join(structural.issues)
        elif str(feature.attributes.get("CDSValidation", "passed")) == "failed":
            failure = str(feature.attributes.get("CDSValidationIssues", "saved CDS validation failed"))
    except (KeyError, TypeError, ValueError, InputFormatError) as error:
        failure = str(error)

    if failure is not None:
        attributes["HomologyCDSPromotion"] = "failed"
        attributes["HomologyCDSPromotionIssues"] = failure
        return Feature(
            feature.seqid,
            feature.source,
            feature.type,
            feature.interval,
            feature.score,
            feature.strand,
            feature.phase,
            attributes,
        ), False

    attributes.update(
        {
            "CDSValidation": "passed",
            "CDSValidationPolicy": "complete-orf-v1",
            "CDSValidationIssues": "none",
            "HomologyCDSPromotion": "passed",
            "PromotionEvidence": "homology-format6",
            "PromotedToCDS": "true",
        }
    )
    return Feature(
        feature.seqid,
        feature.source,
        "CDS",
        feature.interval,
        feature.score,
        feature.strand,
        "0",
        attributes,
    ), True


def _filter(args: argparse.Namespace) -> int:
    if args.min_bitscore < 0 or args.max_evalue < 0:
        raise ValueError("alignment score and e-value thresholds must be non-negative")
    if not 0 <= args.min_identity <= 100:
        raise ValueError("minimum identity must be between 0 and 100")
    if args.min_alignment_length < 0:
        raise ValueError("minimum alignment length must be non-negative")

    hits: dict[str, dict[str, object]] = {}
    with open_text(args.alignments) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise InputFormatError(f"alignment row {line_number} has fewer than 12 columns")
            try:
                identity = float(fields[2])
                alignment_length = int(fields[3])
                evalue = float(fields[10])
                bitscore = float(fields[11])
            except ValueError as error:
                raise InputFormatError(f"invalid numeric value at alignment row {line_number}") from error
            if not all(math.isfinite(value) for value in (identity, evalue, bitscore)):
                raise InputFormatError(f"non-finite numeric value at alignment row {line_number}")
            if not 0 <= identity <= 100 or alignment_length < 0 or evalue < 0:
                raise InputFormatError(f"out-of-range numeric value at alignment row {line_number}")
            if (
                bitscore >= args.min_bitscore
                and evalue <= args.max_evalue
                and identity >= args.min_identity
                and alignment_length >= args.min_alignment_length
            ):
                hit = {
                    "subject": fields[1],
                    "identity": identity,
                    "alignment_length": alignment_length,
                    "evalue": evalue,
                    "bitscore": bitscore,
                }
                query_id = storf_query_id(fields[0])
                previous = hits.get(query_id)
                if previous is None or (bitscore, -evalue) > (
                    previous["bitscore"],
                    -previous["evalue"],
                ):
                    hits[query_id] = hit
    if args.out.exists() and not args.overwrite:
        raise FileExistsError(f"output exists; use --overwrite: {args.out}")
    document = read_gff(args.gff)
    sequences = None
    if args.promote_cds:
        sequences = read_fasta(args.fasta) if args.fasta else document.embedded_fasta
        if sequences is None:
            raise InputFormatError(
                "--promote-cds requires --fasta unless the input GFF3 contains embedded FASTA"
            )
        preflight = validate_inputs(args.gff, args.fasta)
        errors = [issue for issue in preflight if issue.severity == "error"]
        if errors:
            preview = "; ".join(f"[{issue.code}] {issue.message}" for issue in errors[:3])
            raise InputFormatError(f"filter input validation failed: {preview}")

    kept = removed = promoted = promotion_failed = 0
    seen_candidate_ids: set[str] = set()
    feature_number = 0
    with atomic_output(args.out, args.overwrite) as output:
        for raw in document.pre_fasta_lines:
            if raw.startswith("#") or not raw.strip():
                output.write(raw)
                continue
            feature = document.features[feature_number]
            feature_number += 1
            identifier = str(feature.attributes.get("ID", ""))
            is_storf = feature.source == "StORF-Reporter" or identifier.startswith("StORF_")
            if not is_storf:
                output.write(raw)
                continue
            seen_candidate_ids.add(identifier)
            hit = hits.get(identifier)
            if hit is not None:
                assert sequences is not None or not args.promote_cds
                annotated, did_promote = _homology_annotated_feature(
                    feature, hit, args.promote_cds, sequences
                )
                output.write(feature_line(annotated) + "\n")
                kept += 1
                promoted += int(did_promote)
                promotion_failed += int(args.promote_cds and not did_promote and feature.type != "CDS")
            elif args.keep_unmatched:
                output.write(raw)
                kept += 1
            else:
                removed += 1
        if document.embedded_fasta:
            output.write("##FASTA\n")
            write_fasta_records(output, document.embedded_fasta.items())

    if args.promote_cds:
        postflight = validate_inputs(args.out, args.fasta)
        errors = [issue for issue in postflight if issue.severity == "error"]
        if errors:
            preview = "; ".join(f"[{issue.code}] {issue.message}" for issue in errors[:3])
            raise InputFormatError(f"filtered GFF3 failed post-run validation: {preview}")
    unmatched_hits = len(set(hits) - seen_candidate_ids)
    print(
        f"kept {kept} StORFs; removed {removed}; promoted {promoted} to CDS; "
        f"promotion failed for {promotion_failed}; qualifying query IDs absent from GFF3: {unmatched_hits}; "
        f"output: {args.out}"
    )
    return 0


def _validate(args: argparse.Namespace) -> int:
    issues = validate_inputs(args.gff, args.fasta)
    if args.json:
        print(json.dumps([issue.as_dict() for issue in issues], indent=2))
    elif issues:
        for issue in issues:
            print(f"{issue.severity.upper()} [{issue.code}] {issue.message}")
    else:
        print("valid: no GFF3/FASTA issues found")
    return 1 if any(issue.severity == "error" for issue in issues) else 0


def _install_check() -> int:
    print(f"StORF-Reporter: {StORF_Reporter_Version}")
    print(f"Python: {platform.python_version()} ({sys.executable})")
    try:
        import pyrodigal

        print(f"Pyrodigal: {pyrodigal.__version__}")
    except ImportError:
        print("Pyrodigal: not installed")
    scripts_dir = Path(sys.executable).parent
    print("Command wrappers:")
    for command, label in (
        ("storf", "unified CLI"),
        ("StORF-Reporter", "legacy full workflow"),
        ("StORF-Finder", "legacy FASTA/UR candidate finder"),
        ("UR-Extractor", "legacy unannotated-region extractor"),
        ("StORF-Extractor", "legacy reported-feature sequence extractor"),
        ("StORF-Remover", "legacy alignment-evidence filter"),
        ("storf-add", "merge filtered candidates into a reference GFF3"),
        ("storf-uncovered", "extract sequence uncovered after augmentation"),
        ("storf-lengths", "summarise FASTA sequence lengths"),
    ):
        local_command = scripts_dir / command
        resolved = str(local_command) if local_command.exists() else shutil.which(command)
        print(f"  {command} ({label}): {resolved or 'not installed'}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "augment":
            return _augment(args)
        if args.command == "find":
            return _find(args)
        if args.command == "extract-urs":
            return _extract_urs(args)
        if args.command == "extract":
            return _extract(args)
        if args.command == "filter":
            return _filter(args)
        if args.command == "validate":
            return _validate(args)
        if args.command == "install-check":
            return _install_check()
        if args.command == "citation":
            print(CITATION)
            return 0
    except (FileExistsError, FileNotFoundError, InputFormatError, OSError, RuntimeError, ValueError) as error:
        print(f"storf: error: {error}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
