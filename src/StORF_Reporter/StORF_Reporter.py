"""Legacy ``StORF-Reporter`` CLI implemented as a v2 workflow adapter.

The historical command remains available for scripts, but all biological
decisions and coordinates now pass through :mod:`StORF_Reporter.core`.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .core import resolve_gene_ident
from .io import InputFormatError, read_fasta
from .models import RunSummary
from .utils import parse_bool as _bool
from .workflow import (
    AnalysisConfig,
    AnalysisResult,
    analyse,
    discover_pairs,
    ensure_outputs_available,
    is_fasta_path,
    load_annotated_genome,
    logical_stem,
    pyrodigal_features,
    write_augmented_gff,
    write_candidate_outputs,
    write_combined_cds_storf_fasta,
    write_html_report,
    write_manifest,
)


def _annotation_pairs(options: argparse.Namespace) -> list[tuple[str, Path, Path | None]]:
    annotation, mode = options.annotation_type
    path = Path(options.path).resolve()
    if annotation == "Pyrodigal":
        return []
    if options.alt_filename:
        if mode != "Out_Dir":
            raise ValueError("-af is only valid with Prokka/Bakta Out_Dir")
        for extension in (".gff3", ".gff"):
            candidate = path / f"{options.alt_filename}{extension}"
            if candidate.exists():
                return discover_pairs(candidate)
        raise FileNotFoundError(f"no GFF3 named {options.alt_filename!r} in {path}")
    recursive = mode == "Multiple_Out_Dirs"
    pairs = discover_pairs(path, recursive=recursive)

    # Bakta can emit equivalent .gff and .gff3 files. Prefer .gff3 once per stem.
    chosen: dict[tuple[Path, str], tuple[str, Path, Path | None]] = {}
    for pair in pairs:
        gff = pair[1]
        key = (gff.parent, logical_stem(gff))
        previous = chosen.get(key)
        preferred = gff.name.lower().endswith((".gff3", ".gff3.gz"))
        previous_preferred = previous and previous[1].name.lower().endswith((".gff3", ".gff3.gz"))
        if previous is None or (preferred and not previous_preferred):
            chosen[key] = pair
    return [chosen[key] for key in sorted(chosen, key=lambda item: (str(item[0]), item[1]))]


def _gene_ident(options: argparse.Namespace) -> str:
    annotation = options.annotation_type[0]
    if annotation in {"Prokka", "Bakta"} or options.gene_ident == "Prokka":
        return annotation if annotation in {"Prokka", "Bakta"} else "Prokka"
    if annotation == "Ensembl" or options.gene_ident == "ID=gene":
        return "ID=gene"
    return options.gene_ident


def _code(options: argparse.Namespace) -> int:
    if options.genetic_code is not None:
        return options.genetic_code
    stops = frozenset(item.strip().upper() for item in options.stop_codons.split(","))
    if stops == frozenset({"TAA", "TAG"}):
        return 4
    if stops == frozenset({"TAA", "TAG", "TGA"}):
        return 11
    raise ValueError("custom stop sets are ambiguous; select --genetic-code 1, 4 or 11")


def _config(options: argparse.Namespace) -> AnalysisConfig:
    if options.partial_storf or options.whole_contig:
        raise ValueError(
            "partial and run-through modes are disabled in the correctness-first release because "
            "they are not stop-to-stop ORFs"
        )
    requested_type = options.feature_type
    promote = None
    if requested_type == "CDS":
        requested_type = "ORF"
        promote = "exploratory"
    gene_ident = _gene_ident(options)
    selector_mode, selector_values = resolve_gene_ident(gene_ident)
    return AnalysisConfig(
        excluded_types=selector_values if selector_mode == "type" else (),
        gene_ident=gene_ident,
        min_ur=options.minlen,
        max_ur=options.maxlen,
        extension=options.exlen,
        min_orf=options.min_orf,
        max_orf=options.max_orf,
        genetic_code=_code(options),
        ambiguous_threshold=float(options.non_standard),
        include_short=options.short_storfs != "False" or options.short_storfs_only,
        short_only=options.short_storfs_only,
        short_overlap_mode="allow" if options.short_storfs == "Olap" else "exclude",
        short_min=30,
        short_max=120,
        require_start=options.start_filtering,
        include_consecutive=options.con_storfs or options.con_only,
        consecutive_only=options.con_only,
        overlap_mode={"none": "none", "single-strand": "same-strand", "both-strand": "all"}[
            options.olap_filtering
        ],
        overlap_coordinate_basis=getattr(options, "overlap_coordinate_basis", "coding"),
        maximum_overlap=options.overlap_nt,
        maximum_annotation_overlap=options.allowed_overlap,
        feature_type=requested_type,
        gff_coordinate_mode=getattr(options, "gff_coordinate_mode", "stop-to-stop"),
        include_first_stop=bool(options.stop_inclusive),
        protein_coordinate_mode=options.protein_coordinate_mode,
        promote_tier=promote,
    )


def _write_result(
    result: AnalysisResult,
    config: AnalysisConfig,
    options: argparse.Namespace,
    default_directory: Path,
    default_stem: str,
    output_number: int | None = None,
) -> dict[str, Path]:
    output_dir = Path(options.o_dir).resolve() if options.o_dir else default_directory
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = options.o_name or default_stem
    if options.o_name and output_number is not None:
        stem = f"{stem}_{output_number + 1}"
    prefix = output_dir / stem
    gzip_suffix = ".gz" if options.gz else ""
    gff_path = Path(f"{prefix}.gff{gzip_suffix}")
    expected = [
        gff_path,
        Path(f"{prefix}_StORFs.fasta{gzip_suffix}"),
        Path(f"{prefix}_StORFs_aa.fasta{gzip_suffix}"),
        Path(f"{prefix}_StORFs.tsv{gzip_suffix}"),
        Path(f"{prefix}_report.html"),
        Path(f"{prefix}_manifest.json"),
    ]
    write_combined = options.annotation_type[1] in {"Out_Dir", "Multiple_Out_Dirs"} or (
        options.annotation_type[0] == "Pyrodigal" and options.py_fasta
    )
    if write_combined:
        expected.append(Path(f"{prefix}.fasta{gzip_suffix}"))
        if options.translate:
            expected.append(Path(f"{prefix}_aa.fasta{gzip_suffix}"))
    ensure_outputs_available(expected, options.overwrite)
    outputs: dict[str, Path] = {
        "augmented_gff": write_augmented_gff(
            result,
            config,
            gff_path,
            overwrite=options.overwrite,
            embed_fasta=None,
        )
    }
    outputs.update(
        write_candidate_outputs(
            result,
            prefix,
            gzip_output=options.gz,
            overwrite=options.overwrite,
            line_width=60 if options.line_wrap else 0,
            include_leading_stop=options.stop_inclusive,
            keep_terminal_stop=not options.remove_stop,
            protein_coordinate_mode=options.protein_coordinate_mode,
            remove_inframe_stops=options.remove_stop,
        )
    )
    if write_combined:
        combined = Path(f"{prefix}.fasta{gzip_suffix}")
        outputs["combined_nucleotide_fasta"] = write_combined_cds_storf_fasta(
            result,
            combined,
            overwrite=options.overwrite,
            line_width=60 if options.line_wrap else 0,
            include_leading_stop=options.stop_inclusive,
        )
        if options.translate:
            combined_aa = Path(f"{prefix}_aa.fasta{gzip_suffix}")
            outputs["combined_protein_fasta"] = write_combined_cds_storf_fasta(
                result,
                combined_aa,
                amino_acid=True,
                genetic_code=config.genetic_code,
                overwrite=options.overwrite,
                line_width=60 if options.line_wrap else 0,
                include_leading_stop=options.stop_inclusive,
                protein_coordinate_mode=options.protein_coordinate_mode,
                remove_inframe_stops=options.remove_stop,
            )
    report = Path(f"{prefix}_report.html")
    outputs["html_report"] = write_html_report(result, report, options.overwrite)
    manifest = Path(f"{prefix}_manifest.json")
    outputs["manifest"] = manifest
    write_manifest(result, config, manifest, outputs, options.overwrite)
    return outputs


def _run_annotated(options: argparse.Namespace, config: AnalysisConfig) -> list[AnalysisResult]:
    results: list[AnalysisResult] = []
    pairs = _annotation_pairs(options)
    for number, (sample, gff, fasta) in enumerate(pairs):
        document, sequences = load_annotated_genome(gff, fasta)
        result = analyse(
            sequences,
            document.features,
            config,
            document=document,
            input_gff=gff,
            input_fasta=fasta,
            sample=sample,
        )
        if options.storf_order == "strand":
            result.candidates.sort(
                key=lambda item: (item.seqid, item.strand == "-", item.interval.start, item.interval.end)
            )
        stem = f"{sample}_StORF-Reporter_Extended"
        outputs = _write_result(
            result,
            config,
            options,
            gff.parent,
            stem,
            number if len(pairs) > 1 else None,
        )
        if options.verbose:
            print(f"{sample}: {result.summary.candidates} candidates")
            for name, output in outputs.items():
                print(f"  {name}: {output}")
        results.append(result)
    return results


def _run_pyrodigal(options: argparse.Namespace, config: AnalysisConfig) -> list[AnalysisResult]:
    path = Path(options.path).resolve()
    if path.is_file():
        fastas = [path]
    elif path.is_dir():
        fastas = sorted(item for item in path.iterdir() if item.is_file() and is_fasta_path(item))
    else:
        raise FileNotFoundError(f"input does not exist: {path}")
    if not fastas:
        raise FileNotFoundError(f"no FASTA inputs found under {path}")
    results = []
    for number, fasta in enumerate(fastas):
        sequences = read_fasta(fasta)
        features = pyrodigal_features(sequences, options.py_train)
        result = analyse(
            sequences,
            features,
            config,
            input_fasta=fasta,
            sample=logical_stem(fasta),
        )
        if options.storf_order == "strand":
            result.candidates.sort(
                key=lambda item: (item.seqid, item.strand == "-", item.interval.start, item.interval.end)
            )
        fasta_stem = logical_stem(fasta)
        default_stem = f"{fasta_stem}_Pyrodigal_StORF-Reporter_Extended"
        outputs = _write_result(
            result,
            config,
            options,
            fasta.parent,
            default_stem,
            number if len(fastas) > 1 else None,
        )
        if options.py_unstorfed:
            original = AnalysisResult(
                sequences=sequences,
                features=features,
                regions=[],
                candidates=[],
                summary=RunSummary(sequences=len(sequences), annotations=len(features)),
                input_fasta=fasta,
                sample=fasta_stem,
            )
            output_dir = Path(options.o_dir).resolve() if options.o_dir else fasta.parent
            pyrodigal_gff = output_dir / f"{fasta_stem}_Pyrodigal_Predictions.gff"
            write_augmented_gff(original, config, pyrodigal_gff, options.overwrite, False)
            outputs["pyrodigal_gff"] = pyrodigal_gff
        if options.verbose:
            print(f"{fasta.name}: {len(features)} Pyrodigal CDS; {result.summary.candidates} candidates")
        results.append(result)
    return results


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="StORF-Reporter",
        description=f"StORF-Reporter {StORF_Reporter_Version}: augment prokaryotic annotations",
    )
    parser.add_argument(
        "-anno",
        dest="annotation_type",
        nargs=2,
        required=True,
        metavar=("ANNOTATION", "INPUT_MODE"),
        choices=(
            "Prokka",
            "Bakta",
            "Out_Dir",
            "Multiple_Out_Dirs",
            "Single_GFF",
            "Multiple_GFFs",
            "Ensembl",
            "Feature_Types",
            "Single_Genome",
            "Multiple_Genomes",
            "Single_Combined_GFF",
            "Multiple_Combined_GFFs",
            "Pyrodigal",
            "Single_FASTA",
            "Multiple_FASTA",
        ),
    )
    parser.add_argument("-p", dest="path", required=True)
    parser.add_argument("-af", dest="alt_filename")
    parser.add_argument("-oname", dest="o_name")
    parser.add_argument("-odir", dest="o_dir")
    parser.add_argument(
        "-sout",
        dest="storfs_out",
        type=_bool,
        default=False,
        help="Default - False: Print out StORF sequences separately",
    )
    parser.add_argument("-lw", dest="line_wrap", type=_bool, default=True)
    parser.add_argument("-aa", dest="translate", action="store_true")
    parser.add_argument("-gz", dest="gz", type=_bool, default=False)
    parser.add_argument("-py_train", choices=("longest", "individual", "meta"), default="longest")
    parser.add_argument("-py_fasta", dest="py_fasta", type=_bool, default=True)
    parser.add_argument("-py_unstorfed", dest="py_unstorfed", type=_bool, default=False)
    parser.add_argument(
        "-gene_ident",
        dest="gene_ident",
        default="CDS",
        help=(
            "feature types that occupy annotated regions: CDS by default; provide an exact "
            "comma-separated list, Prokka/Bakta, or ID=gene for Ensembl attribute matching"
        ),
    )
    parser.add_argument("-min_len", dest="minlen", type=int, default=30)
    parser.add_argument("-max_len", dest="maxlen", type=int, default=100_000)
    parser.add_argument("-ex_len", dest="exlen", type=int, default=50)
    parser.add_argument(
        "-spos",
        dest="stop_inclusive",
        type=_bool,
        default=False,
        help=("include the upstream delimiting stop in GFF3 and candidate-locus FASTA (default: False)"),
    )
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
    parser.add_argument("-con_storfs", dest="con_storfs", type=_bool, default=False)
    parser.add_argument("-con_only", dest="con_only", type=_bool, default=False)
    parser.add_argument("-ps", dest="partial_storf", type=_bool, default=False)
    parser.add_argument("-wc", dest="whole_contig", type=_bool, default=False)
    parser.add_argument("-short_storfs", choices=("False", "Nolap", "Olap"), default="False")
    parser.add_argument("-short_storfs_only", dest="short_storfs_only", type=_bool, default=False)
    parser.add_argument("-minorf", dest="min_orf", type=int, default=99)
    parser.add_argument("-maxorf", dest="max_orf", type=int, default=60_000)
    parser.add_argument("-codons", dest="stop_codons", default="TAG,TGA,TAA")
    parser.add_argument("--genetic-code", type=int, choices=(1, 4, 11))
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
    parser.add_argument("-f_type", dest="feature_type", default="ORF")
    parser.add_argument(
        "--gff-coordinate-mode",
        choices=("stop-to-stop", "start-to-stop"),
        default="stop-to-stop",
        help="GFF coordinates for unpromoted candidates (default: stop-to-stop)",
    )
    parser.add_argument("-non_standard", dest="non_standard", type=float, default=0.20)
    parser.add_argument("-olap", dest="overlap_nt", type=int, default=50)
    parser.add_argument("-ao", dest="allowed_overlap", type=int, default=50)
    parser.add_argument("-so", dest="storf_order", choices=("start_pos", "strand"), default="start_pos")
    parser.add_argument("-overwrite", "--overwrite", action="store_true")
    parser.add_argument("-verbose", "--verbose", action="store_true")
    parser.add_argument("-v", "--version", action="version", version=StORF_Reporter_Version)
    return parser


def _validate_mode(options: argparse.Namespace) -> None:
    annotation, mode = options.annotation_type
    compatible = {
        "Prokka": {"Out_Dir", "Multiple_Out_Dirs", "Single_GFF", "Multiple_GFFs"},
        "Bakta": {"Out_Dir", "Multiple_Out_Dirs", "Single_GFF", "Multiple_GFFs"},
        "Ensembl": {"Single_Genome", "Multiple_Genomes", "Single_Combined_GFF", "Multiple_Combined_GFFs"},
        "Feature_Types": {
            "Single_Genome",
            "Multiple_Genomes",
            "Single_Combined_GFF",
            "Multiple_Combined_GFFs",
        },
        "Pyrodigal": {"Single_FASTA", "Multiple_FASTA"},
    }
    if mode not in compatible.get(annotation, set()):
        raise ValueError(f"incompatible -anno combination: {annotation} {mode}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    try:
        _validate_mode(options)
        config = _config(options)
        results = (
            _run_pyrodigal(options, config)
            if options.annotation_type[0] == "Pyrodigal"
            else _run_annotated(options, config)
        )
        if not options.verbose:
            total = sum(result.summary.candidates for result in results)
            print(f"complete: {len(results)} genome(s), {total} retained StORF candidates")
        return 0
    except (
        FileExistsError,
        FileNotFoundError,
        InputFormatError,
        OSError,
        RuntimeError,
        ValueError,
    ) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
