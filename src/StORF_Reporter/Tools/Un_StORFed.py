"""Extract sequence still unannotated after StORF augmentation."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    __package__ = "StORF_Reporter.Tools"

from ..core import extract_unannotated_regions, feature_matches_gene_ident
from ..io import read_fasta, read_gff
from ..models import RunSummary
from ..utils import parse_bool as _bool
from ..workflow import AnalysisResult, write_region_outputs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Extract regions uncovered by original annotations and StORFs"
    )
    parser.add_argument("-f", "--fasta_seq", dest="fasta", required=True, help="whole-genome FASTA")
    parser.add_argument("-gff", "--gff", required=True, help="augmented GFF3")
    parser.add_argument("-min_len", "--min-length", dest="minlen", default=50, type=int)
    parser.add_argument(
        "-gene_ident",
        "--gene-ident",
        dest="gene_ident",
        default="CDS",
        help=(
            "original annotations that occupy sequence: exact comma-separated types, "
            "Prokka/Bakta, or an attribute query such as ID=gene (default: CDS)"
        ),
    )
    parser.add_argument("-o", "--output_file", "--output-prefix", dest="out_file")
    parser.add_argument("-gz", "--gzip", dest="gz", nargs="?", const=True, default=False, type=_bool)
    parser.add_argument("--overwrite", action="store_true")
    options = parser.parse_args(argv)

    sequences = read_fasta(options.fasta)
    document = read_gff(options.gff)
    occupying = [
        feature
        for feature in document.features
        if feature_matches_gene_ident(feature, options.gene_ident)
        or feature.source == "StORF-Reporter"
        or str(feature.attributes.get("ID", "")).startswith("StORF_")
    ]
    excluded = {feature.type for feature in occupying}
    regions = extract_unannotated_regions(
        sequences,
        occupying,
        excluded,
        min_length=options.minlen,
        max_length=max(map(len, sequences.values())),
        extension=0,
    )
    result = AnalysisResult(
        sequences=sequences,
        features=document.features,
        regions=regions,
        candidates=[],
        summary=RunSummary(unannotated_regions=len(regions)),
        sample=Path(options.gff).stem,
    )
    prefix = Path(options.out_file) if options.out_file else Path(options.gff).with_suffix("")
    paths = write_region_outputs(result, f"{prefix}_UnStORFed", options.gz, options.overwrite)
    print(f"{len(regions)} uncovered regions: {paths[0]} and {paths[1]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
