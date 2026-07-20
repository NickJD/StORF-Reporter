"""Legacy ``UR-Extractor`` compatibility command backed by the v2 core."""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .core import extract_unannotated_regions, feature_matches_gene_ident
from .io import InputFormatError
from .models import RunSummary
from .utils import parse_bool as _bool
from .workflow import AnalysisResult, load_annotated_genome, logical_stem, write_region_outputs


def _result(options: argparse.Namespace) -> AnalysisResult:
    document, sequences = load_annotated_genome(options.gff, getattr(options, "fasta", None))
    regions = extract_unannotated_regions(
        sequences,
        document.features,
        min_length=int(options.minlen),
        max_length=int(options.maxlen),
        extension=int(options.exlen),
        gene_ident=options.gene_ident,
    )
    return AnalysisResult(
        sequences=sequences,
        features=document.features,
        regions=regions,
        candidates=[],
        summary=RunSummary(
            sequences=len(sequences),
            annotations=len(document.features),
            unannotated_regions=len(regions),
        ),
        document=document,
        input_gff=Path(options.gff).resolve(),
        input_fasta=Path(options.fasta).resolve() if getattr(options, "fasta", None) else None,
        sample=logical_stem(options.gff),
    )


def extractor(options: argparse.Namespace):
    """Compatibility API used by older StORF-Reporter releases.

    When ``options.nout`` is true, return the historical dictionary structure;
    otherwise write standards-compliant UR FASTA and GFF3 outputs.
    """

    result = _result(options)
    if bool(getattr(options, "nout", False)):
        legacy = collections.OrderedDict()
        by_seqid: dict[str, collections.OrderedDict[str, list[str]]] = {}
        for seqid in result.sequences:
            by_seqid[seqid] = collections.OrderedDict()
        for region in result.regions:
            ext_start, ext_end = region.extended.to_gff()
            core_start, core_end = region.core.to_gff()
            by_seqid[region.seqid][f"{ext_start}_{ext_end}"] = [
                f"{core_start}_{core_end}",
                region.sequence,
            ]
        positions: dict[str, list[str]] = {seqid: [] for seqid in result.sequences}
        for feature in result.features:
            if feature_matches_gene_ident(feature, options.gene_ident):
                positions[feature.seqid].append(f"{feature.interval.start + 1}_{feature.interval.end}")
        for seqid, sequence in result.sequences.items():
            legacy[seqid] = (sequence, len(sequence), sorted(positions[seqid]), by_seqid[seqid])
        return legacy

    prefix = getattr(options, "output_prefix", None)
    if prefix is None:
        output_dir = Path(getattr(options, "o_dir", None) or Path(options.gff).parent)
        name = getattr(options, "o_name", None) or logical_stem(options.gff)
        prefix = output_dir / name
    paths = write_region_outputs(
        result,
        prefix,
        gzip_output=bool(getattr(options, "gz", False)),
        overwrite=bool(getattr(options, "overwrite", False)),
        gene_ident=str(options.gene_ident),
    )
    return paths


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="UR-Extractor",
        description=f"StORF-Reporter {StORF_Reporter_Version}: extract unannotated regions",
    )
    parser.add_argument("-gff", required=True, help="GFF3 genome annotation")
    parser.add_argument("-f", dest="fasta", help="matching FASTA; optional for embedded FASTA GFF3")
    parser.add_argument("-ident", default="_UR", help=argparse.SUPPRESS)
    parser.add_argument("-min_len", dest="minlen", type=int, default=30)
    parser.add_argument("-max_len", dest="maxlen", type=int, default=100_000)
    parser.add_argument("-ex_len", dest="exlen", type=int, default=50)
    parser.add_argument(
        "-gene_ident",
        dest="gene_ident",
        default="CDS",
        help=(
            "feature types that occupy annotated regions: CDS by default; provide an exact "
            "comma-separated list, Prokka/Bakta, or ID=gene for attribute-prefix matching"
        ),
    )
    parser.add_argument("-oname", dest="o_name")
    parser.add_argument("-odir", dest="o_dir")
    parser.add_argument("-gz", dest="gz", type=_bool, default=False)
    parser.add_argument("-verbose", dest="verbose", type=_bool, default=False)
    parser.add_argument("-nout", dest="nout", type=_bool, default=False, help=argparse.SUPPRESS)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("-v", "--version", action="version", version=StORF_Reporter_Version)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    try:
        paths = extractor(options)
        if options.verbose and paths:
            print(f"UR FASTA: {paths[0]}\nUR GFF3: {paths[1]}")
        return 0
    except (FileExistsError, FileNotFoundError, InputFormatError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
