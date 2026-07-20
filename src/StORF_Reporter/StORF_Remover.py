"""Legacy ``StORF-Remover`` command for evidence-based GFF3 filtering."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "StORF_Reporter"

from .constants import StORF_Reporter_Version
from .io import InputFormatError, open_text, read_gff, write_fasta_records
from .utils import parse_bool as _bool
from .workflow import atomic_output, logical_stem, storf_query_id


def _load_hits(path: str | Path, minimum: float) -> set[str]:
    hits: set[str] = set()
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise InputFormatError(f"BLAST/DIAMOND row {line_number} has fewer than 12 columns")
            try:
                bitscore = float(fields[11])
            except ValueError as error:
                raise InputFormatError(f"invalid bitscore at row {line_number}: {fields[11]!r}") from error
            if bitscore >= minimum:
                hits.add(storf_query_id(fields[0]))
    return hits


def remover(options: argparse.Namespace) -> tuple[int, int]:
    hits = _load_hits(options.blast, float(options.minscore))
    document = read_gff(options.gff)
    kept = removed = 0
    with atomic_output(options.output, options.overwrite) as handle:
        for raw in document.pre_fasta_lines:
            if raw.startswith("#") or not raw.strip():
                handle.write(raw)
                continue
            fields = raw.rstrip("\n").split("\t")
            identifier = next((item[3:] for item in fields[8].split(";") if item.startswith("ID=")), "")
            is_storf = (
                fields[1] == "StORF-Reporter" or identifier.startswith("StORF_") or "_StORF_" in identifier
            )
            if not is_storf or identifier in hits:
                handle.write(raw)
                kept += int(is_storf)
            else:
                removed += 1
        if document.embedded_fasta:
            handle.write("##FASTA\n")
            write_fasta_records(handle, document.embedded_fasta.items())
    return kept, removed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="StORF-Remover",
        description=f"StORF-Reporter {StORF_Reporter_Version}: retain candidates with alignment evidence",
    )
    parser.add_argument("-gff", required=True)
    parser.add_argument("-blast", required=True, help="BLAST/DIAMOND tabular format 6")
    parser.add_argument("-min_score", dest="minscore", type=float, default=30.0)
    parser.add_argument("-oname", dest="o_name")
    parser.add_argument("-odir", dest="o_dir")
    parser.add_argument("-gz", dest="gz", type=_bool, default=False)
    parser.add_argument("-verbose", dest="verbose", type=_bool, default=False)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("-v", "--version", action="version", version=StORF_Reporter_Version)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    options = parser.parse_args(argv)
    try:
        gff = Path(options.gff).resolve()
        output_dir = Path(options.o_dir).resolve() if options.o_dir else gff.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = options.o_name or f"{logical_stem(gff)}_NoHit_StORFs_Removed"
        options.output = output_dir / f"{stem}.gff{'.gz' if options.gz else ''}"
        kept, removed = remover(options)
        if options.verbose:
            print(f"kept {kept}; removed {removed}; output: {options.output}")
        return 0
    except (FileExistsError, FileNotFoundError, InputFormatError, OSError, ValueError) as error:
        parser.error(str(error))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
