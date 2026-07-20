"""Strict FASTA/GFF3 input and output helpers."""

from __future__ import annotations

import gzip
from collections import OrderedDict
from collections.abc import Iterable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TextIO
from urllib.parse import quote, unquote

from .models import Feature, Interval


class InputFormatError(ValueError):
    """Raised when an input cannot be interpreted without guessing."""


@contextmanager
def open_text(path: str | Path, mode: str = "rt") -> Iterator[TextIO]:
    """Open plain text or gzip data, detecting gzip by magic bytes on read."""

    path = Path(path)
    if "r" in mode:
        with path.open("rb") as probe:
            compressed = probe.read(2) == b"\x1f\x8b"
        if compressed:
            handle = gzip.open(path, mode, encoding="utf-8", newline="")
        else:
            handle = path.open(mode, encoding="utf-8", newline="")
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".gz":
            handle = gzip.open(path, mode, encoding="utf-8", newline="")
        else:
            handle = path.open(mode, encoding="utf-8", newline="")
    try:
        yield handle
    finally:
        handle.close()


def parse_fasta_lines(lines: Iterable[str], source: str = "FASTA") -> OrderedDict[str, str]:
    sequences: OrderedDict[str, str] = OrderedDict()
    current: str | None = None
    chunks: list[str] = []

    def flush() -> None:
        nonlocal chunks
        if current is None:
            return
        if current in sequences:
            raise InputFormatError(f"duplicate sequence identifier {current!r} in {source}")
        sequence = "".join(chunks).replace(" ", "").upper()
        if not sequence:
            raise InputFormatError(f"sequence {current!r} in {source} is empty")
        invalid = sorted(set(sequence) - set("ACGTRYSWKMBDHVN.-"))
        if invalid:
            raise InputFormatError(
                f"sequence {current!r} in {source} contains non-IUPAC characters: {''.join(invalid)}"
            )
        sequences[current] = sequence
        chunks = []

    for line_number, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            current = line[1:].split()[0]
            if not current:
                raise InputFormatError(f"empty FASTA identifier at {source}:{line_number}")
        elif current is None:
            raise InputFormatError(f"sequence before FASTA header at {source}:{line_number}")
        else:
            chunks.append(line)
    flush()
    if not sequences:
        raise InputFormatError(f"no FASTA sequences found in {source}")
    return sequences


def read_fasta(path: str | Path) -> OrderedDict[str, str]:
    path = Path(path)
    with open_text(path) as handle:
        lines = list(handle)
    fasta_index = next(
        (index + 1 for index, line in enumerate(lines) if line.startswith("##FASTA")),
        None,
    )
    if fasta_index is not None:
        lines = lines[fasta_index:]
    return parse_fasta_lines(lines, str(path))


def parse_attributes(text: str) -> dict[str, str]:
    if text == "." or not text:
        return {}
    attributes: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if not item:
            continue
        if "=" not in item:
            raise InputFormatError(f"invalid GFF3 attribute {item!r}")
        key, value = item.split("=", 1)
        attributes[unquote(key)] = unquote(value)
    return attributes


def format_attributes(attributes: dict[str, object] | None) -> str:
    if not attributes:
        return "."
    safe = ":@|+-_."
    return ";".join(
        f"{quote(str(key), safe=safe)}={quote(str(value), safe=safe + ',')}"
        for key, value in attributes.items()
        if value is not None
    )


class GffDocument:
    def __init__(
        self,
        path: Path,
        pre_fasta_lines: list[str],
        features: list[Feature],
        embedded_fasta: OrderedDict[str, str] | None,
    ) -> None:
        self.path = path
        self.pre_fasta_lines = pre_fasta_lines
        self.features = features
        self.embedded_fasta = embedded_fasta


def read_gff(path: str | Path) -> GffDocument:
    path = Path(path)
    with open_text(path) as handle:
        all_lines = list(handle)
    if not all_lines:
        raise InputFormatError(f"GFF3 file is empty: {path}")

    fasta_at = next(
        (index for index, line in enumerate(all_lines) if line.startswith("##FASTA")),
        None,
    )
    before = all_lines if fasta_at is None else all_lines[:fasta_at]
    features: list[Feature] = []
    for line_number, raw in enumerate(before, 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            raise InputFormatError(
                f"expected 9 tab-separated GFF3 columns at {path}:{line_number}, found {len(fields)}"
            )
        seqid, source, feature_type, start, end, score, strand, phase, attrs = fields
        try:
            interval = Interval.from_gff(int(start), int(end))
        except ValueError as error:
            raise InputFormatError(f"{path}:{line_number}: {error}") from error
        if strand not in {"+", "-", ".", "?"}:
            raise InputFormatError(f"invalid strand {strand!r} at {path}:{line_number}")
        if phase not in {"0", "1", "2", "."}:
            raise InputFormatError(f"invalid phase {phase!r} at {path}:{line_number}")
        features.append(
            Feature(
                seqid=seqid,
                source=source,
                type=feature_type,
                interval=interval,
                score=score,
                strand=strand,
                phase=phase,
                attributes=parse_attributes(attrs),
                raw_line=raw.rstrip("\r\n"),
            )
        )
    embedded = None
    if fasta_at is not None:
        embedded = parse_fasta_lines(all_lines[fasta_at + 1 :], str(path))
    return GffDocument(path, before, features, embedded)


def feature_line(feature: Feature) -> str:
    start, end = feature.interval.to_gff()
    return "\t".join(
        (
            feature.seqid,
            feature.source,
            feature.type,
            str(start),
            str(end),
            feature.score,
            feature.strand,
            feature.phase,
            format_attributes(dict(feature.attributes)),
        )
    )


def write_fasta_records(
    handle: TextIO,
    records: Iterable[tuple[str, str]],
    width: int = 60,
) -> None:
    for identifier, sequence in records:
        handle.write(f">{identifier}\n")
        if width <= 0:
            handle.write(sequence + "\n")
        else:
            for offset in range(0, len(sequence), width):
                handle.write(sequence[offset : offset + width] + "\n")
