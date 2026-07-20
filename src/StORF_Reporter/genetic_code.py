"""Small, dependency-free genetic-code helpers.

Codes 1, 4 and 11 cover the standard, Mycoplasma/Spiroplasma and bacterial
codes most relevant to StORF-Reporter. Unsupported codes fail explicitly.
"""

from __future__ import annotations

_STANDARD = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}

_TABLES = {1: _STANDARD, 11: _STANDARD}
_TABLE_4 = dict(_STANDARD)
_TABLE_4["TGA"] = "W"
_TABLES[4] = _TABLE_4

_STOPS = {
    1: frozenset({"TAA", "TAG", "TGA"}),
    4: frozenset({"TAA", "TAG"}),
    11: frozenset({"TAA", "TAG", "TGA"}),
}

_STARTS = {1: ("ATG", "TTG", "CTG"), 4: ("ATG", "GTG", "TTG"), 11: ("ATG", "GTG", "TTG")}

_COMPLEMENT = str.maketrans(
    "ACGTRYMKBDHVNacgtrymkbdhvn",
    "TGCAYRKMVHDBNtgcayrkmvhdbn",
)


def supported_codes() -> tuple[int, ...]:
    return tuple(sorted(_TABLES))


def stop_codons(code: int) -> frozenset[str]:
    try:
        return _STOPS[code]
    except KeyError as error:
        raise ValueError(f"unsupported genetic code {code}; supported codes: {supported_codes()}") from error


def start_codons(code: int) -> tuple[str, ...]:
    try:
        return _STARTS[code]
    except KeyError as error:
        raise ValueError(f"unsupported genetic code {code}; supported codes: {supported_codes()}") from error


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1].upper()


def translate(
    sequence: str,
    code: int = 11,
    keep_terminal_stop: bool = True,
    start_as_methionine: bool = False,
) -> str:
    try:
        table = _TABLES[code]
    except KeyError as error:
        raise ValueError(f"unsupported genetic code {code}; supported codes: {supported_codes()}") from error
    sequence = sequence.upper()
    protein = "".join(table.get(sequence[index : index + 3], "X") for index in range(0, len(sequence) - 2, 3))
    if start_as_methionine and protein and sequence[:3] in start_codons(code):
        protein = "M" + protein[1:]
    if not keep_terminal_stop and protein.endswith("*"):
        protein = protein[:-1]
    return protein
