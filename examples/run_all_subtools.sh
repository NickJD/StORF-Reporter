#!/usr/bin/env bash

set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
storf_python="${STORF_PYTHON:-python}"
output_root="${1:-/tmp/storf-reporter-2.0-inspection}"

if [[ -e "$output_root" ]] && [[ -n "$(find "$output_root" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    printf 'Output directory is not empty: %s\nChoose a new directory or move the existing results.\n' "$output_root" >&2
    exit 2
fi

mkdir -p "$output_root"
export PYTHONPATH="$repo_root/src${PYTHONPATH:+:$PYTHONPATH}"
export STORF_INSPECTION_ROOT="$output_root"

run_module() {
    "$storf_python" -m "$@"
}

section() {
    printf '\n[%s] %s\n' "$(date '+%H:%M:%S')" "$1"
}

section "Runtime and citation"
run_module StORF_Reporter.cli install-check
run_module StORF_Reporter.cli citation

section "Unified augment: Ensembl E. coli, default promotion policy"
run_module StORF_Reporter.cli augment \
    --input "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    --profile ensembl \
    --outdir "$output_root/unified/augment" \
    --overwrite

unified_ecoli="$output_root/unified/augment/E-coli_Combined/E-coli_Combined"
run_module StORF_Reporter.cli validate \
    --gff "${unified_ecoli}_StORF-Reporter.gff" \
    --json > "$output_root/unified/augment-validation.json"

section "Unified augment: Con-StORFs compete with ordinary StORFs"
run_module StORF_Reporter.cli augment \
    --input "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    --profile ensembl \
    --include-consecutive \
    --outdir "$output_root/unified/augment-con-storf" \
    --overwrite

section "Unified augment: genetic code 4 Mycoplasma"
run_module StORF_Reporter.cli augment \
    --input "$repo_root/Test_Datasets/Combined_GFFs/Myco_Combined.gff" \
    --profile ensembl \
    --genetic-code 4 \
    --no-promote-cds \
    --outdir "$output_root/unified/myco-code4" \
    --overwrite

section "Unified extract-urs: multisequence Pseudomonas"
run_module StORF_Reporter.cli extract-urs \
    --gff "$repo_root/Test_Datasets/Multi_Contig_Genomes/Pseudomonas_aeruginosa.gff3" \
    --fasta "$repo_root/Test_Datasets/Multi_Contig_Genomes/Pseudomonas_aeruginosa.fa" \
    --gene-ident Prokka \
    --outdir "$output_root/unified/pseudomonas-urs" \
    --overwrite

section "Unified find: pre-extracted Bakta UR FASTA"
run_module StORF_Reporter.cli find \
    --fasta "$repo_root/Test_Datasets/Bakta_E-coli/E-coli_Combined_UR.fasta" \
    --sample bakta-urs \
    --no-promote-cds \
    --outdir "$output_root/unified/find" \
    --overwrite

section "Unified extract: default candidate-locus nucleotide and protein sequences"
run_module StORF_Reporter.cli extract \
    --gff "${unified_ecoli}_StORF-Reporter.gff" \
    --out "$output_root/unified/extracted-candidates.fna" \
    --overwrite
run_module StORF_Reporter.cli extract \
    --gff "${unified_ecoli}_StORF-Reporter.gff" \
    --amino-acid \
    --out "$output_root/unified/extracted-candidates.faa" \
    --overwrite

section "Unified filter: bundled historical BLAST format-6 evidence"
run_module StORF_Reporter.cli filter \
    --gff "$repo_root/Test_Datasets/StORF_Extractor_And_Remover/Myco_UR_StORF-R.gff" \
    --alignments "$repo_root/Test_Datasets/StORF_Extractor_And_Remover/Myco_URs_StORFs_aa_Swiss.tab" \
    --min-bitscore 50 \
    --max-evalue 1e-5 \
    --out "$output_root/unified/myco-hits-only.gff" \
    --overwrite

section "Legacy StORF-Reporter adapter"
run_module StORF_Reporter.StORF_Reporter \
    -anno Ensembl Single_Combined_GFF \
    -p "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    -odir "$output_root/legacy/reporter" \
    -verbose \
    --overwrite

section "Legacy UR-Extractor adapter"
run_module StORF_Reporter.UR_Extractor \
    -gff "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    -gene_ident ID=gene \
    -odir "$output_root/legacy/ur-extractor" \
    -verbose True \
    --overwrite

section "Legacy StORF-Finder adapter"
run_module StORF_Reporter.StORF_Finder \
    -f "$repo_root/Test_Datasets/Bakta_E-coli/E-coli_Combined_UR.fasta" \
    -ua False \
    -aa True \
    -odir "$output_root/legacy/finder" \
    -oname bakta-urs \
    -verbose True \
    --overwrite

section "Legacy StORF-Extractor adapter"
run_module StORF_Reporter.StORF_Extractor \
    -storf_input Combined \
    -p "${unified_ecoli}_StORF-Reporter.gff" \
    -gff_out True \
    -aa True \
    -stop_ident False \
    -odir "$output_root/legacy/extractor" \
    -oname e-coli-extracted \
    -verbose True \
    --overwrite

section "Legacy StORF-Remover adapter"
run_module StORF_Reporter.StORF_Remover \
    -gff "$repo_root/Test_Datasets/StORF_Extractor_And_Remover/Myco_UR_StORF-R.gff" \
    -blast "$repo_root/Test_Datasets/StORF_Extractor_And_Remover/Myco_URs_StORFs_aa_Swiss.tab" \
    -min_score 50 \
    -odir "$output_root/legacy/remover" \
    -oname myco-hits-only \
    -verbose True \
    --overwrite

section "Legacy Pyrodigal plus StORF-Reporter workflow"
run_module StORF_Reporter.StORF_Reporter \
    -anno Pyrodigal Single_FASTA \
    -p "$repo_root/Test_Datasets/Pyrodigal/E-coli.fa" \
    -py_train longest \
    -py_fasta True \
    -py_unstorfed True \
    -aa \
    -odir "$output_root/legacy/pyrodigal" \
    -verbose \
    --overwrite

section "Auxiliary storf-add"
run_module StORF_Reporter.Tools.StORF_Adder \
    --genome_DNA "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    --reference_annotation "$repo_root/Test_Datasets/Combined_GFFs/E-coli_Combined.gff" \
    --gene_ident ID=gene \
    --additional_annotation "$output_root/legacy/extractor/e-coli-extracted.gff" \
    --overlap 50 \
    --output_file "$output_root/aux/e-coli-added.gff" \
    --overwrite
run_module StORF_Reporter.cli validate --gff "$output_root/aux/e-coli-added.gff"

section "Auxiliary storf-uncovered"
run_module StORF_Reporter.Tools.Un_StORFed \
    --fasta_seq "$repo_root/Test_Datasets/Matching_GFF_FASTA/E-coli.fa" \
    --gff "${unified_ecoli}_StORF-Reporter.gff" \
    --gene-ident ID=gene \
    --min-length 50 \
    --output-prefix "$output_root/aux/e-coli" \
    --overwrite

section "Auxiliary storf-lengths"
run_module StORF_Reporter.Tools.UR_Lenghts \
    --URs "${unified_ecoli}_UR.fasta" > "$output_root/aux/e-coli-ur-lengths.tsv"

section "Cross-command parity and scientific invariants"
"$storf_python" - <<'PY'
import os
from pathlib import Path

from StORF_Reporter.genetic_code import translate
from StORF_Reporter.io import read_gff
from StORF_Reporter.models import Interval

root = Path(os.environ["STORF_INSPECTION_ROOT"])


def storfs(path):
    return [feature for feature in read_gff(path).features if feature.source == "StORF-Reporter"]


def ids_and_candidates(path):
    return {
        (
            feature.attributes.get("ID"),
            feature.attributes.get("CandidateCoordinates"),
            feature.strand,
        )
        for feature in storfs(path)
    }


def ids_and_display(path):
    return {
        (feature.attributes.get("ID"), feature.interval, feature.strand)
        for feature in storfs(path)
    }


def fasta(path):
    records = {}
    identifier = None
    for raw in Path(path).read_text().splitlines():
        if raw.startswith(">"):
            identifier = raw[1:].split()[0]
            records[identifier] = ""
        elif raw and identifier is not None:
            records[identifier] += raw
    return records


unified = root / "unified/augment/E-coli_Combined/E-coli_Combined_StORF-Reporter.gff"
legacy = root / "legacy/reporter/E-coli_Combined_StORF-Reporter_Extended.gff"
assert ids_and_candidates(unified) == ids_and_candidates(legacy)
assert ids_and_display(unified) == ids_and_display(legacy)

promoted = [feature for feature in storfs(unified) if feature.type == "CDS"]
assert promoted
for feature in storfs(unified):
    start, end = map(int, feature.attributes["CandidateCoordinates"].split("-"))
    candidate_interval = Interval.from_gff(start, end)
    expected = (
        Interval(candidate_interval.start + 3, candidate_interval.end)
        if feature.strand == "+"
        else Interval(candidate_interval.start, candidate_interval.end - 3)
    )
    assert feature.interval == expected
assert any(
    feature.interval
    != Interval.from_gff(*map(int, feature.attributes["CodingCoordinates"].split("-")))
    for feature in promoted
)

unified_urs = read_gff(
    root / "unified/augment/E-coli_Combined/E-coli_Combined_UR.gff"
).features
legacy_urs = read_gff(root / "legacy/ur-extractor/E-coli_Combined_UR.gff").features
assert [feature.interval for feature in unified_urs] == [feature.interval for feature in legacy_urs]

unified_find = root / "unified/find/bakta-urs_StORF-Reporter.gff"
legacy_find = root / "legacy/finder/bakta-urs.gff"
assert ids_and_candidates(unified_find) == ids_and_candidates(legacy_find)

expected_aa = fasta(
    root / "unified/augment/E-coli_Combined/E-coli_Combined_StORFs_aa.fasta"
)
expected_nt = fasta(
    root / "unified/augment/E-coli_Combined/E-coli_Combined_StORFs.fasta"
)
features_by_id = {feature.attributes["ID"]: feature for feature in storfs(unified)}
assert set(expected_nt) == set(expected_aa) == set(features_by_id)
assert all(
    expected_aa[identifier]
    == translate(
        sequence,
        int(features_by_id[identifier].attributes["GeneticCode"]),
        keep_terminal_stop=False,
        start_as_methionine=False,
    )
    for identifier, sequence in expected_nt.items()
)
assert fasta(root / "unified/extracted-candidates.faa") == expected_aa
assert fasta(root / "legacy/extractor/e-coli-extracted_aa.fasta") == expected_aa

unified_hits = ids_and_candidates(root / "unified/myco-hits-only.gff")
legacy_hits = ids_and_candidates(root / "legacy/remover/myco-hits-only.gff")
assert unified_hits == legacy_hits
assert len(unified_hits) == 70

con_path = (
    root
    / "unified/augment-con-storf/E-coli_Combined/E-coli_Combined_StORF-Reporter.gff"
)
con_features = [
    feature
    for feature in storfs(con_path)
    if feature.attributes.get("StORF_Type") == "Con-StORF"
]
assert con_features
assert not any(feature.type == "CDS" for feature in con_features)
assert all(int(feature.attributes["MiddleStopCount"]) >= 1 for feature in con_features)

added = read_gff(root / "aux/e-coli-added.gff").features
keys = [(feature.seqid, feature.interval.start, feature.interval.end) for feature in added]
assert keys == sorted(keys)

summary = root / "inspection-summary.tsv"
summary.write_text(
    "check\tvalue\n"
    f"unified_legacy_ecoli_candidates\t{len(ids_and_candidates(unified))}\n"
    f"unified_legacy_urs\t{len(unified_urs)}\n"
    f"unified_legacy_find_candidates\t{len(ids_and_candidates(unified_find))}\n"
    f"historical_blast_hits_retained\t{len(unified_hits)}\n"
    f"con_storfs\t{len(con_features)}\n"
    f"promoted_cds_locus_preserved\t{len(promoted)}\n"
    f"candidate_protein_locus_verified\t{len(expected_aa)}\n"
    f"protein_roundtrip_records\t{len(expected_aa)}\n"
)
print(summary.read_text(), end="")
PY

section "Complete"
printf 'Inspect outputs under: %s\n' "$output_root"
printf 'Compact audit summary: %s\n' "$output_root/inspection-summary.tsv"
