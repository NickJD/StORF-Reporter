# Changelog

All notable changes are documented here. The project follows semantic
versioning from v2 onward.

## 2.0.0 - 20/07/2026

Correctness-first rewrite with one typed analysis core, validated GFF3,
deterministic identifiers, scientific evidence tiers, a unified `storf` CLI,
working legacy adapters, regression tests, release automation, and workflow
integration scaffolding.

Source-checkout execution is supported for every public and auxiliary script
without requiring an editable install or a manually configured `PYTHONPATH`.

`storf augment` now performs mandatory input and output validation, writes UR
intermediates and a validation report, and records complete-ORF structural
validation for every candidate. The default high-tier CDS feature replacement
requires a supported start, triplet structure, terminal stop, no internal stop,
no ambiguous coding base, adequate coding length, and the requested evidence
tier. Promotion does not replace the selected reported locus with coding
coordinates.

Candidate nucleotide and amino-acid FASTA now both default to the stop-derived
locus with its upstream stop omitted. A `CDS` label never changes FASTA
coordinates implicitly; ranked coding protein output requires explicit
`--protein-coordinate-mode coding` or `--sequence-region coding`.

The original `gene_ident` system is restored across all UR workflows. Exact
feature-type lists, the `Prokka`/`Bakta` alias, and attribute-prefix queries
such as `ID=gene` now share one implementation. This fixes Ensembl profiles
that previously ignored known RNA-gene and pseudogene annotations.

`storf filter` now supports independent BLAST/DIAMOND format-6 thresholds,
best-hit provenance, optional retention of unmatched ORFs, and sequence-checked
promotion of homology-supported candidates to CDS.

Con-StORFs now restore maximal runs of two or more qualifying adjacent StORFs,
including runs with several middle stops. Their dedicated validation requires
the ranked start to precede the first middle stop, and they cannot be promoted
to conventional CDS features automatically or by homology alone. Enabling
Con-StORFs intentionally lets them compete with and replace ordinary StORFs.

Alternative starts are ranked using RBS motif/spacing, start-codon prior,
retained coding length, and the size-class minimum. Candidate overlap uses
ranked coding intervals by default, with a historical stop-to-stop option.
GFF3 output retains every candidate as ORF/StORF by default. Literal
stop-to-stop bounds remain in `CandidateCoordinates`; following historical
behavior, displayed GFF3 and nucleotide FASTA omit the upstream stop unless
`--include-first-stop` is selected. `--gff-coordinate-mode start-to-stop`
remains available, and promotion changes feature type without changing the
selected coordinate representation.
`Length` is replaced by `CandidateLength` and `CodingLength`.

Augmented GFF3 files now interleave StORF features with neighbouring original
features by contig and coordinate. `--primary-tier` can restrict the main
GFF3/FASTA/TSV outputs to high or medium-and-above candidates while preserving
lower-tier candidates in separate review GFF3/FASTA/TSV files.

See [RELEASE_2.0.0.md](RELEASE_2.0.0.md) for the complete record,
compatibility notes, validation evidence, and known limitations.

## 1.4.7 - 15/01/2026

Historical release. Fixed duplicate StORFs caused by incomplete flushing
between contigs. Versions before 2.0 used separate coordinate and output paths
for each command and did not have an executable regression suite.
