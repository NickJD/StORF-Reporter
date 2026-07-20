# Scientific interpretation guide

## What StORF-Reporter establishes

StORF-Reporter establishes that an unannotated region contains a reading frame
bounded by in-frame stop codons under the selected genetic code. It also reports
simple, inspectable evidence that may make the candidate more useful for review.

It does not establish that the sequence is translated, functional, conserved,
or correctly initiated.

## Structural CDS validation

`storf augment` validates inputs before UR extraction and validates its final
GFF3 before reporting success. Every candidate is also assessed for possible
CDS representation. A candidate passes `complete-orf-v1` only when its proposed
coding sequence has a supported in-frame start, a valid terminal stop, triplet
length, no internal in-frame stop, no ambiguous nucleotide, and sufficient
coding length for its selected size class.

Passing these checks means that the feature can be represented as a complete
CDS without violating those structural rules. It does not demonstrate that the
feature is a gene. The default high-tier promotion policy requires both
structural validation and the selected heuristic evidence tier;
`--promote-cds` changes that tier and `--no-promote-cds` disables feature-type
replacement. Calls that do not qualify remain `ORF`/`StORF` features. Con-StORFs
are evaluated separately as potential readthrough candidates. They pass
`con-storf-readthrough-v1` only when the ranked start is before the first of one
or more middle stops and the first translated segment is long enough. They
remain ineligible for automatic CDS promotion because their coding path still
contains internal stops.

## Con-StORF interpretation

Con-StORFs are maximal runs of at least two adjacent qualifying ordinary
StORFs. `--include-consecutive` puts those candidates into the same evidence and
overlap ranking as ordinary StORFs; it can therefore replace ordinary calls
rather than merely append new calls. Compare manifests and candidate IDs when
comparing runs with and without this option.

A passed Con-StORF candidate validation establishes only that a plausible start
precedes the first middle stop. It is not evidence that any middle stop is
biologically read through. Protein FASTA records retain `*` at middle stops by
default so downstream tools and reviewers can see that structure.
`--remove-inframe-stops` removes only those protein symbols and does not change
nucleotide loci.

## Start and coordinate interpretation

The ranked start uses supported start-codon priors, local RBS motif/spacing, and
retained coding length. It remains a sequence heuristic rather than an
organism-trained start-site model. Candidate-level homology can support a locus
but does not by itself resolve alternative starts.

`CandidateCoordinates` is the literal full stop-to-stop interval. By historical
default, GFF3 `ORF`/`StORF` columns and candidate nucleotide FASTA omit the
upstream delimiting stop codon; `--include-first-stop` restores the literal
full interval. The optional `--gff-coordinate-mode start-to-stop` changes
feature columns to the ranked coding interval. Promotion changes column 3 to
`CDS` without overriding the selected GFF coordinate mode. The validated
reading frame remains explicit in `CodingCoordinates`; `CandidateCoordinates`,
`CandidateLength`, and `CodingLength` preserve both representations in every
mode.

Default nucleotide and protein FASTA both represent the stop-derived candidate
locus with its upstream stop omitted. The `CDS` feature label does not select a
different FASTA interval. Ranked coding protein output is opt-in through
`--protein-coordinate-mode coding` or `storf extract --sequence-region coding`.

## Primary and review output tiers

`--primary-tier` controls output placement independently of CDS promotion:

- `exploratory` includes all retained candidates in the primary annotation;
- `medium` includes medium and high candidates; and
- `high` includes high candidates only.

Candidates below the threshold are written to separate candidate-only GFF3,
nucleotide FASTA, amino-acid FASTA, and evidence TSV review files. They are not
discarded. This allows conservative downstream annotation without losing the
exploratory set needed for sensitivity analyses. Within the primary augmented
GFF3, StORFs are placed beside other features by contig and coordinate rather
than appended at the end.

## Recommended review order

1. Confirm the genetic code and assembly topology.
2. Inspect overlap with existing CDS, RNA, repeat, CRISPR, and low-complexity
   annotations.
3. Review the selected start codon and RBS component.
4. Search candidate proteins against curated and taxonomically appropriate
   databases.
5. Examine conservation across related genomes and pangenome context.
6. Add independent transcriptomic, proteomic, or Ribo-seq evidence when
   available.
7. Promote to CDS only under a documented project-specific policy.

## Evidence TSV imports

Run discovery once to obtain stable IDs. Create a tab-separated file:

```tsv
candidate_id	homology	proteomics	conservation
StORF_0123456789abcd	0.9	1.0	0.8
```

Values must be between zero and one. Column names are preserved in output. On a
second run, use `--evidence-tsv evidence.tsv`. These values are score inputs;
they are not assumed to be independent probabilities.

## BLAST/DIAMOND evidence after discovery

`storf filter` is a separate post-discovery path for standard 12-column
BLAST/DIAMOND format-6 results. Qualifying hits are selected using configurable
bit-score, e-value, percent-identity, and aligned-length thresholds. The best
qualifying hit for each stable candidate ID is recorded in GFF3 attributes.

Use `--keep-unmatched --promote-cds --fasta genome.fna` to retain the complete
candidate set while promoting homology-supported candidates that also pass an
independent sequence-level structural recheck. Omit `--keep-unmatched` to make
a homology-supported-only annotation. A database hit is evidence of similarity,
not proof of translation, and thresholds should be documented with the database
version and taxonomic scope.

## Pangenome analyses

Do not insert every exploratory candidate as CDS before running Roary, Panaroo,
or another pangenome method. That can inflate the accessory genome and fragment
families. Prefer one of these strategies:

- retain candidates as `ORF`/`StORF` and analyse them separately;
- promote only high-tier candidates with independent evidence; or
- run a sensitivity analysis comparing the original annotation, high-tier
  additions, and all candidates.

Always retain the manifest and evidence TSV with downstream results.
