# StORF-Reporter 2.0.0 release record

Date: 2026-07-19
Primary objective: correctness and scientific usefulness

This document is the authoritative record of the v2 release work. It
describes implementation changes, scientific policy, compatibility decisions,
validation evidence, integration artefacts, known limitations, and tasks that
require a larger calibration study.

## Executive summary

v2 replaces the independent, positional-list implementations used by the old
Reporter, Finder, and UR Extractor with one typed analysis core. The same code
now powers the unified `storf` command and all historical command names.

The most important policy change is that a stop-to-stop ORF is treated as a
**candidate**, not automatically as a missed gene. Unpromoted candidates use
`ORF` (or user-selected `StORF`). The default high-tier feature replacement
changes a candidate to `CDS` only after successful complete-ORF structural
validation. Its reported locus is unchanged by promotion.

The release fixes confirmed silent false negatives, incorrect UR coordinates,
case-sensitive prediction, inactive ambiguous-base filtering, broken gzip and
packaging paths, invalid CDS phase, and numerous CLI modes that previously
crashed or ignored their options.

## Release scope
- coordinate and interval rewrite;
- UR extraction rewrite;
- StORF, Short-StORF, and Con-StORF discovery on both strands;
- alternative genetic codes 1, 4, and 11;
- transparent evidence scoring and candidate tiers;
- deterministic overlap filtering and identifiers;
- conservative CDS promotion;
- mandatory augment input/output validation and candidate CDS validation;
- unified CLI and legacy compatibility adapters;
- GFF3/FASTA validation;
- gzip, path, manifest, and atomic-output handling;
- Pyrodigal `longest`, `individual`, and `meta` modes;
- candidate extraction and independently runnable BLAST/DIAMOND filtering and
  CDS promotion;
- HTML, TSV, FASTA, GFF3, and JSON output products;
- test, lint, wheel-build, and clean-install release gates;
- nonempty quickstart data;
- citation, changelog, and release documentation.


## Confirmed defects corrected

### Silent terminal-candidate loss

The old annotation filler did not flush candidates after the last feature or at
`##FASTA`. A valid candidate after the final annotation could be silently
discarded with exit status 0.

v2 extracts all URs before candidate discovery and writes candidates from a
complete result object. Output no longer depends on encountering a subsequent
feature line. A regression fixture places the only StORF after the final CDS.

### Incorrect UR coordinate model

The old extractor mixed zero- and one-based coordinates, included annotated
bases in some gaps, used the last feature start for a terminal UR, assumed input
feature order, and did not merge overlaps.

v2 uses one convention internally:

```text
zero-based, half-open [start, end)
```

GFF3 is converted at the I/O boundary only. Annotation intervals are validated,
sorted, and merged before their complement is calculated. UR thresholds apply
to the unextended core; flank extensions are stored independently and clipped
to sequence bounds.

Correct edge behaviour is now tested:

- no selected annotations → the whole contig is one UR;
- a selected whole-contig annotation → no UR;
- unsorted and overlapping annotations → deterministic merged complement;
- the final UR starts after the final annotation end;
- out-of-bounds features → explicit failure.

### Collapsed `ID=gene` selector and Ensembl annotation leakage

An early v2 adapter incorrectly converted the historical attribute query
`-gene_ident ID=gene` into the exact feature type `gene`. This changed the
meaning of Ensembl mode: records such as `ncRNA_gene`, `rRNA_gene`,
`tRNA_gene`, and `pseudogene` whose `ID` attributes began with `gene` were no
longer treated as occupied annotation.

The original selector system is restored in the shared core and exposed as
`--gene-ident`:

- `CDS` is the default exact feature type;
- comma-separated values are exact, case-sensitive GFF3 feature types;
- `Prokka` and `Bakta` select the historical broad feature list; and
- `ID=gene` is an attribute-prefix query across feature types.

The unified CLI, legacy Reporter and UR-Extractor, `StORF_Adder`, Galaxy, and
Nextflow now use the same selector. The selected expression is recorded in the
manifest and in generated GFF3 comments. Attribute queries cannot be mixed
with type names in one expression, preventing ambiguous selection semantics.

For the bundled Ensembl E. coli genome, the corrected `ID=gene` policy marks
4,257 annotation records as occupied: 4,053 `gene`, 89 `tRNA_gene`, 61
`ncRNA_gene`, 32 `pseudogene`, and 22 `rRNA_gene` records. Under the final
start-ranking and coding-overlap policy, the corrected run has 2,733 URs and
3,439 retained candidates, compared with 2,613 URs and 3,882 candidates under
the incorrect type-only mask. No retained candidate overlaps
one of the restored RNA/pseudogene annotations by more than the configured
50-nt annotation-overlap allowance.

### Incorrect Con-StORF construction, validation, and promotion

The first v2 implementation generated every span from stop
`i` to stop `i+2`, even when the two component StORFs did not independently
meet the ordinary length policy. It could represent only one middle stop rather
than the multi-stop maximal runs supported by the original algorithm. Because
start selection searched past the middle stop, a downstream start could also
turn a Con-StORF into a duplicate of an ordinary StORF. On the bundled Ensembl
E. coli run, 18 such duplicate Con-StORFs were incorrectly promoted to CDS,
contradicting the documented policy.

The corrected implementation:

- creates one maximal Con-StORF from each run of at least two qualifying
  adjacent ordinary StORFs;
- retains one or more middle-stop offsets, codons, and oriented genomic
  coordinates;
- searches for starts only before the first middle stop;
- validates Con-StORFs with the explicit `con-storf-readthrough-v1` policy;
- prohibits both heuristic and homology-only CDS promotion for Con-StORFs; and
- keeps the documented takeover behaviour: enabled Con-StORFs compete with and
  may replace ordinary StORFs in overlap filtering.

The same audit found that selecting the first available start and filtering on
full stop-to-stop intervals could suppress compatible known coding calls.
Start sites are now ranked using RBS motif/spacing, codon prior, retained coding
length, and the size-class minimum. Candidate overlap uses ranked coding
intervals by default, with an explicit historical candidate-coordinate option.
At the inspected locus this retains the exact Prokka `flhA_1` coding interval
at `249694-251433 (-)` together with a non-overlapping opposite-strand coding
candidate near the independent `lafU` call.

The regenerated Ensembl E. coli Con-StORF run contains 3,231 retained
candidates: 3,019 ordinary StORFs and 212 Con-StORFs. Of the Con-StORFs, 180
pass readthrough-candidate validation, 32 fail it, and zero are promoted to
CDS. Middle-stop counts are 189 with one, 20 with two, and 3 with three. The
ordinary run retains 3,439 candidates; the lower Con-enabled total confirms the
documented takeover/replacement semantics rather than additive semantics.

### Case-dependent results

All nucleotide input is normalised to uppercase once during parsing. Identical
uppercase and lowercase fixtures now produce identical candidates on both
strands.

### Inactive ambiguous-base filter

The legacy `check_non_standard_chars` function survived, but its call was
removed during the v1.4.1 filtering rewrite. v2 calculates the non-ACGT
fraction for every candidate and applies `--ambiguous-threshold` before overlap
selection. The rejected count is recorded in the run manifest.

The comparison is “greater than the threshold,” so a threshold of zero permits
only candidates containing no ambiguous bases.

### Short-StORF range and short-only crash

Short-StORFs now require both `short_min <= length <= short_max`. The default is
30–120 nt. Short-only selection filters typed candidate records and no longer
attempts an unhashable intersection of dictionary items containing lists.

The release audit found a second adapter regression: historical `Nolap` and
`Olap` had been collapsed into one filter, and the full legacy Reporter ignored
`-short_storfs_only`. `Nolap` again filters ordinary and Short-StORFs together;
`Olap` filters the two classes independently and permits cross-class overlap;
short-only now implies short discovery in both legacy adapters. The unified
equivalents are `--short-only` and `--short-overlap exclude|allow`.

### Candidate-versus-coding sequence extraction

Candidate classification and FASTA representation are independent. Default
nucleotide and amino-acid FASTA both derive from `CandidateCoordinates` with
the upstream delimiting stop omitted. This applies equally to `ORF` and
promoted `CDS` features; neither validation nor promotion silently substitutes
the ranked coding interval.

`--protein-coordinate-mode coding` explicitly requests ranked coding protein
FASTA during `augment`/`find`. `storf extract --sequence-region coding` and the
legacy extractor's matching option provide the same explicit post-run choice.
Discovery retains candidates without `CodingCoordinates` at their candidate
locus rather than silently dropping them from a coding-mode FASTA.
All modes honour per-feature `GeneticCode` and orient reverse-strand sequence
before translation. The inspection run requires the default protein FASTA to
round-trip exactly through both extractors.

### Historical BLAST query IDs did not join to GFF3 IDs

Historical StORF FASTA headers appended `;UR=...` and other attributes directly
to the identifier. BLAST format 6 retained the entire header as `qseqid`, while
the GFF3 `ID` ended at the first semicolon. Both filter implementations
therefore found qualifying hits but retained zero candidates from the bundled
Mycoplasma example.

Modern stable IDs remain unchanged; historical query names are now normalized
to their GFF3 ID token. Unified `filter` and legacy `StORF-Remover` each retain
the same 70 candidates at bit score 50 from the bundled example.

### Auxiliary merge and legacy batch adapter regressions

`storf-add` previously appended accepted features at the bottom, tested overlap
against each annotation independently rather than total merged annotation
coverage, and had silently changed its historical default selector from `CDS`
to `CDS,rRNA,tRNA`. It now restores the `CDS` default, sums overlap against the
union of selected annotations, rejects duplicate IDs, and inserts accepted
features by contig and coordinate while preserving sorted reference lines.

The legacy Reporter now honours line wrapping, Short-StORF modes, strand output
ordering, and `-py_fasta`; multiple inputs with one custom `-oname` are numbered
instead of overwriting one another. `storf-uncovered` accepts both historical
`-gz True` and flag-style `--gzip` syntax and exposes the original gene-ident
selector for the annotations that occupy sequence.

### Alternative genetic codes

Codes 1, 4, and 11 are explicit and validated. Code 4 treats `TGA` as
tryptophan, so it cannot accidentally terminate a Mycoplasma/Spiroplasma
candidate. Translation and stop discovery use the same table.

Unsupported codes fail rather than silently using the standard table.

### Broken Pyrodigal individual training

The historical `indvidual` typo and untrained finder path are removed.

- `meta`: one metagenomic finder;
- `longest`: train on the longest contig when it is at least 20 kb, otherwise
  use metagenomic mode;
- `individual`: train per contig when at least 20 kb, otherwise use metagenomic
  mode for that contig.

Pyrodigal features are converted directly from its documented one-based,
end-inclusive coordinates into v2 intervals.


### Invalid or misleading GFF3

v2 emits the canonical `##gff-version 3` pragma, nine tab-separated columns,
escaped attributes, unique stable IDs, bounded coordinates, and valid phase.

- `ORF` and `StORF` use phase `.`.
- a promoted `CDS` retains the user-selected displayed candidate locus and uses
  phase `0`; its ranked start-to-stop interval remains explicit in
  `CodingCoordinates`;
- a stop-to-stop region without a start is never emitted as CDS; and
- a candidate with an internal stop, ambiguous coding base, invalid terminal
  stop, non-triplet length, or inadequate coding length is not promoted.

`storf validate` checks parsing, sequence-ID agreement, coordinate bounds,
duplicate IDs, StORF IDs, CDS phase, and sequence alphabet. For v2 candidates
it also checks candidate/coding coordinate parsing and containment, explicit
lengths, displayed-coordinate mode, delimiting and middle stops, start-codon
agreement, and saved CDS-validation consistency. For
StORF-Reporter-produced CDS features it additionally checks the actual oriented
sequence for a supported start, terminal stop, triplet length, ambiguity, and
internal stops under the declared genetic code.

### Broken compression

Plain text and gzip input are detected by gzip magic bytes rather than filename
guessing. Every output ending `.gz` is written through `gzip.open`; compression
is tested by magic bytes and decompression.

### Broken installed commands

The `StORF-Remover` entry point now targets
`StORF_Reporter.StORF_Remover:main`. All command versions exit successfully
without requiring analysis inputs. Wheel smoke tests cover:

- `storf`;
- `StORF-Reporter`;
- `StORF-Finder`;
- `UR-Extractor`;
- `StORF-Extractor`;
- `StORF-Remover`;
- `storf-add`;
- `storf-uncovered`; and
- `storf-lengths`.

### Unsafe and inconsistent CLI Booleans

No command uses `type=eval`. Historical `True`/`False` values use a strict
parser; the unified CLI uses normal flags and `BooleanOptionalAction` where
appropriate.

`-gff False`, amino-acid-only output, gzip, ordinary FASTA input, output
directories, and version flags are covered by subprocess tests.

### Maintenance cleanup before release

The final source audit deliberately removed only code with no source, test,
documentation, entry-point, or package-export references:

- six identical legacy Boolean parsers were replaced by one strict shared
  parser;
- candidate-locus nucleotide selection and candidate/coding protein
  translation were centralized so unified and historical writers cannot drift;
- orphaned v1 helpers for iterative reverse complementation, frame translation,
  non-standard-character checking, embedded-FASTA detection, ORF sorting, and
  unused post-hoc evidence rescoring were removed; and
- an unused text-buffer helper, imports, and tuple bindings were removed, while
  type-only imports and import ordering were normalized.

Supported programmatic compatibility entry points (`STORF_Finder`,
`UR_Extractor.extractor`, and `StORF_Remover.remover`) and comments that define
active reproducibility contracts were retained. These include the historical
`STORF_Finder(..., split_index=3)` call signature, `gene_ident` semantics,
source-checkout import handling, pre-v2 GFF validation, and the maximal
Con-StORF definition. No commented-out implementation blocks, stale
TODO/FIXME markers, placeholder branches, broad lint suppressions, or duplicate
package-version literals remain in the Python implementation. The canonical
package version now lives in `constants.py`; Setuptools reads it directly and
the release test checks distribution metadata against it.

## Candidate definition

For each UR and each strand, v2 examines the three oriented frames. A standard
StORF spans two adjacent in-frame stop codons and includes both stop triplets in
its internal candidate sequence. Candidates entirely within an extension flank
are rejected; every retained candidate must intersect the original UR core.

Con-StORFs restore the original maximal-run definition: at least two adjacent
ordinary StORFs must each meet the ordinary length policy, and the maximal run
is represented by one Con-StORF. It can therefore contain one or several middle
in-frame stops. `--include-consecutive` places these candidates into the same
overlap ranking as ordinary StORFs, so a Con-StORF can take over overlapping
ordinary calls; the result is not necessarily a strict superset.

Partial and run-through modes are deliberately disabled in v2.0.0. Those
sequences are not stop-to-stop ORFs and need a separate, evidence-aware model.
The legacy flags now fail with an explanation instead of emitting incorrect
coordinates or lengths.

## Evidence model

`EvidenceScore` is an interpretable prioritisation heuristic. It is **not** a
posterior probability, confidence estimate, false-discovery rate, or substitute
for biological validation.

Every supported in-frame start that leaves the required size-class coding
length is ranked. Selection uses RBS motif strength and spacing, start-codon
prior, and retained coding length, with genomic order only as the final tie
break. RBS context may extend before the delimiting stop. Con-StORF start search
ends at the first middle stop. The current evidence components are:

| Component | Contribution |
|---|---:|
| Coding length | `min(25, max(0, (coding_nt - 30) / 12))` |
| UR-core coverage | `20 × core_fraction` |
| Ambiguous bases | `-50 × ambiguous_fraction` |
| ATG start | 20 |
| GTG start | 15 |
| TTG start | 12 |
| Other table-supported start | 8 |
| Simple Shine–Dalgarno motif | 15 |
| Each external 0–1 evidence value | `20 × value` |

The total is clipped to 0–100. External values are deliberately visible by
name and are never interpreted as independent probabilities.

Tiers:

- `high`: score at least 65;
- `medium`: score at least 40;
- `exploratory`: score below 40.

The evidence TSV contains the component dictionary for every retained
candidate. GFF3 attributes state `ScoreModel=heuristic-v1-not-a-probability`.

The RBS component is a deliberately simple motif search, not a trained
organism-specific initiation model. It should be replaced or calibrated before
being used for strong claims.

## Complete-ORF validation and CDS promotion

The `complete-orf-v1` policy separates structural validity from biological
evidence. A proposed coding interval must:

1. begin with a start codon supported by the selected NCBI table;
2. have a length divisible by three and meet the size-class coding threshold;
3. end with a stop codon supported by that table;
4. contain no internal in-frame stop; and
5. contain no ambiguous nucleotide.

### Feature type and locus are independent

CDS promotion is a feature-classification decision, not a coordinate rewrite.
For every retained ordinary StORF, the configured evidence tier and structural
validation can replace column 3 from `ORF`/`StORF` to `CDS`, but columns 4–5
remain in the user-selected locus representation. Under the default
`stop-to-stop` mode they describe the stop-derived locus with the upstream
delimiting stop omitted, exactly as historical `-spos False` output did. This
is true whether `CDSValidation` passes or fails and whether the feature remains
`ORF` or becomes `CDS`.

`CandidateCoordinates` always preserves the literal two-stop boundary;
`CodingCoordinates` records the ranked start-to-terminal-stop interval used for
structural validation and explicitly requested coding extraction.
`--include-first-stop` adds the
upstream stop to default GFF3 and candidate nucleotide FASTA output.
`--gff-coordinate-mode start-to-stop` is the explicit way to make GFF3 columns
use the ranked coding locus. Neither automatic tier promotion nor homology
promotion overrides these coordinate choices.

The GFF3 records `CandidateCoordinates`, `CodingCoordinates`,
`CandidateLength`, `CodingLength`, candidate-validation fields, and CDS-
validation fields for every candidate. `CandidateCoordinates` and
`CandidateLength` preserve the literal full stop-to-stop candidate. Historical
default output omits the upstream delimiting stop codon from GFF3 columns and
candidate nucleotide FASTA; `--include-first-stop` restores the literal full
interval. `--gff-coordinate-mode start-to-stop` selects the ranked coding
interval. A qualifying promotion changes only the feature type to `CDS`; it
does not override either selected locus representation. The independently
validated coding sequence remains explicit in `CodingCoordinates`.

Con-StORFs use `con-storf-readthrough-v1`: a selected start must precede the
first middle stop, the first translated segment must meet the size threshold,
and all candidate boundaries and sequence characters must be valid. One or
more expected middle stops remain incompatible with conventional complete-CDS
validation. Con-StORFs are explicitly excluded from automatic and homology-
only CDS promotion even when their readthrough-candidate validation passes.

`--promote-cds {high,medium,exploratory}` requires both this structural
validation and the requested evidence tier. The unified workflow defaults to
the high tier; `--no-promote-cds` disables feature-type replacement. Neither
choice changes the selected reported coordinates.
Structural validity is necessary for representing a complete CDS but is not
evidence that the locus is expressed, functional, or conserved.

Explicit coding-mode protein output treats a table-supported alternative
initiation codon (for example GTG under table 11) as initiator methionine.
Default candidate-locus translation does not reinterpret an internal codon as
the sequence start. Internal occurrences retain their ordinary amino-acid
translation.
Expected Con-StORF middle stops are emitted as `*` by default;
`--remove-inframe-stops` removes those protein symbols only.

## End-to-end augment validation

`storf augment` is the complete annotation-aware workflow requested for both
single and batch use. It accepts:

- one GFF3 plus a separate FASTA;
- one combined GFF3 with an embedded `##FASTA` section;
- a directory containing multiple unambiguous GFF3/FASTA or combined-GFF3
  samples, optionally recursively; or
- a sample manifest for explicit cross-directory pairing.

For each sample, augment now performs input validation, UR extraction, StORF
discovery and prioritisation, candidate structural validation, output writing,
and post-write GFF3/FASTA validation. In addition to the augmented GFF3 and
candidate products, it writes UR GFF3/FASTA intermediates and a validation JSON
report. The report gives input/output issues and counts candidates that passed
the CDS structural policy, failed it, and were actually promoted. Input errors
stop analysis; warnings are reported without being silently discarded.

## Coordinate-aware GFF3 insertion and evidence-tier partitioning

The first v2 writer copied the complete original feature section and then
appended every StORF immediately before `##FASTA`. Although legal GFF3, that
made loci difficult to inspect and could surprise downstream software expecting
coordinate-ordered features.

The augmented writer now inserts candidates by sequence ID and output genomic
coordinate. For already sorted input GFF3, original lines, metadata comments,
and feature-group separators are preserved while StORFs are inserted between
their genomic neighbours. If an input feature section is unsorted, the output
is explicitly normalized by contig and coordinate and records that fact in a
GFF3 comment.

`--primary-tier {high,medium,exploratory}` independently controls which
candidates enter the primary augmented GFF3 and matching candidate FASTA/TSV
files. The default `exploratory` retains all candidates for compatibility.
`medium` includes medium and high calls; `high` includes high calls only.
Candidates below the chosen threshold are not discarded: they are written to
candidate-only `_StORFs_below-<tier>.gff`, nucleotide FASTA, amino-acid FASTA,
and evidence TSV review files. The HTML report lists all retained candidates
and their destination, while validation JSON records primary and review counts
and validates both GFF3 partitions.

Tier placement and CDS promotion remain separate decisions. The recommended
conservative combination is `--primary-tier high --promote-cds high`: only
high-tier candidates enter the main annotation, structurally valid high-tier
calls become CDS, and all lower-tier candidates remain available for review.

## Annotation-overlap policy

UR extensions allow discovery of candidates crossing an annotation boundary.
v2 records how many candidate nucleotides lie outside the UR core. Candidates
over the configurable `maximum_annotation_overlap` are filtered before
candidate-candidate overlap selection.

Candidate overlap filtering is deterministic and ranks by:

1. evidence score;
2. candidate length;
3. sequence ID and genomic coordinates;
4. strand; and
5. stable candidate ID.

Modes are `none`, `same-strand`, and `all`. The manifest separately reports
annotation-overlap and candidate-overlap removals. Conflicts use ranked coding
intervals by default so long pre-start regions do not suppress compatible
coding calls. `--overlap-coordinate-basis candidate` restores full stop-to-stop
conflict measurement.

## Stable identifiers and reproducibility

Candidate IDs are derived from sequence ID, zero-based interval, strand,
genetic code, and candidate type using SHA-256. They do not depend on an input
directory or output filename.

Each JSON manifest records:

- tool and schema versions;
- sample name;
- input paths and SHA-256 hashes;
- every analysis parameter;
- sequence, annotation, UR, raw candidate, retained candidate, rejection, and
  tier counts;
- output paths; and
- the evidence-score scientific notice.

Writes are atomic within the destination directory. Existing outputs fail
unless overwrite is explicitly enabled.

## Unified interface

```text
storf augment
storf find
storf extract-urs
storf extract
storf filter
storf validate
storf install-check
storf citation
```

### Input discovery

`augment` accepts one GFF3, a directory, recursive directory discovery, or a
TSV manifest. Matching FASTA stems are preferred. Ambiguous pairings fail and
request `--fasta` or a manifest.

Named `default`, `prokka`, `bakta`, `ensembl`, and `cds-only` profiles make the
annotated-region policy explicit. `--gene-ident` overrides a profile using the
historical selector syntax: exact comma-separated types, `Prokka`/`Bakta`, or
an attribute-prefix query such as `ID=gene`. `--exclude-types` remains an exact
type-list compatibility alias.

Bakta directories containing equivalent `.gff` and `.gff3` outputs are
deduplicated by the legacy adapter, with `.gff3` preferred.

### Output products

Each augmented sample receives:

- augmented GFF3;
- unannotated-region GFF3 and FASTA;
- StORF nucleotide FASTA;
- StORF amino-acid FASTA;
- evidence TSV;
- HTML review report;
- JSON validation report; and
- JSON manifest.

When `--primary-tier` excludes any retained candidates, the sample also
receives a candidate-only below-tier GFF3 plus matching nucleotide FASTA,
amino-acid FASTA, and evidence TSV review files.

Prokka/Bakta output-directory mode also writes a combined original-CDS plus
StORF FASTA. Original CDS extraction respects strand and `transl_table` where
present.

### Independent alignment-evidence filter

`storf filter` reads standard 12-column BLAST or DIAMOND format-6 output after
discovery. Qualifying hits must satisfy configurable bit-score, e-value,
percent-identity, and aligned-length thresholds; the strongest qualifying hit
per candidate is stored in explicit GFF3 homology attributes.

The default mode preserves the historical filtering meaning by removing
StORFs without qualifying hits. `--keep-unmatched` instead retains them at
their existing type. `--promote-cds` requires the matching FASTA unless it is
embedded, revalidates the oriented coding sequence, and promotes only supported
complete ORFs. This creates a separately reproducible route from the initial
mixed candidate annotation to a homology-supported mixed ORF/CDS annotation.
Homology support is deliberately not labelled experimental confirmation.

## Legacy compatibility

The old command names remain, but incorrect semantics were not preserved.

| Historical behaviour | v2 behaviour |
|---|---|
| Default StORFs written as CDS without validation | Only structurally valid high-tier candidates default to `CDS`; others remain `ORF`/`StORF` |
| `-f_type CDS` promoted every StORF | Only structurally valid complete candidates are CDS |
| `type=eval` Boolean parsing | Strict `True`/`False` parser |
| Custom stop strings silently mixed with fixed translation | Standard sets map to code 11 or 4; other sets fail |
| Partial/run-through modes emitted questionable records | Explicitly disabled pending a dedicated model |
| Different commands used different coordinate logic | All commands call the same v2 core |
| `-odir` often required trailing slash | All paths use `pathlib` joins |
| Output names/IDs depended on filesystem paths | IDs are coordinate/content-policy stable |
| Empty input could succeed | Empty/malformed input fails |
| Errors were often caught as gzip failures | Parsing and compression errors remain distinct |

The user-edited legacy `-sout` help text is preserved as:
“Default - False: Print out StORF sequences separately”.

## Auxiliary-tool status

- `storf-add` replaces the removed ORForise dependency with a local,
  standards-aware GFF3 merge filtered by annotation overlap.
- `storf-uncovered` calculates sequence remaining after original annotations
  and reported ORFs from the whole genome, rather than manipulating fragile UR
  identifiers.
- `storf-lengths` now includes the final FASTA record, supports gzip through
  shared I/O, handles one-record input, and no longer requires NumPy.
- `StORF-Extractor` accepts any v2 source/ID rather than only a hard-coded old
  source string.
- `StORF-Remover` accepts floating-point bit scores, validates 12-column
  BLAST/DIAMOND rows, preserves embedded FASTA, and reports kept/removed counts.

The `aux_tools/` directory is publication-analysis material, not installed
product code. Its status and reproducibility boundary are documented in
`aux_tools/README.md`; hard-coded historical scripts are not part of release
validation.

## Packaging and distribution

- Package version is PEP 440 compliant: `2.0.0`.
- Runtime NumPy dependency removed.
- Pyrodigal constrained to the supported major line: `>=3,<4`.
- Isolated builds require only Setuptools and wheel.
- Every console entry point is present in the built wheel.
- Every public and auxiliary entry module also supports direct source-tree
  execution (`python cli.py`, legacy `python SCRIPT.py`, and
  `python Tools/SCRIPT.py`) without an installed package or `PYTHONPATH`.
- Direct execution establishes package context before loading shared modules;
  it does not maintain a second set of fallback imports that could diverge
  from the wheel/conda implementation.
- The repository contains staging artifacts for Bioconda, Galaxy, Nextflow,
  bio.tools, and RO-Crate metadata. Container, Apptainer, Snakemake, and
  WorkflowHub packaging remain release tasks rather than being claimed present.
- The Galaxy augment wrapper defaults to no CDS promotion; users must select a
  tier explicitly, matching the unified CLI's conservative default.
- `CITATION.cff`, changelog, CI, contribution guidance, and release automation
  are present.

The final release must use one immutable Git tag for GitHub, PyPI, containers,
and Zenodo. Do not repeat the v1.4.7 state in which PyPI had no matching latest
GitHub release/tag.

## Validation evidence

### Automated suite

On 2026-07-20:

```text
ruff check src tests
All checks passed

pytest -q
79 passed
```

Coverage includes:

- interval/GFF round trips;
- no annotation and whole-contig annotation;
- unsorted and overlapping annotations;
- terminal UR extraction;
- bounded extensions;
- exact feature-type, Prokka/Bakta alias, and `ID=gene` attribute selectors;
- unified/profile/legacy parity for attribute-based UR extraction;
- tier-controlled CDS feature replacement and explicit opt-out;
- out-of-bounds rejection;
- lowercase equivalence;
- reverse-strand coordinate mapping;
- ambiguous-base rejection;
- exact Short-StORF bounds;
- legacy Short-StORF-only and `Nolap`/`Olap` adapter semantics;
- genetic code 4;
- alternative initiation-codon translation as methionine at the selected start;
- maximal Con-StORF typing with one or several middle stops;
- start-before-first-middle-stop readthrough validation;
- minimum-length and RBS-aware alternative-start ranking, including context
  before the delimiting stop;
- evidence components and CDS promotion;
- complete-ORF promotion rejection for Con-StORFs with internal stops;
- homology-only CDS-promotion rejection for Con-StORFs;
- deterministic coding- and candidate-coordinate overlap filtering;
- stop-derived and start-to-stop GFF coordinate modes, opt-in inclusion of the
  upstream delimiting stop, and explicit candidate and coding lengths;
- second-line GFF3 features;
- duplicate FASTA IDs;
- non-IUPAC input rejection;
- gzip magic detection;
- final-UR end-to-end augmentation;
- real gzip outputs;
- manifest hashes and notices;
- Pyrodigal individual mode;
- all version entry points;
- all nine direct source-script entry points with `PYTHONPATH` removed;
- legacy default output directories;
- standalone Finder gzip/amino-acid-only/ordinary-FASTA paths;
- `-gff False`; and
- unified CLI dry-run and augmentation;
- coordinate-interleaved StORF placement between neighbouring annotations;
- high-tier primary output with lower-tier GFF3/FASTA review partitioning;
- automatic augment preflight/postflight validation and UR/validation outputs;
- invalid input rejection before augment analysis;
- independent homology-supported CDS promotion with FASTA revalidation;
- StORF CDS start-codon validation;
- v2 candidate/coding coordinate and length validation;
- exact unified/legacy extractor protein round trips;
- historical semicolon-suffixed BLAST query-ID joins;
- total-overlap and coordinate-order behavior in `storf-add`;
- collision-free custom output names for multiple legacy inputs;
- directory discovery that excludes prior UR/Finder/Reporter derivatives;
- compressed GFF3/FASTA directory pairing.

### Full subtool inspection run

`examples/run_all_subtools.sh` was run from a fresh output directory on
2026-07-20. It executed every unified command, all five historical adapters,
Pyrodigal mode, and all three auxiliary commands on bundled data. The final
cross-command assertions reported:

```text
unified_legacy_ecoli_candidates  3439
unified_legacy_urs               2733
unified_legacy_find_candidates   3191
historical_blast_hits_retained     70
con_storfs                        212
promoted_cds_locus_preserved      130
candidate_protein_locus_verified 3439
protein_roundtrip_records        3439
```

The unified and legacy full workflows have identical candidate IDs,
`CandidateCoordinates`, and displayed loci. Their feature types can differ
because unified `augment` applies validated high-tier CDS replacement by
default, but that replacement does not alter coordinates. Unified/legacy UR extraction
and FASTA-only discovery are exact. Both extractors reproduce every workflow
protein, both filters retain the same historical evidence-supported set, all
Con-StORFs have one or more middle stops and zero are CDS, and the `storf-add`
result is coordinate ordered and independently valid. The runnable commands
and manual inspection checklist are documented in
`docs/SUBTOOL_INSPECTION_RUNS.md`.

The final locus-policy rerun additionally asserted the reported interval of
every E. coli StORF against its strand-aware stop-derived expectation and
confirmed that 130 promoted CDS features retained those loci. It also required
at least one promoted feature to have a distinct `CodingCoordinates` interval,
so the check cannot pass through accidental equality when the selected start is
immediately after the upstream stop.

The final candidate-FASTA-policy rerun translated every one of the 3,439
default nucleotide records directly from its stop-derived candidate locus and
required exact equality with `_StORFs_aa.fasta`. Both extractors reproduced the
same 3,439 records. This prevents either the `CDS` type or the presence of
`CodingCoordinates` from silently selecting coding-coordinate FASTA.

### Real bundled dataset

The Prokka E. coli fixture completed with the broad Prokka selector and final
start/overlap defaults:

```text
sequences                         1
annotations                   4362
unannotated regions           2800
candidates before filtering   7110
filtered by annotation overlap 200
filtered by candidate overlap 3687
retained candidates           3223
high tier                       28
medium tier                   1383
exploratory tier              1812
```

The augmented output passed `storf validate` with no errors. The low high-tier
count and large exploratory set demonstrate why candidate/CDS separation is
scientifically necessary.

These counts are a v2 baseline, not a claim that every candidate is a gene and
not a direct comparison with v1 because the coordinate and filtering semantics
have intentionally changed.


## Known limitations and deliberate deferrals

### Calibration

Evidence-tier thresholds are transparent but not empirically calibrated. Before
the final 2.0 paper/application release, evaluate precision–recall and false
discovery against held-out known genes, curated small proteins, shuffled
negatives, and independent proteomic/Ribo-seq evidence.

### Circular replicons

Origin-spanning annotations and candidates require a circular interval model.
v2.0.0 fails on out-of-bounds features rather than guessing. Circular
support should be introduced with dedicated origin-crossing fixtures and GFF3
representation tests.

### Additional genetic codes

Only tables 1, 4, and 11 are implemented. Add further tables from authoritative
NCBI definitions together with translation and stop-discovery tests.

### RBS and coding potential

The RBS motif is simple and the length term is heuristic. Pyrodigal candidate
coding scores, codon usage, GC deviation, organism-specific initiation models,
and conservation are not yet computed automatically. External values can be
imported through the evidence TSV without being hidden inside the score.


