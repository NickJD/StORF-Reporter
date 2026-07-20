# StORF-Reporter 2.0 subtool inspection runs

This is a release-oriented set of reproducible runs for every supported
command. It uses only files already present under `Test_Datasets/`, writes to a
new output directory, and finishes with cross-command parity checks. The runs
prioritise biological and coordinate correctness over minimum runtime.

## Run everything

From the repository root:

```bash
bash examples/run_all_subtools.sh /tmp/storf-reporter-2.0-inspection
```

Set `STORF_PYTHON` if the intended release environment is not the `python`
found on `PATH`:

```bash
STORF_PYTHON=/path/to/environment/bin/python \
  bash examples/run_all_subtools.sh /tmp/storf-reporter-2.0-inspection
```

The target directory must be absent or empty. The script does not delete or
replace an existing inspection run. It uses `python -m` so the source tree is
tested even before installation; these modules are the same implementations
used by the installed console commands.

Expected runtime on a laptop is approximately one to two minutes. The two
E. coli augment runs and Pyrodigal run dominate runtime. A successful run ends
with `inspection-summary.tsv` and the message `Complete`.

## What each run exercises

| Command | Example input | Output to inspect | Correctness question |
|---|---|---|---|
| `storf install-check` | active environment | terminal | Is the intended Python/package/dependency active, and which wrappers are installed? |
| `storf citation` | none | terminal | Is the publication citation correct? |
| `storf augment` | embedded Ensembl E. coli GFF3 | `unified/augment/E-coli_Combined/` | Are URs, candidates, CDS promotion, interleaved GFF3, validation, report, and manifest coherent? |
| `storf augment --include-consecutive` | embedded Ensembl E. coli GFF3 | `unified/augment-con-storf/` | Do Con-StORFs replace competing calls, retain middle stops, and remain non-CDS? |
| `storf augment --genetic-code 4` | Mycoplasma combined GFF3 | `unified/myco-code4/` | Is TGA treated according to translation table 4? |
| `storf extract-urs` | multisequence Pseudomonas GFF3/FASTA | `unified/pseudomonas-urs/` | Are all contigs represented with bounded, non-overlapping UR cores? |
| `storf find` | Bakta UR FASTA | `unified/find/` | Does FASTA-only discovery work without annotation context? Coordinates are local to each FASTA record. |
| `storf extract` | new augmented E. coli GFF3 | `unified/extracted-candidates.*` | Do nucleotide and protein sequences both use the stop-derived locus minus its upstream stop by default? |
| `storf filter` | historical Mycoplasma GFF3 plus BLAST format 6 | `unified/myco-hits-only.gff` | Are the 70 qualifying historical query IDs matched despite semicolon-suffixed FASTA headers? |
| `StORF-Reporter` | same E. coli GFF3 | `legacy/reporter/` | Does the compatibility CLI retain the same candidate identity and candidate coordinates? |
| `UR-Extractor` | same E. coli GFF3 | `legacy/ur-extractor/` | Is `-gene_ident ID=gene` identical to the unified Ensembl profile? |
| `StORF-Finder` | same Bakta UR FASTA | `legacy/finder/` | Does the adapter retain the same calls as unified `find` when both use local FASTA coordinates? |
| `StORF-Extractor` | new augmented E. coli GFF3 | `legacy/extractor/` | Does the legacy extractor reproduce all ranked protein sequences exactly? |
| `StORF-Remover` | historical Mycoplasma GFF3 plus BLAST format 6 | `legacy/remover/` | Does the compatibility filter retain the same 70 candidates as unified `filter`? |
| legacy Pyrodigal mode | raw E. coli FASTA | `legacy/pyrodigal/` | Are Pyrodigal CDS, UR-derived StORFs, combined FASTAs, and original Pyrodigal GFF all written? |
| `storf-add` | reference E. coli plus candidate-only GFF3 | `aux/e-coli-added.gff` | Are accepted candidates interleaved by coordinate and filtered by total annotation overlap? |
| `storf-uncovered` | genome FASTA plus augmented GFF3 | `aux/e-coli_UnStORFed_UR.*` | What sequence remains outside selected original annotations and reported StORFs? |
| `storf-lengths` | E. coli UR FASTA | `aux/e-coli-ur-lengths.tsv` | Are count and length statistics plausible, including the final FASTA record? |

`find` has no whole-genome annotation or coordinate frame. Its GFF3 coordinates
therefore refer to each input FASTA record. Use `augment` when genome-coordinate
insertion beside existing annotations is required. For the legacy Finder run,
`-ua False` deliberately requests the same local-coordinate interpretation.

## Installed command equivalents

After installing the package, the principal module invocations correspond to:

```text
python -m StORF_Reporter.cli                   -> storf
python -m StORF_Reporter.StORF_Reporter        -> StORF-Reporter
python -m StORF_Reporter.StORF_Finder          -> StORF-Finder
python -m StORF_Reporter.UR_Extractor          -> UR-Extractor
python -m StORF_Reporter.StORF_Extractor       -> StORF-Extractor
python -m StORF_Reporter.StORF_Remover         -> StORF-Remover
python -m StORF_Reporter.Tools.StORF_Adder     -> storf-add
python -m StORF_Reporter.Tools.Un_StORFed      -> storf-uncovered
python -m StORF_Reporter.Tools.UR_Lenghts      -> storf-lengths
```

Lowercase aliases are also installed for the five historical mixed-case
commands. `storf`, `StORF`, and `STORF` are aliases of the unified CLI; they are
not separate biological subtools.

## Manual inspection checklist

1. Open the augmented E. coli GFF3 and confirm StORF records occur beside
   neighbouring features rather than in one block at the bottom.
2. Compare `CandidateCoordinates`, `CodingCoordinates`, `CandidateLength`, and
   `CodingLength` against GFF columns and the extracted FASTAs. Confirm that CDS
   promotion changes column 3 only: it must not replace default columns 4–5
   with `CodingCoordinates`. The default displayed/NT locus is three bases
   shorter than `CandidateCoordinates` because the upstream stop is omitted.
3. Inspect the validation JSON. It should have status `passed` and no error
   entries; warnings about a historical input header are acceptable.
4. Open the HTML report and inspect candidates across high, medium, and
   exploratory tiers. `EvidenceScore` is a ranking heuristic, not a
   probability.
5. In the Con-StORF GFF3, confirm every Con-StORF has at least one middle stop,
   uses `con-storf-readthrough-v1`, and is not a CDS.
6. Confirm at least one promoted CDS has GFF columns that differ from its
   `CodingCoordinates`; this is intentional because promotion changes type but
   does not override the selected stop-derived locus.
7. Confirm default protein extraction exactly reproduces `_StORFs_aa.fasta`
   from the candidate locus. Run `--sequence-region coding` separately only
   when inspecting the ranked coding translation.
8. Inspect `myco-hits-only.gff` and the legacy remover output. Both should
   retain 70 StORFs from the bundled BLAST table.
9. Inspect the manifest JSON for absolute inputs, SHA-256 hashes, full
   configuration, output paths, and candidate counts.

## Release tests after visual inspection

```bash
ruff check src tests
PYTHONPATH=src python -m pytest -q
python -m build
python -m twine check dist/*
```

For the final distribution check, install the newly built wheel into a clean
environment and rerun `storf install-check` plus `--version` for every wrapper.
