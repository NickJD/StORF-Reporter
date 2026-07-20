# Migrating from v1 to v2

## Recommended command

Replace a mode-coupled invocation such as:

```bash
StORF-Reporter -anno Prokka Out_Dir -p prokka-output/
```

with:

```bash
storf augment --input prokka-output/ --profile prokka --outdir storf-results/
```

The historical command still works and calls the same v2 core.

## Deliberate output changes

- Unpromoted candidates are `ORF`, not automatically CDS. The default high-tier
  policy changes structurally valid qualifying features to `CDS`;
  `--promote-cds` selects the evidence threshold and `--no-promote-cds`
  disables feature-type replacement.
- `-gene_ident` is available as `--gene-ident`: exact type lists and `Prokka`
  retain their historical meaning, while `ID=gene` correctly queries the
  attribute column across feature types.
- Candidate IDs no longer contain filesystem paths.
- GFF3 phase is standards-compliant.
- UR boundaries can differ because annotations are sorted and merged before
  their complement is calculated.
- Lowercase and uppercase FASTA now produce the same output.
- Ambiguous-base filtering is active.
- Candidate overlap prefers evidence score before length.
- Candidate overlap uses ranked coding intervals by default; use
  `--overlap-coordinate-basis candidate` to reproduce full stop-to-stop
  conflict measurement.
- `CandidateCoordinates` preserves the literal full stop-to-stop interval.
  Historical default output omits the upstream delimiting stop from GFF3
  columns and candidate nucleotide FASTA; `--include-first-stop` restores the
  literal full interval. Use `--gff-coordinate-mode start-to-stop` for ranked
  coding coordinates. CDS promotion changes feature type only and does not
  override the selected reported-locus mode.
- Candidate nucleotide and protein FASTA both use the stop-derived locus by
  default. Coding-coordinate protein FASTA now requires explicit
  `--protein-coordinate-mode coding` or `storf extract --sequence-region coding`.
- The ambiguous `Length` attribute is replaced by `CandidateLength` and
  `CodingLength`.
- `--include-consecutive` creates maximal multi-stop Con-StORFs and lets them
  compete with ordinary StORFs. Results can therefore be replacements rather
  than a strict superset of an ordinary run.
- All `.gz` outputs are actually compressed.
- Empty or mismatched inputs fail rather than producing empty success output.

If a downstream tool requires CDS features, use `--promote-cds high` initially
and document the promotion threshold. The release record explains the score and
remaining calibration limitations.

## Comparing versions

Do not compare only the total number of GFF3 rows. Compare:

- UR core intervals;
- raw and filtered candidate counts;
- reasons for removal;
- candidate sequences and stable loci;
- evidence tiers; and
- downstream results using the original, high-tier, and full candidate sets.

Store the v1 command, v2 manifest, and both input hashes with the comparison.
