# Quickstart fixture

This deliberately small genome contains one known CDS at positions 1–100 and
one 99-nt, positive-strand stop-to-stop ORF after that annotation. It exists as
an executable regression fixture rather than as an empty filename placeholder.

```bash
storf augment --input genome.gff3 --fasta genome.fna --outdir results --overlap-mode none
storf validate --gff results/genome/genome_StORF-Reporter.gff --fasta genome.fna
```

Expected summary: one sequence, one UR, and at least one retained positive-strand
candidate. Exact evidence scores are versioned in the automated tests.
