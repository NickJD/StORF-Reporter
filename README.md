# Intergenic-Region-Extractor and StORF
Python3 script to extract Intergenic Regions from DNA sequences uses FASTA and GFF files as input.

For Help: python3 IR_Extractor.py -h  
Example: python2 IR_Extractor.py -f genomes/E-coli.fa -gff genomes/E-coli.gff -o genomes/E-coli_IR
```python
usage: IR_Extractor.py [-h] -f FASTA -gff GFF [-ident IDENT] [-min_len MINLEN]
                       [-ex_len EXLEN] [-gene_ident GENE_IDENT] -o OUT_PREFIX

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta_seq FASTA
                        FASTA file for Intergenic Region seq extraction
  -gff GFF              GFF annotation file for the FASTA
  -ident IDENT          Identifier given for Intergenic Region output
                        sequences: Default Input_IR
  -min_len MINLEN       Minimum IR Length: Default 30
  -ex_len EXLEN         IR Extension Length: Default 50
  -gene_ident GENE_IDENT
                        Identifier used for extraction: Default = ID=gene:
  -o OUT_PREFIX, --output_prefix OUT_PREFIX
                        Output file prefix - Without filetype

```

Python3 script to detect Stop - Stop Codon ORFs from Fasta sequences.  
For Help: python3 StORF.py -h  
Example: python3 StORF.py -seq E-coli_IR.fa -o E-coli_IR_Storf
```python
usage: StORF.py [-h] -seq FASTA [-ir {True,False}] [-wc {True,False}]
                [-ps {True,False}] [-filt [{none,soft,hard}]]
                [-aa {True,False}] [-minorf MIN_ORF] [-maxorf MAX_ORF]
                [-codons STOP_CODONS] [-olap OVERLAP_NT] [-gff {True,False}]
                -o OUTPUT

optional arguments:
  -h, --help            show this help message and exit
  -seq FASTA            Input Sequence File
  -ir {True,False}      Default - Treat input as Intergenic: Use "-ir False"
                        for standard fasta
  -wc {True,False}      Default - False: StORFs reported across entire
                        sequence
  -ps {True,False}      Default - False: Partial StORFs reported
  -filt [{none,soft,hard}]
                        Default - Hard: Filtering level none is not
                        recommended, soft for single strand filtering and hard
                        for both-strand longest-first tiling
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 99999: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ("TAG,TGA,TAA"): List Stop Codons to use
  -olap OVERLAP_NT      Default - 20: Maximum number of nt of a StORF which
                        can overlap another StORF.
  -gff {True,False}     Default - True: StORF Output a GFF file
  -o OUTPUT             Output file prefix - Without filetype

```
