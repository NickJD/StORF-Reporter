# StORF-Reporter - Preprint: https://www.biorxiv.org/content/10.1101/2022.03.31.486628v1

Two Python3 scripts which extracts unannotated regions and then (St)ORFs using Stop-Stop codons.

# NO Need to clone repo (Currently contains testing data) - Just download UR_Extractor.py and StORF-Finder.py for now.

To run, you need: (/genomes contain all needed files for demo)  
* Input FASTA file  
* Input GFF file with predicted ORFs   
* Python 3 and the two scripts; IR_Extractor.py and StORF.py

# UR_Extractor.py
Python3 script to extract Unannotated Regions from DNA sequences uses FASTA and GFF files as input.

For Help: python3 UR_Extractor.py -h  
Example: python3 UR_Extractor.py -f genomes/E-coli.fasta.gz -gff genomes/E-coli.gff -o genomes/E-coli_UR -gz True
```python
usage: UR_Extractor.py [-h] -f FASTA -gff GFF [-ident IDENT] [-min_len MINLEN]
                       [-max_len MAXLEN] [-ex_len EXLEN]
                       [-gene_ident GENE_IDENT] -o OUT_PREFIX
                       [-gz {True,False}]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta_seq FASTA
                        FASTA file for Unannotated Region seq extraction
  -gff GFF              GFF annotation file for the FASTA
  -ident IDENT          Identifier given for Unannotated Region output
                        sequences: Default "Input"_UR
  -min_len MINLEN       Minimum UR Length: Default 30
  -max_len MAXLEN       Maximum UR Length: Default 100,000
  -ex_len EXLEN         UR Extension Length: Default 50
  -gene_ident GENE_IDENT
                        Identifier used for extraction of "genic" regions
                        "CDS,rRNA,tRNA": Default for Ensembl_Bacteria =
                        "ID=gene"
  -o OUT_PREFIX, --output_prefix OUT_PREFIX
                        Output file prefix - Without filetype
  -gz {True,False}      Default - False: Output as .gz
```
# StORF-Finder.py
Python3 script to extract Stop - Stop Codon (St)ORFs from Fasta sequences.  

For Help: python3 StORF_Finder.py -h  
Example: python3 StORF_Finder.py -seq genomes/E-coli_UR.fasta.gz -o genomes/E-coli_UR_StORF -gz True
```python
usage: StORF_Finder.py [-h] -f FASTA [-ua {True,False}] [-wc {True,False}]
                       [-ps {True,False}] [-filt [{none,soft,hard}]]
                       [-aa {True,False}] [-con_storfs {True,False}]
                       [-aa_only {True,False}] [-con_only {True,False}]
                       [-stop_ident {True,False}] [-minorf MIN_ORF]
                       [-maxorf MAX_ORF] [-codons STOP_CODONS]
                       [-olap OVERLAP_NT] [-gff {True,False}] [-o OUT_PREFIX]
                       [-lw {True,False}] [-gz {True,False}] [-v {True,False}]

StORF Run Parameters.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA              Input FASTA File
  -ua {True,False}      Default - Treat input as Unannotated: Use "-ua False"
                        for standard fasta
  -wc {True,False}      Default - False: StORFs reported across entire
                        sequence
  -ps {True,False}      Default - False: Partial StORFs reported
  -filt [{none,soft,hard}]
                        Default - Hard: Filtering level none is not
                        recommended, soft for single strand filtering and hard
                        for both-strand longest-first tiling
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -aa_only {True,False}
                        Default - False: Only output Amino Acid Fasta
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -stop_ident {True,False}
                        Default - True: Identify Stop Codon positions with '*'
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 50kb: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ('TAG,TGA,TAA'): List Stop Codons to use
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which
                        can overlap another StORF.
  -gff {True,False}     Default - True: StORF Output a GFF file
  -o OUT_PREFIX         Default - False/Same as input name with '_StORF-R':
                        Output filename prefix - Without filetype
  -lw {True,False}      Default - False: Line wrap FASTA sequence output at 60
                        chars
  -gz {True,False}      Default - False: Output as .gz
  -v {True,False}       Default - False: Print out runtime status
```
