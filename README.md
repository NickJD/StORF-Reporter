# Intergenic-Region-Extractor and StORF

Two Python3 scripts which extracts intergenic regions and then (St)ORFs using Stop-Stop codons.

To run, you need: (/genomes contain all needed files for demo)  
* Input FASTA file  
* Input GFF file with predicted ORFs   
* Python 3 and the two scripts; IR_Extractor.py and StORF.py

# IR_Extractor.py
Python3 script to extract Intergenic Regions from DNA sequences uses FASTA and GFF files as input.

For Help: python3 IR_Extractor.py -h  
Example: python3 IR_Extractor.py -f genomes/E-coli.fasta.gz -gff genomes/E-coli.gff -o genomes/E-coli_IR -gz True
```python
usage: IR_Extractor.py [-h] -f FASTA -gff GFF [-ident IDENT] [-min_len MINLEN]
                       [-ex_len EXLEN] [-gene_ident GENE_IDENT] -o OUT_PREFIX
                       [-gz {True,False}]

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA, --fasta_seq FASTA
                        FASTA file for Intergenic Region seq extraction
  -gff GFF              GFF annotation file for the FASTA
  -ident IDENT          Identifier given for Intergenic Region output sequences: Default "Input"_IR
  -min_len MINLEN       Minimum IR Length: Default 30
  -ex_len EXLEN         IR Extension Length: Default 50
  -gene_ident GENE_IDENT
                        Identifier used for extraction: Default = ID=gene
  -o OUT_PREFIX, --output_prefix OUT_PREFIX
                        Output file prefix - Without filetype
  -gz {True,False}      Default - False: Output as .gz



```
# StORF.py
Python3 script to extract Stop - Stop Codon (St)ORFs from Fasta sequences.  

For Help: python3 StORF.py -h  
Example: python3 StORF.py -seq genomes/E-coli_IR.fasta.gz -o genomes/E-coli_IR_StORF -gz True
```python
usage: StORF.py [-h] -seq FASTA [-ir {True,False}] [-wc {True,False}] 
                    [-ps {True,False}] [-filt [{none,soft,hard}]] [-aa {True,False}] [-con_storfs {True,False}] 
                    [-aa_only {True,False}] [-con_only {True,False}] [-stop_ident {True,False}] [-minorf MIN_ORF] 
                    [-maxorf MAX_ORF] [-codons STOP_CODONS] [-olap OVERLAP_NT] [-gff {True,False}] -o OUT_PREFIX 
                    [-gz {True,False}]


optional arguments:
  -h, --help            show this help message and exit
  -seq FASTA            Input Sequence File
  -ir {True,False}      Default - Treat input as Intergenic: Use "-ir False" for standard fasta
  -wc {True,False}      Default - False: StORFs reported across entire sequence
  -ps {True,False}      Default - False: Partial StORFs reported
  -filt [{none,soft,hard}]
                        Default - Hard: Filtering level none is not recommended, soft for single strand filtering and hard for both-strand longest-first tiling
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -aa_only {True,False}
                        Default - False: Only output Amino Acid Fasta
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -stop_ident {True,False}
                        Default - True: Identify Stop Codon positions with "*"
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 99999: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ("TAG,TGA,TAA"): List Stop Codons to use
  -olap OVERLAP_NT      Default - 20: Maximum number of nt of a StORF which can overlap another StORF.
  -gff {True,False}     Default - True: StORF Output a GFF file
  -o OUT_PREFIX         Output file prefix - Without filetype
  -gz {True,False}      Default - False: Output as .gz




```
