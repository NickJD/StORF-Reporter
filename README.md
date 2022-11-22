# StORF-Reporter - Preprint: https://www.biorxiv.org/content/10.1101/2022.03.31.486628v3

### StORF-Reporter, a toolkit that takes as input an annotated genome and returns missed CDS genes from the Unannotated Regions (URs).

# Please use `pip3 install StORF-Reporter' to install StORF-Reporter.
### This will also install numpy and the 'ORForise' package from https://github.com/NickJD/ORForise to allow for additional functionality. 

### Consider using '--no-cache-dir' with pip to ensure the download of the newest version of the package.

## Please Note: To report Con-StORFs (Pseudogenes and genes that reuse stop codons), use "-con_storfs True" and "-con_only" to disable reporting of StORFs (missed complete genes). 

#############################################################

# UR-Extractor
### Subpackage to extract Unannotated Regions from DNA sequences using FASTA and GFF files as input.

### Menu - (UR-Extractor -h):  
```console
UR-Extractor -f genomes/E-coli.fasta.gz -gff genomes/E-coli.gff -o genomes/E-coli_UR -gz True
```

```python
usage: UR_Extractor.py [-h] -f FASTA -gff GFF [-ident IDENT] [-min_len MINLEN] [-max_len MAXLEN] [-ex_len EXLEN] [-gene_ident GENE_IDENT] [-o OUT_FILE] [-gz {True,False}] [-v {True,False}]

StORF-Reporter v0.6.0: UR-Extractor Run Parameters.

Required Arguments:
  -f FASTA              FASTA file for Unannotated Region seq extraction
  -gff GFF              GFF annotation file for the FASTA

Optional Arguments:
  -ident IDENT          Identifier given for Unannotated Region output sequences - Do not modify if output is to be used by StORF-Finder: Default "Sequence-ID"_UR
  -min_len MINLEN       Minimum UR Length: Default 30
  -max_len MAXLEN       Maximum UR Length: Default 100,000
  -ex_len EXLEN         UR Extension Length on 5' and 3': Default 50
  -gene_ident GENE_IDENT
                        Identifier used for extraction of Unannotated Regions "CDS,rRNA,tRNA": Default for Ensembl_Bacteria = "ID=gene" or "-gene_ident CDS" for "most" genome annotations

Output:
  -o OUT_FILE           Output file - Without filetype - default appends "_UR" to end of input gff filename (replaces '.gff')
  -gz {True,False}      Default - False: Output as .gz

Misc:
  -v {True,False}       Default - False: Print out runtime status

```
## StORF-Finder
### Subpackage to extract Stop - Stop Codon (St)ORFs from Fasta sequences - Worked directly with the output of UR-Extractor.  

### Menu - (StORF-Finder -h):   
```console
StORF-Finder -f genomes/E-coli_UR.fasta.gz -o genomes/E-coli_UR_StORF 
```

```python
usage: StORF_Finder.py [-h] -f FASTA [-ua {True,False}] [-wc {True,False}] [-ps {True,False}] [-olap_filt [{none,single-strand,both-strand}]] [-start_filt {True,False}] [-con_storfs {True,False}] [-con_only {True,False}]
                       [-short_storfs {False,Nolap,Olap}] [-short_storfs_only {True,False}] [-stop_ident {True,False}] [-type [{StORF,CDS,ORF}]] [-minorf MIN_ORF] [-maxorf MAX_ORF] [-codons STOP_CODONS] [-olap OVERLAP_NT] [-s SUFFIX]
                       [-so [{start_pos,strand}]] [-spos {True,False}] [-o OUT_FILE] [-af AFFIX] [-gff {True,False}] [-aa {True,False}] [-aa_only {True,False}] [-lw {True,False}] [-gff_fasta {True,False}] [-gz {True,False}]
                       [-v {True,False}]

StORF-Reporter v0.6.0: StORF-Finder Run Parameters.

Required Arguments:
  -f FASTA              Input FASTA File - (UR_Extractor output)

Optional Arguments:
  -ua {True,False}      Default - Treat input as Unannotated: Use "-ua False" for standard fasta
  -wc {True,False}      Default - False: StORFs reported across entire sequence
  -ps {True,False}      Default - False: Partial StORFs reported
  -olap_filt [{none,single-strand,both-strand}]
                        Default - "both-strand": Filtering level "none" is not recommended, "single-strand" for single strand filtering and both-strand for both-strand longest-first tiling
  -start_filt {True,False}
                        Default - False: Filter out StORFs without at least one of the 3 common start codons (best used for short-storfs).
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -short_storfs {False,Nolap,Olap}
                        Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs which
                        areoverlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs. Overlap is defined by "-olap".
  -short_storfs_only {True,False}
                        Default - True. Only report Short-StORFs?
  -stop_ident {True,False}
                        Default - True: Identify Stop Codon positions with '*'
  -type [{StORF,CDS,ORF}]
                        Default - "StORF": Which GFF feature type for StORFs to be reported as in GFF
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 50kb: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ('TAG,TGA,TAA'): List Stop Codons to use
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which can overlap another StORF.
  -s SUFFIX             Default - Do not append suffix to genome ID
  -so [{start_pos,strand}]
                        Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.
  -spos {True,False}    Default - False: Print out StORF positions inclusive of first stop codon
  -o OUT_FILE           Default - False: Without filetype - default appends '_StORF-R' to end of input gff filename (replaces '.gff')
  -af AFFIX             Default - None: '-af Con-StORFs' can be used to append an identifier to output filename to distinguish Con-StORF from StORF runs)

Output:
  -gff {True,False}     Default - True: Output a GFF file
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -aa_only {True,False}
                        Default - False: Only output Amino Acid Fasta
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60 chars
  -gff_fasta {True,False}
                        Default - False: Report all gene sequences (nt) at the bottom of GFF files in PROKKA output mode
  -gz {True,False}      Default - False: Output as .gz

Misc:
  -v {True,False}       Default - False: Print out runtime status

```

##########
## StORF-Reporter

### Subpackage to run UR-Extractor and StORF-Finder together on a PROKKA annotation and produce output which can be used in tools such as Roary that by default accepts PROKKA's output format: 
#### For use on a single PROKKA output directory - 
```console
StORF-Reporter -anno PROKKA -pd ../PROKKA_04062022/
```
#### For use on a directory containing multiple PROKKA output gffs
```console
StORF-Reporter -anno PROKKA -p_gff ../PROKKA_Outputs/
```

### Menu - (StORF-Reporter -h):
```python
usage: StORF_Reporter.py [-h] -anno [{PROKKA,Ensembl,Gene}] -dt [{PROKKA_Out,PROKKA_GFFs}] -pd PROKKA_DIR [-comb COMBINED_GFFS] [-spos {True,False}] [-rs {True,False}] [-sout {True,False}] [-con_storfs {True,False}]
                         [-con_only {True,False}] [-short_storfs {False,Nolap,Olap}] [-short_storfs_only {True,False}] [-min_len MINLEN] [-max_len MAXLEN] [-ex_len EXLEN] [-minorf MIN_ORF] [-olap_filt [{none,single-strand,both-strand}]]
                         [-start_filt {True,False}] [-type [{StORF,CDS,ORF}]] [-olap OVERLAP_NT] [-ao ALLOWED_OVERLAP] [-lw {True,False}] [-aa {True,False}] [-gff_fasta {True,False}] [-gz {True,False}] [-v {True,False}]

StORF-Reporter v0.6.0: StORF-Reporter Run Parameters.

Required Arguments:
  -anno [{PROKKA,Ensembl,Gene}]
                        Annotation type to be StORF-Reported - Options: PROKKA = "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR";Ensembl = "ID=gene" ;Gene = "gene"
  -dt [{PROKKA_Out,PROKKA_GFFs}]
                        Type of directory to be loaded by StORF-Reported: Options: PROKKA_Out = Single PROKKA output directory; PROKKA_GFFs = Directory containing multiple PROKKA GFF files
  -pd PROKKA_DIR        PROKKA output directory to be used - Produces a new GFF and FASTA containing all Coding and Non-Coding Seqs for each PROKKA output GFF

Optional Arguments:
  -comb COMBINED_GFFS   Provide directory containing GFFs with sequences combined into single file to be StORFed - Only produces modified GFFs
  -spos {True,False}    Default - False: Print out StORF positions inclusive of first stop codon
  -rs {True,False}      Default - True: Remove stop "*" from StORF amino acid sequences
  -sout {True,False}    Default - False: Print out StORF sequences separately?
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -short_storfs {False,Nolap,Olap}
                        Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs which
                        areoverlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs. Overlap is defined by "-olap".
  -short_storfs_only {True,False}
                        Default - True. Only report Short-StORFs?
  -min_len MINLEN       Default - 30: Minimum UR Length
  -max_len MAXLEN       Default - 100,000: Maximum UR Length
  -ex_len EXLEN         Default - 50: UR Extension Length
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -olap_filt [{none,single-strand,both-strand}]
                        Default - "both-strand": Filtering level "none" is not recommended, "single-strand" for single strand filtering and both-strand for both-strand longest-first tiling
  -start_filt {True,False}
                        Default - False: Filter out StORFs without at least one of the 3 common start codons (best used for short-storfs).
  -type [{StORF,CDS,ORF}]
                        Default - "CDS": Which GFF feature type for StORFs to be reported as in GFF - "CDS" is probably needed for use in tools such as Roary
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which can overlap another StORF.
  -ao ALLOWED_OVERLAP   Default - 50 nt: Maximum overlap between a StORF and an original gene.

Output::
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60 chars
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -gff_fasta {True,False}
                        Default - False: Report all gene sequences (nt) at the bottom of GFF files in PROKKA output mode
  -gz {True,False}      Default - False: Output as .gz

Misc::
  -v {True,False}       Default - False: Print out runtime status

```


## Test Datasets: 
### The direcotory 'Genomes' contains GFF and FASTA datasets to test a users installation and use of StORF-Reporter. Example output files are also available for comparison.