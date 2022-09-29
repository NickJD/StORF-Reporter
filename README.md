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
usage: UR-Extractor [-h] -f FASTA -gff GFF [-ident IDENT] [-min_len MINLEN]
                       [-max_len MAXLEN] [-ex_len EXLEN]
                       [-gene_ident GENE_IDENT] [-o OUT_FILE]
                       [-gz {True,False}] [-v {True,False}]

StORF-Reporter v0.5.55: UR-Extractor Run Parameters.

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
                        Identifier used for extraction of "unannotated"
                        regions "CDS,rRNA,tRNA": Default for Ensembl_Bacteria
                        = "ID=gene"
  -o OUT_FILE, --output_file OUT_FILE
                        Output file - Without filetype - default appends "_UR"
                        to end of input gff filename (replaces '.gff')
  -gz {True,False}      Default - False: Output as .gz
  -v {True,False}       Default - False: Print out runtime status



```
## StORF-Finder
### Subpackage to extract Stop - Stop Codon (St)ORFs from Fasta sequences - Worked directly with the output of UR-Extractor.  

### Menu - (StORF-Finder -h):   
```console
StORF-Finder -seq genomes/E-coli_UR.fasta.gz -o genomes/E-coli_UR_StORF -gz True
```

```python
usage: StORF-Finder [-h] -f FASTA [-ua {True,False}] [-wc {True,False}]
                       [-ps {True,False}] [-filt [{none,soft,hard}]]
                       [-aa {True,False}] [-aa_only {True,False}]
                       [-con_storfs {True,False}] [-con_only {True,False}]
                       [-stop_ident {True,False}] [-type [{StORF,CDS,ORF}]]
                       [-minorf MIN_ORF] [-maxorf MAX_ORF]
                       [-codons STOP_CODONS] [-olap OVERLAP_NT]
                       [-gff {True,False}] [-s SUFFIX]
                       [-so [{start_pos,strand}]] [-spos {True,False}]
                       [-o OUT_FILE] [-af AFFIX] [-lw {True,False}]
                       [-gz {True,False}] [-v {True,False}]

StORF-Reporter v0.5.55: StORF-Finder Run Parameters.

optional arguments:
  -h, --help            show this help message and exit
  -f FASTA              Input FASTA File
  -ua {True,False}      Default - Treat input as Unannotated: Use "-ua False"
                        for standard fasta
  -wc {True,False}      Default - False: StORFs reported across entire
                        sequence
  -ps {True,False}      Default - False: Partial StORFs reported
  -filt [{none,soft,hard}]
                        Default - "hard": Filtering level "none" is not
                        recommended, "soft" for single strand filtering and
                        hard for both-strand longest-first tiling
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -aa_only {True,False}
                        Default - False: Only output Amino Acid Fasta
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -stop_ident {True,False}
                        Default - True: Identify Stop Codon positions with '*'
  -type [{StORF,CDS,ORF}]
                        Default - "StORF": Which GFF feature type for StORFs
                        to be reported as in GFF
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 50kb: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ('TAG,TGA,TAA'): List Stop Codons to use
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which
                        can overlap another StORF.
  -gff {True,False}     Default - True: StORF Output a GFF file
  -s SUFFIX             Default - Do not append suffix to genome ID
  -so [{start_pos,strand}]
                        Default - Start Position: How should StORFs be ordered
                        when >1 reported in a single UR.
  -spos {True,False}    Default - False: Print out StORF positions inclusive
                        of first stop codon
  -o OUT_FILE           Default - False: Without filetype - default appends
                        '_StORF-R' to end of input gff filename (replaces
                        '.gff')
  -af AFFIX             Default - None: '-af Con-StORFs' can be used to append
                        an identifier to output filename to distinguish Con-
                        StORF from StORF runs)
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60
                        chars
  -gz {True,False}      Default - False: Output as .gz
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
usage: StORF-Reporter [-h] -anno [{PROKKA,Ensembl,Gene}] [-pd PROKKA_DIR]
                         [-p_gff PROKKA_GFFS] [-col GENOME_COLLECTION]
                         [-comb COMBINED_GFFS] [-spos {True,False}]
                         [-rs {True,False}] [-sout {True,False}]
                         [-con_storfs {True,False}] [-con_only {True,False}]
                         [-min_len MINLEN] [-max_len MAXLEN] [-ex_len EXLEN]
                         [-minorf MIN_ORF] [-type [{StORF,CDS,ORF}]]
                         [-olap OVERLAP_NT] [-ao ALLOWED_OVERLAP]
                         [-lw {True,False}] [-aa {True,False}]
                         [-gz {True,False}] [-v {True,False}]

StORF-Reporter v0.5.55: StORF-Reporter Run Parameters.

optional arguments:
  -h, --help            show this help message and exit
  -anno [{PROKKA,Ensembl,Gene}]
                        Default - PROKKA: Annotation type to be StORF-
                        Reported:Options: PROKKA =
                        "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR";Ensembl =
                        "ID=gene" ;Gene = "gene"
  -pd PROKKA_DIR, --PROKKA_dir PROKKA_DIR
                        PROKKA output directory to be used if PROKKA chosen -
                        Produces a new GFF and FASTA containing all Coding and
                        Non-Coding Seqs
  -p_gff PROKKA_GFFS, --PROKKA_GFFs PROKKA_GFFS
                        Provide directory contain GFFs to be StORFed - Only
                        produces modified GFFs
  -col GENOME_COLLECTION
                        Provide directories containing a collection of
                        matching FASTA and GFF files(list .fa then .gff
                        containing directories separated by commas -
                        ./FA,./GFF)
  -comb COMBINED_GFFS   Provide directory containing GFFs with sequences
                        combined into single file to be StORFed - Only
                        produces modified GFFs
  -spos {True,False}    Default - False: Print out StORF positions inclusive
                        of first stop codon
  -rs {True,False}      Default - True: Remove stop "*" from StORF amino acid
                        sequences
  -sout {True,False}    Default - False: Print out StORF sequences separately?
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -min_len MINLEN       Default - 30: Minimum UR Length
  -max_len MAXLEN       Default - 100,000: Maximum UR Length
  -ex_len EXLEN         Default - 50: UR Extension Length
  -minorf MIN_ORF       Default - 100: Minimum StORF size in nt
  -type [{StORF,CDS,ORF}]
                        Default - "CDS": Which GFF feature type for StORFs to
                        be reported as in GFF - "CDS" is probably needed for
                        use in tools such as Roary
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which
                        can overlap another StORF.
  -ao ALLOWED_OVERLAP   Default - 50 nt: Maximum overlap between a StORF and
                        an original gene.
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60
                        chars
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -gz {True,False}      Default - False: Output as .gz
  -v {True,False}       Default - False: Print out runtime status


```


