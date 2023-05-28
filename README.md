# StORF-Reporter - Preprint: https://www.biorxiv.org/content/10.1101/2022.03.31.486628v3

### StORF-Reporter, a toolkit that returns missed CDS genes from the Unannotated Regions (URs) of prokaryotic genomes.

# Please use `pip3 install StORF-Reporter' to install StORF-Reporter v0.7.5.
### This will also install the python-standard library numpy (>=1.22.0,<1.24.0), Pyrodigal - (https://github.com/althonos/pyrodigal) and ORForise (https://github.com/NickJD/ORForise). 

### Consider using '--no-cache-dir' with pip to ensure the download of the newest version of StORF-Reporter.

## Please Note: To report Con-StORFs (Pseudogenes and genes that have alternative use of stop codons), use "-con_storfs True". To disable the reporting of StORFs use "-con_only". 

### The directory "Test_Datasets" is provided to confirm functionality of StORF-Reporter.

#############################################################
# StORF-Reporter:
## Most common use cases - 
### Supplement a current annotation from a tool such as Prokka or Bakta. A new GFF file will be created compatible with downstream pangenome analysis tools such as Roary and Panaroo.

#### For use on a single Prokka/Bakta output directory - Will also create a new fasta file with Prokka/Bakta gene and StORF sequences.
```console
StORF-Reporter -anno Prokka Out_Dir -p .../Test_Datasets/Prokka_E-coli/
```
#### For use on a directory containing multiple Prokka/Bakta output gffs - Only produces new GFF files. 
```console
StORF-Reporter -anno Prokka Multiple_GFFs -p .../Test_Datasets/Prokka_Outputs/
```

#### For use on a GFF file from a CDS prediction tool such as Prodigal - Provide a GFF file and StORF-Reporter will find the matching .fa/.fasta/.fna (must have the same name). 
```console
StORF-Reporter -anno Feature_Types Single_Genome -p .../Test_Datasets/Matching_GFF_FASTA/Myco.gff
```

#### For use on a directory containing multiple GFF files from a CDS prediction tool such as Prodigal - StORF-Reporter will find the matching .fa/.fasta/.fna (must have the same name). 
```console
StORF-Reporter -anno Feature_Types Multiple_Genomes -p .../Test_Datasets/Matching_GFF_FASTA/
```

#### For use on a directory containing multiple GFF files with embedded FASTA. 
```console
StORF-Reporter -anno Feature_Types Multiple_Combined_GFFs -p .../Test_Datasets/Combined_GFFs/
```

#### To perform a fresh end-to-end annotation of a genome without an annotation, StORF-Reporter will use Pyrodigal to predict CDS genes and then supplement with StORFs. 
```console
StORF-Reporter -anno Pyrodigal Single_FASTA -p .../Test_Datasets/Pyrodigal/E-coli.fa
```

### Menu - (StORF-Reporter -h):
```python
usage: StORF_Reporter.py [-h]
                         [-anno [{Prokka,Bakta,Out_Dir,Single_GFF,Multiple_GFFs,Ensembl,Feature_Types,Single_Genome,Multiple_Genomes,Single_Combined_GFF,Multiple_Combined_GFFs,Pyrodigal,Single_FASTA,Multiple_FASTA} ...]]
                         [-p PATH] [-oname O_NAME] [-odir O_DIR] [-sout {True,False}] [-lw {True,False}] [-aa {True,False}]
                         [-gz {True,False}] [-py_train [{longest,individual,meta}]] [-py_fasta {True,False}]
                         [-py_unstorfed {True,False}] [-gene_ident GENE_IDENT] [-min_len MINLEN] [-max_len MAXLEN]
                         [-ex_len EXLEN] [-spos {True,False}] [-rs {True,False}] [-con_storfs {True,False}]
                         [-con_only {True,False}] [-ps {True,False}] [-wc {True,False}] [-short_storfs {False,Nolap,Olap}]
                         [-short_storfs_only {True,False}] [-minorf MIN_ORF] [-maxorf MAX_ORF] [-codons STOP_CODONS]
                         [-olap_filt [{none,single-strand,both-strand}]] [-start_filt {True,False}] [-so [{start_pos,strand}]]
                         [-f_type [{StORF,CDS,ORF}]] [-olap OVERLAP_NT] [-ao ALLOWED_OVERLAP] [-overwrite {True,False}]
                         [-verbose {True,False}] [-v]

StORF-Reporter v0.7.5: StORF-Reporter Run Parameters.

Required Options:
  -anno [{Prokka,Bakta,Out_Dir,Single_GFF,Multiple_GFFs,Ensembl,Feature_Types,Single_Genome,Multiple_Genomes,Single_Combined_GFF,Multiple_Combined_GFFs,Pyrodigal,Single_FASTA,Multiple_FASTA} ...]
                        Select Annotation and Input options for one of the 3 options listed below
                        ### Prokka/Bakta Annotation Option 1: 
                        	Prokka = Report StORFs for a Prokka annotation; 
                        	Bakta = Report StORFs for a Bakta annotation; 
                        --- Prokka/Bakta Input Options: 
                        	Out_Dir = To provide the output directory of either a Prokka or Bakta run (will produce a new GFF and FASTA file); 
                        	Single_GFF = To provide a single Prokka or Bakta GFF file (will not provide new FASTA file); 
                        	Multiple_GFFs = To provide a directory containing multiple GFF files in Prokka/Bakta format (will not provide a new FASTA file); 
                        
                        ### Standard GFF Annotation Option 2: 
                        	Ensembl = Report StORFs for an Ensembl Bacteria annotation (ID=gene); 
                        	Feature_Types = Used in conjunction with -gene_ident to define features such as CDS,rRNA,tRNA for UR extraction (default CDS); 
                        --- Standard GFF Input Options: 
                        	Single_Genome = To provide a single Genome - accompanying FASTA must share same name as given gff file (can be .fna, fa or .fasta); 
                        	Multiple_Genomes = To provide a directory containing multiple accompanying GFF and FASTA files - files must share the same name (fasta can be .fna, fa or .fasta); 
                        	Single_Combined_GFF = To provide a GFF file with embedded FASTA at the bottom; 
                        	Multiple_Combined_GFFs = To provide a directory containing multiple GFF files with embedded FASTA at the bottom; 
                        
                        ### Complete Annotation Option 3: 
                        	Pyrodigal = Run Pyrodigal then Report StORFs (provide path to single FASTA or directory of multiple FASTA files ;
                        --- Complete Annotation Input Options: 
                        	Single_FASTA = To provide a single FASTA file; 
                        	Multiple_FASTA = To provide a directory containing multiple FASTA files (will detect .fna,.fa,.fasta); 
                        
  -p PATH               Provide input file or directory path

StORF-Reporter Options:
  -oname O_NAME         Default - Appends '_StORF-R' to end of input FASTA filename - Multiple_* runs will be numbered
  -odir O_DIR           Default - Same directory as input FASTA
  -sout {True,False}    Default - False: Print out StORF sequences separately from Prokka/Bakta annotations
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60 chars
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -gz {True,False}      Default - False: Output as .gz

Pyrodigal Options:
  -py_train [{longest,individual,meta}]
                        Default - longest: Type of model training to be done for Pyrodigal CDS prediction: Options: longest =
                        Trains on longest contig; individual = Trains on each contig separately - runs in meta mode if contig is
                        < 20KB; meta = Runs in meta mode for all sequences
  -py_fasta {True,False}
                        Default - False: Output Pyrodigal+StORF predictions in FASTA format
  -py_unstorfed {True,False}
                        Default - False: Provide GFF containing original Pyrodigal predictions

UR-Extractor Options:
  -gene_ident GENE_IDENT
                        Identifier used for extraction of Unannotated Regions "CDS,rRNA,tRNA" - To be used with "-anno
                        Feature_Types"
  -min_len MINLEN       Default - 30: Minimum UR Length
  -max_len MAXLEN       Default - 100,000: Maximum UR Length
  -ex_len EXLEN         Default - 50: UR Extension Length

StORF-Finder Options:
  -spos {True,False}    Default - False: Output StORF positions inclusive of first stop codon
  -rs {True,False}      Default - True: Remove stop "*" from StORF amino acid sequences
  -con_storfs {True,False}
                        Default - False: Output Consecutive StORFs
  -con_only {True,False}
                        Default - False: Only output Consecutive StORFs
  -ps {True,False}      Default - False: Partial StORFs reported
  -wc {True,False}      Default - False: StORFs reported across entire sequence
  -short_storfs {False,Nolap,Olap}
                        Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt
                        that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs
                        which areoverlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs. Overlap is
                        defined by "-olap".
  -short_storfs_only {True,False}
                        Default - True. Only report Short-StORFs?
  -minorf MIN_ORF       Default - 99: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 60kb: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ('TAG,TGA,TAA'): List Stop Codons to use
  -olap_filt [{none,single-strand,both-strand}]
                        Default - "both-strand": Filtering level "none" is not recommended, "single-strand" for single strand
                        filtering and both-strand for both-strand longest-first tiling
  -start_filt {True,False}
                        Default - False: Filter out StORFs without at least one of the 3 common start codons (best used for
                        short-storfs).
  -so [{start_pos,strand}]
                        Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.
  -f_type [{StORF,CDS,ORF}]
                        Default - "CDS": Which GFF feature type for StORFs to be reported as in GFF - "CDS" is probably needed
                        for use in tools such as Roary and Panaroo
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which can overlap another StORF.
  -ao ALLOWED_OVERLAP   Default - 50 nt: Maximum overlap between a StORF and an original gene.

Misc:
  -overwrite {True,False}
                        Default - False: Overwrite StORF-Reporter output if already present
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit
```

###################################

# UR-Extractor:
### Subpackage to extract Unannotated Regions from DNA sequences using FASTA and GFF files as input.

### Menu - (UR-Extractor -h):  
```console
UR-Extractor -f .../Test_Datasets/Matching_GFF_FASTA/E-coli.fa -gff .../Test_Datasets/Matching_GFF_FASTA/E-coli.gff
```

```python
usage: UR_Extractor.py [-h] [-f FASTA] [-gff GFF] [-ident IDENT] [-min_len MINLEN] [-max_len MAXLEN] [-ex_len EXLEN] [-gene_ident GENE_IDENT] [-oname O_NAME] [-odir O_DIR] [-gz {True,False}] [-verbose {True,False}] [-v]

StORF-Reporter v0.7.5: UR-Extractor Run Parameters.

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
  -oname O_NAME         Default - Appends '_UR' to end of input GFF filename
  -odir O_DIR           Default - Same directory as input GFF
  -gz {True,False}      Default - False: Output as .gz

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit

```
## StORF-Finder:
### Subpackage to extract StORFs from Fasta sequences - Works directly with the output of UR-Extractor.  

### Menu - (StORF-Finder -h):   
```console
StORF-Finder -f .../Test_Datasets/Matching_GFF_FASTA/E-coli_UR.fa 
```

```python
usage: StORF_Finder.py [-h] [-f FASTA] [-ua {True,False}] [-wc {True,False}] [-ps {True,False}] [-olap_filt [{none,single-strand,both-strand}]] [-start_filt {True,False}] [-con_storfs {True,False}] [-con_only {True,False}] [-short_storfs {False,Nolap,Olap}] [-short_storfs_only {True,False}]
                       [-stop_ident {True,False}] [-f_type [{StORF,CDS,ORF}]] [-minorf MIN_ORF] [-maxorf MAX_ORF] [-codons STOP_CODONS] [-olap OVERLAP_NT] [-s SUFFIX] [-so [{start_pos,strand}]] [-spos {True,False}] [-oname O_NAME] [-odir O_DIR] [-gff {True,False}] [-aa {True,False}] [-aa_only {True,False}]
                       [-lw {True,False}] [-gff_fasta {True,False}] [-gz {True,False}] [-verbose {True,False}] [-v]

StORF-Reporter v0.7.5: StORF-Finder Run Parameters.

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
                        Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs which areoverlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs.
                        Overlap is defined by "-olap".
  -short_storfs_only {True,False}
                        Default - True. Only report Short-StORFs?
  -stop_ident {True,False}
                        Default - True: Identify Stop Codon positions with '*'
  -f_type [{StORF,CDS,ORF}]
                        Default - "StORF": Which GFF feature type for StORFs to be reported as in GFF
  -minorf MIN_ORF       Default - 99: Minimum StORF size in nt
  -maxorf MAX_ORF       Default - 60kb: Maximum StORF size in nt
  -codons STOP_CODONS   Default - ('TAG,TGA,TAA'): List Stop Codons to use
  -olap OVERLAP_NT      Default - 50: Maximum number of nt of a StORF which can overlap another StORF.
  -s SUFFIX             Default - Do not append suffix to genome ID
  -so [{start_pos,strand}]
                        Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.
  -spos {True,False}    Default - False: Print out StORF positions inclusive of first stop codon

Output:
  -oname O_NAME         Default - Appends '_StORF-R' to end of input FASTA filename
  -odir O_DIR           Default - Same directory as input FASTA
  -gff {True,False}     Default - True: Output a GFF file
  -aa {True,False}      Default - False: Report StORFs as amino acid sequences
  -aa_only {True,False}
                        Default - False: Only output Amino Acid Fasta
  -lw {True,False}      Default - True: Line wrap FASTA sequence output at 60 chars
  -gff_fasta {True,False}
                        Default - False: Report all gene sequences (nt) at the bottom of GFF files in Prokka output mode
  -gz {True,False}      Default - False: Output as .gz

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit

```
## StORF-Extractor
Subpackage to extract sequences reported by StORF-Reporter from a genome annotation.

### Menu - (StORF-Extractor -h):   
```console
StORF-Extractor -storf_input Combined -p .../Test_Datasets/Combined_GFFs/E-coli_Combined_StORF-Reporter_Extended.gff 
```

```python
usage: StORF_Extractor.py [-h] [-storf_input {Combined,Separate}] [-p PATH] [-gff_out {True,False}] [-oname O_NAME] [-odir O_DIR] [-gz {True,False}] [-verbose {True,False}] [-v]

StORF-Reporter v0.7.5: StORF-Extractor Run Parameters.

Required Arguments:
  -storf_input {Combined,Separate}
                        Are StORFs to be extracted from Combined GFF/FASTA or Separate GFF/FASTA files?
  -p PATH               Provide input file or directory path

Output:
  -gff_out {True,False}
                        Default - False: Output StORFs in GFF format
  -oname O_NAME         Default - Appends '_Extracted_StORFs' to end of input GFF filename
  -odir O_DIR           Default - Same directory as input FASTA
  -gz {True,False}      Default - False: Output as .gz

Misc:
  -verbose {True,False}
                        Default - False: Print out runtime messages
  -v                    Default - False: Print out version number and exit

```

## Test Datasets: 
### The directory 'Test_Datasets' contains GFF and FASTA files to test the installation and use of StORF-Reporter. 