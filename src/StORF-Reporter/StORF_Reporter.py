import argparse
from argparse import Namespace
import pathlib
import collections



try:
    from UR_Extractor import extractor
    from StORF_Finder import StORF_Reported
except ImportError:
    from .UR_Extractor import extractor
    from .StORF_Finder import StORF_Reported


###################
gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}
############################

def GFF_StoRF_write(StORF_options,track_contig,gff_outfile, StORF,StORF_Num): # Consistency in outfile
    ID = StORF_options.gff.split('/')[-1].replace('.gff','') + '_StORF_'+str(StORF_Num)
    #ID = track_contig + '_Stop-ORF:' + str(StORF[1]) + '-' + str(StORF[2])
    ### Write out new GFF entry -
    gff_outfile.write(track_contig + '\tStORF_Reporter\t' + StORF_options.feature_type + '\t' +  str(StORF[1]) + '\t' + str(StORF[2]) + '\t.\t' +
        StORF[7] + '\t.\tID=' + ID + ';locus_tag=' + ID + ';INFO=Additional_Annotation_StORF-Reporter;UR_Position=' + StORF[0] + ';Name=StORF_' + str(StORF_Num) + ';StORF_Num=' + str(
            StORF_Num) + ';StORF_Num_In_UR=' + str(StORF[9]) +
        ';StORF_Length=' + str(StORF[8]) + ';StORF_Frame=' + str(StORF[5]) + ';UR_StORF_Frame=' + str(
            StORF[6]) + '\n')

def FASTA_StoRF_write(track_contig, fasta_outfile, StORF):  # Consistency in outfile
    ID = track_contig + '_Stop-ORF:' + str(StORF[1]) + '-' + str(StORF[2])
    ### Wrtie out new FASTA entry
    fasta_outfile.write('>'+ID+'\n'+StORF[10]+'\n')

def FASTA_Load(faa_infile,ffn_infile):
    ### Coding sequences first
    PROKKA_CoDing = collections.defaultdict()
    first = True
    infile = open(faa_infile)
    seq = ''
    ##load in fasta file and make dict of seqs and ids
    for line in infile:
        line = line.strip()
        if line.startswith('>') and first == False:
            PROKKA_CoDing.update({id:seq})
            id = line.replace('>','')
            id = id.split(' ')[0]
            seq = ''
        elif line.startswith('>'):
            id = line.replace('>','')
            id = id.split(' ')[0]
            first = False
        else:
            seq += line
    PROKKA_CoDing.update({id:seq})
    ### Non-Coding sequences next
    PROKKA_Non_CoDing = collections.defaultdict()
    first = True
    infile = open(ffn_infile)
    seq = ''
    ##load in fasta file and make dict of seqs and ids
    for line in infile:
        line = line.strip()
        if line.startswith('>') and first == False:
            PROKKA_Non_CoDing.update({id:seq})
            id = line.replace('>','')
            id = id.split(' ')[0]
            seq = ''
        elif line.startswith('>'):
            id = line.replace('>','')
            id = id.split(' ')[0]
            first = False
        else:
            seq += line
    PROKKA_Non_CoDing.update({id:seq})

    return PROKKA_CoDing,PROKKA_Non_CoDing


def run_UR_Extractor_PROKKA_DIR(Reporter_options): # When given a single PROKKA Directory
    fasta = list(pathlib.Path(Reporter_options.prokka_dir).glob('*.fna'))
    fasta = str(fasta[0])
    Reporter_options.fasta = fasta
    gff_list = list(pathlib.Path(Reporter_options.prokka_dir).glob('*.gff'))
    for gff in gff_list: # Poor checking method - must change
        if '_StORF_Reporter' not in str(gff): # Might fall over
            gff = str(gff)
            Reporter_options.gff = gff
    URs = extractor(Reporter_options)
    return URs,Reporter_options


def run_UR_Extractor_PROKKA_GFFs(Reporter_options,gff): # When given a directory with multiple PROKKA GFFs but without accompianing .fna
    if '_StORF_Reporter' not in str(gff): #Might fall over
        gff = str(gff)
        Reporter_options.gff = gff
    Reporter_options.fasta = gff
    URs = extractor(Reporter_options)
    return URs,Reporter_options

def find_prev_StORFs(StORF_options,Contig_URs,track_current_start,track_prev_stop, track_contig):
    StORFs_to_del = []
    StORFs = []
    is_break = False
    current_Contig_URs = Contig_URs[track_contig]
    for StORF_Num, data in current_Contig_URs.items():
        if is_break == True:
            break
        if data != None:
         #   for UR_StORF, UR_StORF_data in UR_StORFs.items():
            ur_pos = data[7]
            StORF_Pos_In_UR = data[0]
            StORF_Start_In_UR = int(StORF_Pos_In_UR.split(',')[0])
            StORF_Stop_In_UR = int(StORF_Pos_In_UR.split(',')[1])
            ur_frame = data[2]
            strand = data[3]
            key_start = int(ur_pos.split('_')[0])
            StORF_start = key_start + StORF_Start_In_UR
            key_stop = int(ur_pos.split('_')[1])
            StORF_stop = key_start + StORF_Stop_In_UR
            allow_start = StORF_start + StORF_options.allowed_overlap
            allow_stop = StORF_stop - StORF_options.allowed_overlap
            StORF_Seq = data[1]
            StORF_Length = data[4]
            StORF_UR_Num = data[6]
            if allow_stop <= track_current_start:# and allow_start >= track_prev_stop:
                gff_start = str(StORF_Start_In_UR + int(ur_pos.split('_')[-1].split('_')[0])) #Here - check
                gff_stop = str(StORF_Stop_In_UR + int(ur_pos.split('_')[-1].split('_')[0]))
                if strand == '+': # Get original genome frame
                    frame = (int(gff_start) % 3) + 1
                elif strand == '-':
                    frame = (int(gff_stop) % 3) + 4
                StORFs.append([ur_pos,StORF_start,StORF_stop,StORF_Start_In_UR,StORF_Stop_In_UR,frame,ur_frame,strand,
                               StORF_Length,StORF_UR_Num,StORF_Seq])
                if StORF_Num not in StORFs_to_del:
                    StORFs_to_del.append(StORF_Num)
            elif allow_start > track_current_start: # Check
                is_break = True
                break

    for StORF in StORFs_to_del:
        del current_Contig_URs[StORF]
    Contig_URs[track_contig] = current_Contig_URs # Just incase shallow copy does not remove StORF for us
### Might need to force remove more ur storfs regions

    return StORFs, Contig_URs

def find_after_StORFs(options,Contig_URs,track_current_start,track_current_stop, track_contig):
    current_Contig_URs = Contig_URs[track_contig]
    StORFs = []
    is_break = False
    for StORF_Num, data in current_Contig_URs.items():
        if is_break == True:
            break
        if data != None:
            #for UR_StORF, UR_StORF_data in UR_StORFs.items():
            ur_pos = data[7]
            StORF_Pos_In_UR = data[0]
            StORF_Start_In_UR = int(StORF_Pos_In_UR.split(',')[0])
            StORF_Stop_In_UR = int(StORF_Pos_In_UR.split(',')[1])
            ur_frame = data[2]
            strand = data[3]
            key_start = int(ur_pos.split('_')[0])
            StORF_start = key_start + StORF_Start_In_UR
            StORF_stop = key_start + StORF_Stop_In_UR
            allow_start = StORF_start + StORF_options.allowed_overlap
            StORF_Seq = data[1]
            StORF_Length = data[4]
            StORF_UR_Num = data[6]
        if allow_start >= track_current_stop:
            gff_start = str(StORF_Start_In_UR + int(ur_pos.split('_')[-1].split('_')[0]))
            gff_stop = str(StORF_Stop_In_UR + int(ur_pos.split('_')[-1].split('_')[0]))
            if strand == '+':  # Get original genome frame
                frame = (int(gff_start) % 3) + 1
            elif strand == '-':
                frame = (int(gff_stop) % 3) + 4
            StORFs.append([ur_pos,StORF_start, StORF_stop, StORF_Start_In_UR, StORF_Stop_In_UR, frame, ur_frame, strand,
                           StORF_Length, StORF_UR_Num, StORF_Seq])

    return StORFs




def StORF_Filler(StORF_options,Reported_StORFs):
    gff_in = open(StORF_options.gff, 'r')
    outfile = open(StORF_options.gff.replace('.gff','') + '_StORF_Reporter_Combined.gff', 'w')
    track_prev_contig, track_contig = '',''
    StORF_Num = 0
    end = False
    written_line = None
    Contig_URS = collections.defaultdict()
    for Contig, URs in Reported_StORFs.items():
        UR_StORF_Num = 0
        Contig_StORFs = collections.OrderedDict()
        for UR in URs:
            for StORF, data in UR.items():
                data.insert(0, StORF)
                Contig_StORFs.update({UR_StORF_Num: data})
                UR_StORF_Num += 1

        Contig_URS[Contig] = Contig_StORFs

 #   if StORF_options.prokka_dir == True:
    ### Rename .gff to .fasta / faa and load in fasta file
    if StORF_options.prokka_dir == True: # If normal PROKKA Dir run
        fasta_outfile = StORF_options.gff.replace('.gff','')
        fasta_outfile = open(fasta_outfile + '_StORF_Reporter_Combined.fasta','w')
        faa_infile = StORF_options.gff.replace('.gff', '.faa')
        ffn_infile = StORF_options.gff.replace('.gff', '.ffn')
        PROKKA_CoDing,PROKKA_Non_CoDing = FASTA_Load(faa_infile,ffn_infile)
    ##this will be a dict of ids and seqs
    for line in gff_in:
        if not line.startswith('#') and end == False:
            data = line.split('\t')
            track_current_start = int(data[3])
            track_current_stop = int(data[4])
            track_contig = data[0]
            if track_contig != track_prev_contig:  # End of current contig
                if track_prev_contig != '': # Get last StORF on Contig - If Present
                    StORFs = find_after_StORFs(StORF_options, Contig_URS, track_prev_start, track_prev_stop,track_prev_contig) # Changed to prev stop because we are switching from previous contig
                    if StORFs:
                        for StORF in StORFs:
                            GFF_StoRF_write(StORF_options, track_prev_contig, outfile, StORF,
                                            StORF_Num)  # To keep consistency
                            if StORF_options.prokka_dir == True:
                                FASTA_StoRF_write(track_contig, fasta_outfile, StORF)
                            StORF_Num += 1
                track_prev_start, track_prev_stop = 0, 0
            track_prev_contig = track_contig
            if track_current_start == track_prev_start and track_current_stop == track_prev_stop:  # `duplicate' entry in GFF
                if StORF_options.verbose == True:
                    print("skip")
            else:
                StORFs, Contig_URS = find_prev_StORFs(StORF_options,Contig_URS, track_current_start, track_prev_stop, track_contig)
            track_prev_start = track_current_start
            track_prev_stop = track_current_stop
            ##### Print out PROKKA Protein
            if 'gene' in data[2] and StORF_options.prokka_dir == True:
                PROKKA_ID = data[8].split(';')[0].replace('ID=','').replace('_gene','')
                try:
                    PROKKA_Seq = PROKKA_CoDing[PROKKA_ID]
                except KeyError:
                    try:
                        PROKKA_Seq = PROKKA_Non_CoDing[PROKKA_ID]
                    except KeyError:
                        print("PROKKA seq not found")
                fasta_outfile.write('>'+PROKKA_ID+'\n'+PROKKA_Seq+'\n')
            if StORFs:
                for StORF in StORFs:  # ([ur_pos,StORF_start, StORF_stop, StORF_Start_In_UR, StORF_Stop_In_UR, frame, ur_frame, strand, StORF_Length, StORF_UR_Num, StORF_Seq)]
                    GFF_StoRF_write(StORF_options, track_contig, outfile, StORF, StORF_Num)  # To keep consistency
                    if StORF_options.prokka_dir == True:
                        FASTA_StoRF_write(track_contig, fasta_outfile, StORF)
                    StORF_Num += 1
            if line != written_line:
                outfile.write(line)
                written_line = line
            StORFs = None
        elif line.startswith('##FASTA'):
            end = True
            StORFs = find_after_StORFs(StORF_options, Contig_URS, track_prev_start, track_prev_stop, track_prev_contig)  # Changed to prev stop because we are switching from previous contig
            if StORFs:
                for StORF in StORFs:
                    GFF_StoRF_write(StORF_options, track_prev_contig, outfile, StORF, StORF_Num)  # To keep consistency
                    if StORF_options.prokka_dir == True:
                        FASTA_StoRF_write(track_contig, fasta_outfile, StORF)
                    StORF_Num += 1
            if line != written_line:
                outfile.write(line)
                written_line = line

        else:
            if line != written_line:
                outfile.write(line)
                written_line = line
##############################################################


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF-Reporter v0.5.0: StORF_Reporter Run Parameters.')
    parser.add_argument('-anno', action='store', dest='annotation_type', default='PROKKA', const='PROKKA', required=True,
                        choices=['PROKKA', 'Ensembl', 'CDS'], nargs='?',
                        help='Default - PROKKA: Annotation type to be StORF-Reported:'
                             'Options: PROKKA = "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR"'
                             ';Ensembl = "ID=gene" ;CDS = "CDS"')
    parser.add_argument('-pd', '--PROKKA_dir', action='store', dest='prokka_dir', default='', required=False,
                        help='PROKKA output directory to be used if PROKKA chosen - Produces a new GFF and FASTA containing '
                             'all Coding and Non-Coding Seqs')
    parser.add_argument('-p_gff', '--PROKKA_GFFs', action='store', dest='prokka_gffs', default='', required=False,
                        help='Provide directory contain GFFs to be StORFed - Only produces modified GFFs')
    parser.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                        help='Default - 30: Minimum UR Length')
    parser.add_argument('-max_len', action='store', dest='maxlen', default='100000', type=int,
                        help='Default - 100,000: Maximum UR Length')
    parser.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                        help='Default - 50: UR Extension Length')
    parser.add_argument('-type', action='store', dest='feature_type', default='CDS', const='CDS', nargs='?',
                        choices=['StORF', 'CDS', 'ORF'],
                        help='Default - "CDS": Which GFF feature type for StORFs to be reported as in GFF - '
                             '"CDS" is probably needed for use in tools such as Roary')
    parser.add_argument('-olap', action="store", dest='overlap_nt', default=50, type=int,
                        help='Default - 50: Maximum number of nt of a StORF which can overlap another StORF.')
    parser.add_argument('-ao', action="store", dest='allowed_overlap', default=50, type=int,
                        help='Default - 50 nt: Maximum overlap between a StORF and an original gene.')


    parser.add_argument('-lw', action="store", dest='line_wrap', default=True, type=eval, choices=[True, False],
                        help='Default - True: Line wrap FASTA sequence output at 60 chars')
    parser.add_argument('-aa', action="store", dest='translate', default=False, type=eval, choices=[True, False],
                        help='Default - False: Report StORFs as amino acid sequences')

    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-v', action='store', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')

    Reporter_options = parser.parse_args()

    ######### Setup for PROKKA
    if Reporter_options.annotation_type == 'PROKKA' and Reporter_options.prokka_dir != "":
        Reporter_options.gene_ident = "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR"
        Reporter_options.nout = True
        Contigs, Reporter_options = run_UR_Extractor_PROKKA_DIR(Reporter_options)
        ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
        StORF_options = Namespace(reporter=True, gff=Reporter_options.gff, stop_codons="TGA,TAA,TAG",
                                  partial_storf=False, whole_contig=False, storf_order='start_pos',
                                  con_storfs=False, con_only=False, max_orf=50000, filtering='hard',
                                  feature_type=Reporter_options.feature_type, overlap_nt=Reporter_options.overlap_nt,
                                  allowed_overlap=Reporter_options.allowed_overlap, prokka_dir=True, prokka_gffs='',
                                  minlen=30, maxlen=100000, min_orf=100, verbose=False, nout=True)

        Reporter_StORFs = StORF_Reported(Contigs, StORF_options)
        StORF_Filler(StORF_options,Reporter_StORFs)  # Append StORFs to current GFF file from provided genome annotation
        print("Finished: " + Reporter_options.gff)
        ####
    elif Reporter_options.annotation_type == 'PROKKA' and Reporter_options.prokka_gffs != "":
        Reporter_options.gene_ident = "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR"
        Reporter_options.nout = True

        gff_list = list(pathlib.Path(Reporter_options.prokka_gffs).glob('*.gff'))
        for gff in gff_list:
            Contigs, Reporter_options = run_UR_Extractor_PROKKA_GFFs(Reporter_options,gff)
            ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
            StORF_options = Namespace(reporter=True, gff=Reporter_options.gff, stop_codons="TGA,TAA,TAG",
                                      partial_storf=False, whole_contig=False,
                                      con_storfs=False, con_only=False, max_orf=50000, filtering='hard',
                                      feature_type=Reporter_options.feature_type, storf_order='start_pos',
                                      overlap_nt=Reporter_options.overlap_nt, prokka_dir='', prokka_gffs=True,
                                      allowed_overlap=Reporter_options.allowed_overlap,
                                      minlen=30, maxlen=100000, min_orf=100, verbose=False, nout=True)
            Reporter_StORFs = StORF_Reported(Contigs, StORF_options)
            #sequence_id, all_StORFs = StORF_Reported(URs, StORF_options)
            StORF_Filler(StORF_options, Reporter_StORFs)  # Append StORFs to current GFF file from provided genome annotation
            print("Finished: " + Reporter_options.gff) # Will add number of additional StORFs here
    else:
        print("Please provide correct PROKKA output directory")









