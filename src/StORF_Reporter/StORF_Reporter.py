import argparse
import pathlib
import collections
import textwrap
import hashlib
import pyrodigal
from tempfile import NamedTemporaryFile
import os
import sys
import gzip


try:
    from .UR_Extractor import extractor
    from .StORF_Finder import StORF_Reported
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from UR_Extractor import extractor
    from StORF_Finder import StORF_Reported
    from Constants import *


class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return argparse.HelpFormatter._split_lines(self, text, width)

# class StORF_obj: # TODO
#     def __init__(self, name, roll):
#         self.name = name
#         self.roll = roll

################### We are currently fixed using Table 11
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

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate
############################
def get_directory_names(path):
    directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    return directories
############################

def compute_hash(StORF, Reporter_options,track_contig):
    ### Compute the hash/locus tag here
    ID = track_contig + '_UR_' + StORF[0] + '_' + StORF[10] + '_' + str(StORF[9])
    try:
        to_hash = Reporter_options.gff.split('.')[0] + '_' + ID  # create unique hash from outfile name and ID
    except AttributeError: # Will go here for Pyrodigal
        to_hash = Reporter_options.fasta.split('.')[
                      0] + '_' + ID  # create unique hash from outfile name and ID
    if len(StORF[11].split(',')) >= 3:  # fix for con-storfs - TBD
        StORF_Type = 'Con-StORF_'
    else:
        StORF_Type = 'StORF_'
    StORF_Hash = StORF_Type + hashlib.shake_256(to_hash.encode()).hexdigest(8)
    return StORF_Hash, ID


############################
def get_outfile_name(Reporter_options):
    if Reporter_options.o_dir == None and Reporter_options.o_name == None and Reporter_options.alt_filename == None:
        if Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs']:
            output_file = Reporter_options.path
            split_path = output_file.split(os.sep)
            directory_name = split_path[-1]
            output_file = os.path.join(output_file, f"{directory_name}_StORF-Reporter_Extended")
        else:
            if os.path.isdir(Reporter_options.path):
                output_file = Reporter_options.path
            else:
                tmp_extension = Reporter_options.path.split('.')[-1]  # could be .fa/.fasta etc
                output_file = Reporter_options.path.replace('.' + tmp_extension, '')
            if 'Multiple' in Reporter_options.annotation_type[1]:
                if Reporter_options.annotation_type[0] == 'Pyrodigal':
                    output_file = os.path.join(output_file, "_Pyrodigal_StORF-Reporter_Extended")
                else:
                    output_file = os.path.join(output_file, "_StORF-Reporter_Extended")
            else:
                if Reporter_options.annotation_type[0] == 'Pyrodigal':
                    output_file = output_file + '_Pyrodigal_StORF-Reporter_Extended'
                else:
                    output_file = output_file + '_StORF-Reporter_Extended'

    elif Reporter_options.o_dir == None and Reporter_options.o_name == None and Reporter_options.alt_filename != None:
        output_file = Reporter_options.path
        output_file = os.path.join(output_file, f"{Reporter_options.alt_filename}_StORF-Reporter_Extended")

    elif Reporter_options.o_dir != None and Reporter_options.o_name != None:
        Reporter_options.o_dir = os.path.normpath(Reporter_options.o_dir)
        Reporter_options.o_dir = os.path.realpath(Reporter_options.o_dir)
        output_file = Reporter_options.o_dir
        output_file = os.path.join(output_file, Reporter_options.o_name)
    elif Reporter_options.o_dir != None:
        Reporter_options.o_dir = os.path.normpath(Reporter_options.o_dir)
        Reporter_options.o_dir = os.path.realpath(Reporter_options.o_dir)
        if Reporter_options.alt_filename != None:
            directory_name = Reporter_options.alt_filename
        else:
            split_path = Reporter_options.path.split(os.sep)
            directory_name = split_path[-1]
        if Reporter_options.annotation_type[1] == 'Out_Dir':
            output_file = os.path.join(Reporter_options.o_dir, f"{directory_name}_StORF-Reporter_Extended")
        elif Reporter_options.annotation_type[0] == 'Pyrodigal':
            split_path = Reporter_options.path.split(os.sep)
            filename = split_path[-1].split('.')[0]
            output_file = os.path.join(Reporter_options.o_dir, f"{filename}_Pyrodigal_StORF-Reporter_Extended")
        else:
            split_path = Reporter_options.path.split(os.sep)
            filename = split_path[-1].split('.')[0]
            output_file = os.path.join(Reporter_options.o_dir, f"{filename}_StORF-Reporter_Extended")

    elif Reporter_options.o_name != None:
        if 'Multiple' in Reporter_options.annotation_type[1]:
            output_file = os.path.join(Reporter_options.path, Reporter_options.o_name)
        elif Reporter_options.annotation_type[1] == 'Out_Dir':
            output_file = os.path.join(Reporter_options.path, f"{Reporter_options.o_name}")
        else:
            split_path = Reporter_options.path.split(os.sep)
            filename = split_path[-1]
            output_file = Reporter_options.path.replace(filename,Reporter_options.o_name)

    return output_file

############################

def GFF_StORF_write(Reporter_options, track_contig, gff_out, StORF, StORF_Num, StORF_Hash, ID): # Consistency in outfile
    ### Write out new GFF entry -
    strand = StORF[7]
    start = StORF[3]
    stop = StORF[4]
    UR_Start = int(StORF[0].split('_')[0])
    start_stop = StORF[12][0:3]
    if len(StORF[11].split(',')) >= 3: # fix for con-storfs - TBD
        StORF_Type = 'ID=Con-StORF_'
        if StORF[7] == '-':
            mid = int(StORF[11].split(',')[2]) - int(StORF[11].split(',')[1])
            mid_stop = StORF[12][mid:mid + 3]
        else:
            mid = int(StORF[11].split(',')[1]) - int(StORF[11].split(',')[0])
            mid_stop = StORF[12][mid:mid + 3]
    else:
        mid_stop = 'N/A'
        StORF_Type = 'ID=StORF_'
    end_stop = StORF[12][-3:]

    if strand == '+':
        gff_start = max(start + 1 + UR_Start, 1) # this may need to change to account for different exlen
        gff_stop = max(stop + UR_Start, 1)
        if UR_Start == 1:
            gff_start -= 1
            gff_stop -= 1

        if Reporter_options.stop_inclusive == False:  # To remove the start and stop codon positions.
            gff_start = gff_start + 3
        frame = (int(gff_stop) % 3) + 1
    elif strand == '-':
        gff_start = max(start - 2 + UR_Start, 1)
        gff_stop = max(stop - 3 + UR_Start, 1)
        if Reporter_options.stop_inclusive == False:  # To remove the start and stop codon positions.
            gff_stop = gff_stop - 3
        frame = (int(gff_stop) % 3) + 4

    StORF_length = int(gff_stop) + 1 - int(gff_start) # +1 to adjust for base-1

    gff_out.write(track_contig + '\tStORF-Reporter\t' + Reporter_options.feature_type + '\t' +  str(gff_start) + '\t' + str(gff_stop) + '\t.\t' +
        StORF[7] + '\t0\tID=' + StORF_Hash + ';locus_tag=' + ID + ';INFO=Additional_Annotation_StORF-Reporter;UR_Stop_Locations=' + StORF[11].replace(',','-') + ';Name=' +
           StORF[10] + '_' + str(StORF_Num) + ';' + StORF[10] + '_Num_In_UR=' + str(StORF[9]) + ';' + StORF[10] + '_Length=' + str(StORF_length) + ';' + StORF[10] +
           '_Frame=' + str(frame) + ';UR_' + StORF[10] + '_Frame=' + str(StORF[6]) +  ';Start_Stop=' + start_stop + ';Mid_Stops=' + mid_stop  + ';End_Stop='
           + end_stop + ';StORF_Type=' + StORF[10] + '\n')


def FASTA_StORF_write(Reporter_options, fasta_out, StORF, StORF_Hash):  # Consistency in outfile

    ### Wrtie out new FASTA entry - Currently only write out as nt
    fasta_out.write('>'+StORF_Hash+'\n')
    sequence = StORF[-1]
    if Reporter_options.translate == True:
        sequence = translate_frame(sequence[0:])
        if Reporter_options.remove_stop == True:
            sequence = sequence.strip('*')
    if Reporter_options.line_wrap == True:
        wrapped = textwrap.wrap(sequence, width=60)
        for wrap in wrapped:
            fasta_out.write(wrap + '\n')
    else:
        fasta_out.write(sequence+'\n')
    if Reporter_options.storfs_out == True and Reporter_options.annotation_type[1] == 'Out_Dir': # Not efficient
        if Reporter_options.gz == False:
            storf_fasta_outfile = open(fasta_out.name.replace('.fasta','_StORFs_Only.fasta'),'a')
        else:
            storf_fasta_outfile = open(fasta_out.name.replace('.fasta.gz','_StORFs_Only.fasta.gz'),'a')
        storf_fasta_outfile.write('>' + StORF_Hash + '\n')
        if Reporter_options.line_wrap == True:
            wrapped = textwrap.wrap(sequence, width=60)
            for wrap in wrapped:
                storf_fasta_outfile.write(wrap+'\n')
        else:
            storf_fasta_outfile.write(sequence+'\n')

def FASTA_Load(faa_infile,ffn_infile):
    ### Coding sequences first
    Prokka_AA = collections.defaultdict()
    first = True
    infile = open(faa_infile)
    seq = ''
    ##load in fasta file and make dict of seqs and ids
    for line in infile:
        line = line.strip()
        if line.startswith('>') and first == False:
            Prokka_AA.update({id:seq})
            id = line.replace('>','')
            id = id.split(' ')[0]
            seq = ''
        elif line.startswith('>'):
            id = line.replace('>','')
            id = id.split(' ')[0]
            first = False
        else:
            seq += line
    Prokka_AA.update({id:seq})
    ### Coding and Non-Coding sequences in NT form next
    Prokka_NT = collections.defaultdict()
    first = True
    infile = open(ffn_infile)
    seq = ''
    ##load in fasta file and make dict of seqs and ids
    for line in infile:
        line = line.strip()
        if line.startswith('>') and first == False:
            Prokka_NT.update({id:seq})
            id = line.replace('>','')
            id = id.split(' ')[0]
            seq = ''
        elif line.startswith('>'):
            id = line.replace('>','')
            id = id.split(' ')[0]
            first = False
        else:
            seq += line
    Prokka_NT.update({id:seq})
    return Prokka_AA,Prokka_NT

def read_pyrodigal_fasta(fasta_in,sequences):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: seq})
            seq = ''
            sequence_name = line.split(' ')[0].replace('>', '')
        elif line.startswith('>'):
            seq = ''
            sequence_name = line.split(' ')[0].replace('>', '')
        else:
            seq += str(line)
            first = False
        sequences.update({sequence_name: seq})
    return sequences

def pyrodigal_predict(fasta,Reporter_options):
    sequences = collections.OrderedDict()
    try:  # Detect whether fasta/gff files are .gz or text and read accordingly
        fasta_in = gzip.open(fasta, 'rt')
        sequences = read_pyrodigal_fasta(fasta_in, sequences)
    except:
        fasta_in = open(fasta, 'r', encoding='unicode_escape')  # Not sure if needed in long term
        sequences = read_pyrodigal_fasta(fasta_in, sequences)
    pyrodigal_hold = NamedTemporaryFile(mode='w+', delete=False)

    orf_finder = pyrodigal.GeneFinder()
    longest = (max(sequences.values(), key=len)) # Get longest contig for training

    if Reporter_options.py_train == 'meta':
        orf_finder = pyrodigal.GeneFinder(meta=True)
    elif Reporter_options.py_train == 'longest':
        orf_finder.train(longest)
    for sequence_name, seq in sequences.items():
        if Reporter_options.py_train == 'indvidual':
            try:  # Catch if sequence is less than 20kb - will run in meta mode
                orf_finder.train(seq)
                break
            except ValueError:
                orf_finder = pyrodigal.GeneFinder(meta=True)  # If NO contig is >= 20kb#
        else:
            genes = orf_finder.find_genes(seq)
        if Reporter_options.py_unstorfed == True:
            genes.write_genes(Reporter_options.pyrodigal_out_file_fasta, sequence_id=sequence_name)
            genes.write_gff(Reporter_options.pyrodigal_out_file_gff, sequence_id=sequence_name)
        genes.write_gff(pyrodigal_hold, sequence_id=sequence_name)
    if Reporter_options.py_unstorfed == True:
        Reporter_options.pyrodigal_out_file_fasta.close()
        Reporter_options.pyrodigal_out_file_gff.close()
    pyrodigal_hold.close()
    Reporter_options.gff = pyrodigal_hold
    Reporter_options.fasta = fasta
    URs = extractor(Reporter_options)
    return URs, Reporter_options


def run_UR_Extractor_Directory(Reporter_options): # When given a complete Prokka/Bakta Directory
    if Reporter_options.alt_filename:
        identifier = Reporter_options.alt_filename
    else:
        split_path = Reporter_options.path.split(os.sep)
        identifier = split_path[-1]
    Reporter_options.fasta = os.path.join(Reporter_options.path, f"{identifier}.fna")
    if Reporter_options.annotation_type[0] == 'Prokka':
        Reporter_options.gff = Reporter_options.fasta.replace('.fna','.gff')
    if Reporter_options.annotation_type[0] == 'Bakta':
        Reporter_options.gff = Reporter_options.fasta.replace('.fna', '.gff3')
    URs = extractor(Reporter_options)
    return URs,Reporter_options


def run_UR_Extractor_Extended_GFFs(Reporter_options,gff): # When given a directory with multiple GFFs but without accompianing .fna
    if '_StORF-Reporter_Extended' not in str(gff): #Might fall over - put a break
        gff = str(gff)
        Reporter_options.gff = gff
    Reporter_options.fasta = gff
    URs = extractor(Reporter_options)
    return URs,Reporter_options

def run_UR_Extractor_Matched(Reporter_options,gff): # When given a directory with multiple GFFs but with accompianing .fna
    if '_StORF-Reporter_Extended' not in str(gff): #Might fall over - put a break
        gff = str(gff)
        Reporter_options.gff = gff
        fasta = gff.replace('.gff','.fasta')
        fna = gff.replace('.gff', '.fna')
        fa = gff.replace('.gff', '.fa')
        if os.path.isfile(fasta):
            Reporter_options.fasta = fasta
        elif os.path.isfile(fna):
            Reporter_options.fasta = fna
        elif os.path.isfile(fa):
            Reporter_options.fasta = fa
        else:
            sys.exit('No matching FASTA file found for ' + gff)
    URs = extractor(Reporter_options)
    return URs,Reporter_options

def run_UR_Extractor_GFF(Reporter_options, gff): # When given a directory with multiple GFFs but without accompianing .fna
    if '_StORF-Reporter_Extended' not in gff: #Might fall over - put a break
        gff = gff
        Reporter_options.gff = gff
    fasta = str(gff)
    Reporter_options.fasta = fasta
    URs = extractor(Reporter_options)
    return URs,Reporter_options

def find_prev_StORFs(StORF_options, Reporter_options, Contig_URs, track_current_start, track_prev_stop, track_contig):
    StORFs_to_del = []
    StORFs = []
    is_break = False
    try:
        current_Contig_URs = Contig_URs[track_contig]
        for StORF_Num, data in current_Contig_URs.items():
            if is_break == True:
                break
            if data != None:
                ur_pos = data[8] # extended UR
                ur_start = int(ur_pos.split('_')[0])
                StORF_Stop_Pos = data[0]
                StORF_Start_In_UR = int(StORF_Stop_Pos.split(',')[0])
                StORF_Stop_In_UR = int(StORF_Stop_Pos.split(',')[-1])
                ur_frame = data[2]
                strand = data[3]
                gff_start = str(StORF_Start_In_UR + 1 + max(ur_start - Reporter_options.exlen, 0))
                gff_stop = str(StORF_Stop_In_UR + max(ur_start - Reporter_options.exlen, 0))
                allow_start = int(gff_start) + StORF_options.allowed_overlap
                allow_stop = int(gff_stop) - StORF_options.allowed_overlap
                StORF_Seq = data[1]
                StORF_Length = data[4]
                StORF_UR_Num = data[6]
                if allow_stop <= track_current_start:# and allow_start >= track_prev_stop:
                    if strand == '+': # Get original genome frame
                        frame = (int(gff_start) % 3) + 1
                    elif strand == '-':
                        frame = (int(gff_stop) % 3) + 4
                    StORFs.append([ur_pos,gff_start,gff_stop,StORF_Start_In_UR,StORF_Stop_In_UR,frame,ur_frame,strand,
                                   StORF_Length,StORF_UR_Num,data[5],StORF_Stop_Pos,StORF_Seq])
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

    except KeyError:
        if StORF_options.verbose == True:
            print("GFF formatting error - Something wrong with seqeunce region " + track_contig)
        return StORFs, Contig_URs


def find_after_StORFs(StORF_options,Contig_URs,track_current_start,track_current_stop, track_contig):
    StORFs = []
    is_break = False
    try:
        current_Contig_URs = Contig_URs[track_contig]
        for StORF_Num, data in current_Contig_URs.items():
            if is_break == True:
                break
            if data != None:
                # for UR_StORF, UR_StORF_data in UR_StORFs.items():
                ur_pos = data[8] # extended UR
                StORF_Stop_Pos = data[0]
                StORF_Start_In_UR = int(StORF_Stop_Pos.split(',')[0])
                StORF_Stop_In_UR = int(StORF_Stop_Pos.split(',')[-1])
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
                StORFs.append(
                    [ur_pos, StORF_start, StORF_stop, StORF_Start_In_UR, StORF_Stop_In_UR, frame, ur_frame, strand,
                     StORF_Length, StORF_UR_Num, data[5], StORF_Stop_Pos, StORF_Seq])

        return StORFs
    except KeyError:
        if StORF_options.verbose == True:
            print("GFF formatting error - Something wrong with sequence region " + track_contig)
        return StORFs

def StORF_Filler(Reporter_options, Reported_StORFs):
    if Reporter_options.annotation_type[0] == 'Pyrodigal':
        with open(Reporter_options.gff.name, "r") as gff:
            gff_in = gff.read()
            gff_name = Reporter_options.path.split('.')[0]
    else:
        with open(Reporter_options.gff, 'r') as gff:
            gff_in = gff.read()
            gff_name = Reporter_options.gff
    try:
        if Reporter_options.gz == False:
            Reporter_options.gff_outfile = open(Reporter_options.output_file + '.gff', 'w', newline='\n', encoding='utf-8')
        elif Reporter_options.gz == True:
            Reporter_options.gff_outfile = gzip.open(Reporter_options.output_file + '.gff.gz', 'wt', newline='\n', encoding='utf-8')
    except FileNotFoundError:
        sys.exit('Out file location error - Try providing full path.')

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
    ### Rename .gff to .fasta / faa and load in fasta file
    if Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs']: # or Reporter_options.annotation_type[0] == 'Bakta': # If normal Prokka/Bakta Dir run
        faa_infile = Reporter_options.gff.replace('.gff3', '.faa').replace('.gff','.faa')
        ffn_infile = Reporter_options.gff.replace('.gff3', '.ffn').replace('.gff','.ffn')
        Original_AA,Original_NT = FASTA_Load(faa_infile,ffn_infile)
    if Reporter_options.annotation_type[0] == 'Pyrodigal':
        split_path = Reporter_options.fasta.split(os.sep)
        file_name = split_path[-1]
        Reporter_options.gff_outfile.write('##Pyrodigal annotation and StORF-Reporter extended GFF annotation of ' + file_name + '\n') # Windows uses '\'
    else:
        split_path = Reporter_options.gff.split(os.sep)
        file_name = split_path[-1]
        Reporter_options.gff_outfile.write('##StORF-Reporter extended annotation of ' + file_name + '\n')

    #Handling separate output of FASTA sequences
    if Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs'] or Reporter_options.storfs_out == True:
        if Reporter_options.translate == True and Reporter_options.gz == False:
            fasta_outfile = open(Reporter_options.output_file + '_aa.fasta', 'w', newline='\n', encoding='utf-8')
        elif Reporter_options.translate == False and Reporter_options.gz == False:
            fasta_outfile = open(Reporter_options.output_file + '.fasta', 'w', newline='\n', encoding='utf-8')
        elif Reporter_options.translate == True and Reporter_options.gz == True:
            fasta_outfile = gzip.open(Reporter_options.output_file + '_aa.fasta.gz', 'wt', newline='\n', encoding='utf-8')
        elif Reporter_options.translate == False and Reporter_options.gz == True:
            fasta_outfile = gzip.open(Reporter_options.output_file + '.fasta.gz', 'wt', newline='\n', encoding='utf-8')

    Reporter_options.gff_outfile.write('##StORF-Reporter ' + StORF_Reporter_Version + '\n')

    first_region = True
    for line in gff_in.splitlines( ):
        if not line.startswith('#') and end == False:
            data = line.split('\t')
            track_current_start = int(data[3])
            track_current_stop = int(data[4])
            track_contig = data[0]
            if data[2] == 'region':
                Reporter_options.gff_outfile.write(line.strip() + '\n')
                first_region = False
                continue
            if track_contig != track_prev_contig:  # End of current contig
                # if track_prev_contig != '': # Get last StORF on Contig - If Present
                #     StORFs = find_after_StORFs(StORF_options, Contig_URS, track_prev_start, track_prev_stop,track_prev_contig) # Changed to prev stop because we are switching from previous contig
                #     if StORFs:
                #         for StORF in StORFs:
                #             GFF_StORF_write( Reporter_options, track_prev_contig, outfile, StORF,
                #                             StORF_Num)  # To keep consistency
                #             if StORF_options.path == True:
                #                 FASTA_StORF_write(Reporter_options, track_contig, fasta_outfile, StORF)
                #             StORF_Num += 1
                track_prev_start, track_prev_stop = 0, 0
            track_prev_contig = track_contig
            if track_current_start == track_prev_start and track_current_stop == track_prev_stop:  # `duplicate' entry in GFF
                tracked = True
            else:
                StORFs, Contig_URS = find_prev_StORFs(Reporter_options, Reporter_options, Contig_URS, track_current_start, track_prev_stop, track_contig)
                tracked = False
            track_prev_start = track_current_start
            track_prev_stop = track_current_stop
            ##### Print out Prokka/Bakta Gene
            if tracked == False and (Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs'] or Reporter_options.storfs_out == True):
                if ('gene' in data[2] or 'ID=' in data[8]) and Reporter_options.annotation_type[0] in ('Prokka', 'Bakta') and Reporter_options.annotation_type[1] == 'Out_Dir':
                    Original_ID = data[8].split(';')[0].replace('ID=','').replace('_gene','')
                    try:
                        Original_Seq = Original_NT[Original_ID] # used to be AA but changed for to NT. CHECK
                    except (KeyError,UnboundLocalError):
                        try:
                            Original_Seq = Original_NT[Original_ID]
                        except (KeyError,UnboundLocalError):
                            #if Reporter_options.verbose == True:
                            sys.exit("Error: Original seq " + Original_ID + " not found")
                    fasta_outfile.write('>'+Original_ID+'\n')
                    if Reporter_options.translate == True:
                        Original_Seq = translate_frame(Original_Seq[0:])
                        if Reporter_options.remove_stop == True:
                            Original_Seq = Original_Seq.strip('*')
                    if Reporter_options.line_wrap == True:
                        wrapped = textwrap.wrap(Original_Seq, width=60)
                        for wrap in wrapped:
                            fasta_outfile.write(wrap + '\n')
                    else:
                        fasta_outfile.write(Original_Seq+'\n')
            if StORFs:
                for StORF in StORFs:
                    ###Compute hash/locus tag
                    StORF_Hash, ID = compute_hash(StORF,Reporter_options, track_contig)
                    GFF_StORF_write(Reporter_options, track_contig, Reporter_options.gff_outfile, StORF, StORF_Num, StORF_Hash, ID)  # To keep consistency
                    if (Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs']  or Reporter_options.storfs_out == True):
                        FASTA_StORF_write(Reporter_options, fasta_outfile, StORF, StORF_Hash)
                    StORF_Num += 1
            if line != written_line:
                Reporter_options.gff_outfile.write(line.strip()+'\n')
                written_line = line
            StORFs = None
        elif line.startswith('##sequence-region') and first_region != True:
            StORFs = find_after_StORFs(Reporter_options, Contig_URS, track_prev_start, track_prev_stop, track_prev_contig)  # Changed to prev stop because we are switching from previous contig
            if StORFs:
                for StORF in StORFs:
                    ###Compute hash/locus tag
                    StORF_Hash, ID = compute_hash(StORF,Reporter_options, track_contig)
                    GFF_StORF_write(Reporter_options, track_prev_contig, Reporter_options.gff_outfile, StORF, StORF_Num, StORF_Hash, ID)  # To keep consistency
                    if Reporter_options.annotation_type[1] in ['Out_Dir', 'Multiple_Out_Dirs'] or Reporter_options.storfs_out == True:
                        FASTA_StORF_write(Reporter_options, fasta_outfile, StORF, StORF_Hash)
                    StORF_Num += 1
            Reporter_options.gff_outfile.write(line.strip() + '\n')

        elif line.startswith('##FASTA'):
            Reporter_options.gff_outfile.write(line.strip() + '\n')
            end = True
            # StORFs = find_after_StORFs(StORF_options, Contig_URS, track_prev_start, track_prev_stop, track_prev_contig)  # Changed to prev stop because we are switching from previous contig
            # if StORFs:
            #     for StORF in StORFs:
            #         GFF_StORF_write(Reporter_options, track_prev_contig, outfile, StORF)  # To keep consistency
            #         if StORF_options.path == True:
            #             FASTA_StORF_write(Reporter_options, track_contig, fasta_outfile, StORF)
            #         StORF_Num += 1
            # if line != written_line:
            #     outfile.write(line.strip()+'\n')
            #     written_line = line

        else: # Now we just print out the remaining lines in the GFF
            Reporter_options.gff_outfile.write(line.strip()+'\n')
            written_line = line
##############################################################



def main():
    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': StORF-Reporter Run Parameters.', formatter_class=SmartFormatter)
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Options')
    required.add_argument('-anno', action='store', dest='annotation_type', required=False,
                        choices=['Prokka', 'Bakta', 'Out_Dir', 'Multiple_Out_Dirs', 'Single_GFF', 'Multiple_GFFs', 'Ensembl', 'Feature_Types',
                                  'Single_Genome', 'Multiple_Genomes','Single_Combined_GFF', 'Multiple_Combined_GFFs',
                                 'Pyrodigal', 'Single_FASTA', 'Multiple_FASTA'], nargs='*',
                        help='R|Select Annotation and Input options for one of the 3 options listed below\n'
                             '### Prokka/Bakta Annotation Option 1: \n'
                             '\tProkka = Report StORFs for a Prokka annotation; \n'
                             '\tBakta = Report StORFs for a Bakta annotation; \n'
                             '--- Prokka/Bakta Input Options: \n'
                             '\tOut_Dir = To provide the output directory of either a Prokka or Bakta run (will produce a new GFF and FASTA file containing original and extended annotations); \n'
                             '\tMultiple_Out_Dirs = To provide a directory containing multiple Prokka/Bakta standard output directories - Will run on each sequentially; \n' 
                             '\tSingle_GFF = To provide a single Prokka or Bakta GFF - searches for accompanying ".fna" file (will provide a new extended GFF); \n'
                             '\tMultiple_GFFs = To provide a directory containing multiple Prokka or Bakta GFF files - searches for accompanying ".fna" files (will provide a new extended GFF); \n\n'
                             
                             '### Standard GFF Annotation Option 2: \n'
                             '\tEnsembl = Report StORFs for an Ensembl Bacteria annotation (ID=gene); \n'
                             '\tFeature_Types = Used in conjunction with -gene_ident to define features such as CDS,rRNA,tRNA for UR extraction (default CDS); \n'
                             '--- Standard GFF Input Options: \n'
                             '\tSingle_Genome = To provide a single Genome - accompanying FASTA must share same name as given gff file (can be .fna, .fa or .fasta); \n'
                             '\tMultiple_Genomes = To provide a directory containing multiple accompanying GFF and FASTA files - files must share the same name (fasta can be .fna, .fa or .fasta); \n'
                             '\tSingle_Combined_GFF = To provide a GFF file with embedded FASTA at the bottom; \n'
                             '\tMultiple_Combined_GFFs = To provide a directory containing multiple GFF files with embedded FASTA at the bottom; \n\n'
                             
                             '### Complete Annotation Option 3: \n'
                             '\tPyrodigal = Run Pyrodigal then Report StORFs (provide path to single FASTA or directory of multiple FASTA files ;\n'
                             '--- Complete Annotation Input Options: \n'
                             '\tSingle_FASTA = To provide a single FASTA file; \n'
                             '\tMultiple_FASTA = To provide a directory containing multiple FASTA files (will detect .fna,.fa,.fasta); \n\n')

    required.add_argument('-p', action='store', dest='path', default='', required=False,
                        help='Provide input file or directory path')
    ###
    StORF_Reporter_args = parser.add_argument_group('StORF-Reporter Options')
    ### Add options to redirect into new output directory and filename - and output StORFs on their own in their own GFF -
    StORF_Reporter_args.add_argument('-af', action="store", dest='alt_filename', required=False,
                        help='Default - Prokka/Bakta output directory share the same prefix with their gff/fna files - Use this option when Prokka/Bakta output directory name is different '
                             'from the gff/fna files within and StORF-Reporter will search for the gff/fna with the given prefix (MyProkkaDir/"altname".gff) - Does not work with "Multiple_Out_Dirs" option')
    StORF_Reporter_args.add_argument('-oname', action="store", dest='o_name', required=False,
                        help='Default - Appends \'_StORF-Reporter_Extended\' to end of input filename - Takes the directory name of Prokka/Bakta output if given as input or the input for -af if given'
                             '  - Multiple_* runs will be numbered')
    StORF_Reporter_args.add_argument('-odir', action="store", dest='o_dir', required=False,
                        help='Default -  Same directory as input')
    StORF_Reporter_args.add_argument('-sout', action="store", dest='storfs_out', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out StORF sequences separately from Prokka/Bakta annotations')
    StORF_Reporter_args.add_argument('-lw', action="store", dest='line_wrap', default=True, type=eval, choices=[True, False],
                        help='Default - True: Line wrap FASTA sequence output at 60 chars')
    StORF_Reporter_args.add_argument('-aa', action="store", dest='translate', default=False, type=eval, choices=[True, False],
                        help='Default - False: Report StORFs as amino acid sequences')
    # output.add_argument('-gff_fasta', action="store", dest='gff_fasta', default=False, type=eval, choices=[True, False],
    #                     help='Default - False: Report all gene sequences (nt) at the bottom of GFF files in Prokka output mode')
    StORF_Reporter_args.add_argument('-gz', action='store', dest='gz', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    ###
    prodigal_options = parser.add_argument_group('Pyrodigal Options')
    prodigal_options.add_argument('-py_train', action='store', dest='py_train', required=False,
                        choices=['longest', 'individual', 'meta'], nargs='?', default='longest',
                        help='Default - longest: Type of model training to be done for Pyrodigal CDS prediction: '
                             'Options: longest = Trains on longest contig;'
                             ' individual = Trains on each contig separately - runs in meta mode if contig is < 20KB;'
                             ' meta = Runs in meta mode for all sequences')
    prodigal_options.add_argument('-py_fasta', action="store", dest='py_fasta', default=True, type=eval, choices=[True, False],
                                  help='Default - False: Output Pyrodigal+StORF predictions in FASTA format')
    prodigal_options.add_argument('-py_unstorfed', action="store", dest='py_unstorfed', default=False, type=eval, choices=[True, False],
                                  help='Default - False: Provide GFF containing original Pyrodigal predictions')

    ###
    UR_Ex_args = parser.add_argument_group('UR-Extractor Options')
    UR_Ex_args.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                          help='Identifier used for extraction of Unannotated Regions such as "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"'
                               ' - To be used with "-anno Feature_Types" - "-gene_ident Prokka" will select features present in Prokka annotations')
    UR_Ex_args.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                          help='Default - 30: Minimum UR Length')
    UR_Ex_args.add_argument('-max_len', action='store', dest='maxlen', default='100000', type=int,
                          help='Default - 100,000: Maximum UR Length')
    UR_Ex_args.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                          help='Default - 50: UR Extension Length')
    ###
    StORF_Finder_args = parser.add_argument_group('StORF-Finder Options')
    StORF_Finder_args.add_argument('-spos', action="store", dest='stop_inclusive', default=False, type=eval,
                          choices=[True, False],
                          help='Default - False: Output StORF positions inclusive of first stop codon')
    StORF_Finder_args.add_argument('-rs', action="store", dest='remove_stop', default=True, type=eval, choices=[True, False],
                        help='Default - True: Remove stop "*" from StORF amino acid sequences')

    StORF_Finder_args.add_argument('-con_storfs', action="store", dest='con_storfs', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output Consecutive StORFs')
    StORF_Finder_args.add_argument('-con_only', action="store", dest='con_only', default=False, type=eval, choices=[True, False],
                        help='Default - False: Only output Consecutive StORFs')
    StORF_Finder_args.add_argument('-ps', action="store", dest='partial_storf', default=False, type=eval, choices=[True, False],
                        help='Default - False: Partial StORFs reported')
    StORF_Finder_args.add_argument('-wc', action="store", dest='whole_contig', default=False, type=eval, choices=[True, False],
                        help='Default - False: StORFs reported across entire sequence')
    StORF_Finder_args.add_argument('-short_storfs', action="store", dest='short_storfs', default=False, type=str, choices=['False', 'Nolap', 'Olap'],
                        help='Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt '
                             'that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs which are'
                             'overlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs. Overlap is defined by "-olap".')
    StORF_Finder_args.add_argument('-short_storfs_only', action="store", dest='short_storfs_only', default=False, type=eval, choices=[True, False],
                        help='Default - True. Only report Short-StORFs?')
    StORF_Finder_args.add_argument('-minorf', action="store", dest='min_orf', default=99, type=int,
                          help='Default - 99: Minimum StORF size in nt')
    StORF_Finder_args.add_argument('-maxorf', action="store", dest='max_orf', default=60000, type=int,
                        help='Default - 60kb: Maximum StORF size in nt')
    StORF_Finder_args.add_argument('-codons', action="store", dest='stop_codons', default="TAG,TGA,TAA",
                          help='Default - (\'TAG,TGA,TAA\'): List Stop Codons to use')
    StORF_Finder_args.add_argument('-olap_filt', action='store', dest='olap_filtering', default='both-strand',
                          const='both-strand', nargs='?',
                          choices=['none', 'single-strand', 'both-strand'],
                          help='Default - "both-strand": Filtering level "none" is not recommended, "single-strand" for single strand filtering '
                               'and both-strand for both-strand longest-first tiling')
    StORF_Finder_args.add_argument('-start_filt', action="store", dest='start_filtering', default=False, type=eval,
                          choices=[True, False],
                          help='Default - False: Filter out StORFs without at least one of the 3 common start codons (best used for short-storfs).')
    StORF_Finder_args.add_argument('-so', action="store", dest='storf_order', default='start_pos', nargs='?', choices=['start_pos','strand'],
                        required=False,
                        help='Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.')
    StORF_Finder_args.add_argument('-f_type', action='store', dest='feature_type', default='CDS', const='CDS', nargs='?',
                          choices=['StORF', 'CDS', 'ORF'],
                          help='Default - "CDS": Which GFF feature type for StORFs to be reported as in GFF - '
                               '"CDS" is probably needed for use in tools such as Roary and Panaroo')
    StORF_Finder_args.add_argument('-olap', action="store", dest='overlap_nt', default=50, type=int,
                          help='Default - 50: Maximum number of nt of a StORF which can overlap another StORF.')
    StORF_Finder_args.add_argument('-ao', action="store", dest='allowed_overlap', default=50, type=int,
                          help='Default - 50 nt: Maximum overlap between a StORF and an original gene.')
    # optional.add_argument('-col', action='store', dest='genome_collection', default='', required=False,
    #                     help='Provide directories containing a collection of matching FASTA and GFF files'
    #                          '(list .fa then .gff containing directories separated by commas - ./FA,./GFF)')
    #optional.add_argument('-comb', action='store', dest='combined_gffs', default='', required=False,
    #                    help='Provide directory containing GFFs with sequences combined into single file to be StORFed - Only produces modified GFFs')

    misc = parser.add_argument_group('Misc')


    misc.add_argument('-overwrite', action='store', dest='overwrite', default=False, type=eval, choices=[True, False],
                        help='Default - False: Overwrite StORF-Reporter output if already present')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages')
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Print out version number and exit')
    misc.add_argument('-nout', action='store', dest='nout', default=True, type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)
    misc.add_argument('-nout_pyrodigal', action='store', dest='nout_pyrodigal', default=True, type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)

    Reporter_options = parser.parse_args()
    Reporter_options.gff_fasta = None # To be done
    Reporter_options.reporter = True

    if Reporter_options.annotation_type == None or Reporter_options.path == None:
        if Reporter_options.version:
            sys.exit(StORF_Reporter_Version)
        else:
            exit('StORF-Reporter: error: the following arguments are required: -anno, -p')

    if Reporter_options.translate == True and Reporter_options.storfs_out == False and Reporter_options.annotation_type[1] not in ['Out_Dir', 'Multiple_Out_Dirs']:
        exit('StORF-Reporter: "-sout True" is required when "-aa True" is selected')


    print("Thank you for using StORF-Reporter -- A detailed user manual can be found at https://github.com/NickJD/StORF-Reporter\n"
          "Please report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")

    ##############
    ## TODO: Print out user chosen annotation options

    # # Get the used options provided by the user
    # used_options = {k: v for k, v in vars(Reporter_options).items() if v is not None}
    # # Print the used options
    # print("Options used by the user:")
    # for key, value in used_options.items():
    #     print(f"{key}: {value}")

    ##############
    Reporter_options.path = os.path.normpath(Reporter_options.path)
    Reporter_options.path = os.path.realpath(Reporter_options.path)

    ##############
    ## Incompatible argument catching
    if len(Reporter_options.annotation_type) != 2:
        parser.error('Please select two compatible options for required argument -anno.')
    ##############


    #############
    ## User selects Prokka 'standard' Features
    if Reporter_options.gene_ident == 'Prokka':
        Reporter_options.gene_ident = 'misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo'
    #############

    if Reporter_options.annotation_type[1] != 'Multiple_Out_Dirs':
        output_file = get_outfile_name(Reporter_options)



# ========================================================================================================================

    ############## Setup for Prokka/Bakta output directory
    if Reporter_options.annotation_type[0] in ('Prokka','Bakta') and Reporter_options.annotation_type[1] == 'Out_Dir':
        Reporter_options.output_file = output_file
        #### Checking and cleaning
        try:
            for fname in os.listdir(Reporter_options.path):
                if '_StORF-Reporter_Extended' in fname and Reporter_options.overwrite == False:
                    parser.error(
                        'Prokka/Bakta directory not clean and already contains a StORF-Reporter output. Please delete or use "-overwrite True" and try again.')
                elif '_StORF-Reporter_Extended' in fname and Reporter_options.overwrite == True:
                    file_path = os.path.join(Reporter_options.path, fname)
                    os.remove(file_path)
                    if Reporter_options.verbose == True:
                        print('StORF-Reporter output ' + fname + ' will be overwritten.')
        except FileNotFoundError:
            sys.exit("Incorrect file path '" + Reporter_options.path + "' - Please check input")
        ####
        Reporter_options.gene_ident = "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"
        Contigs, Reporter_options = run_UR_Extractor_Directory(Reporter_options)
        ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
        Reporter_StORFs = StORF_Reported(Reporter_options, Contigs)
        StORF_Filler(Reporter_options, Reporter_StORFs)
        if Reporter_options.verbose == True:
            print("Finished: " + Reporter_options.gff.split(os.sep)[-1])

    ############## Setup for Multi[;e Prokka/Bakta output directories
    elif Reporter_options.annotation_type[0] in ('Prokka','Bakta') and Reporter_options.annotation_type[1] == 'Multiple_Out_Dirs':
        fixed_path = Reporter_options.path # So we can modify the path variable later
        directories = get_directory_names(fixed_path)
        for directory in directories:
            directory_path = os.path.join(fixed_path, directory)
            #### Checking and cleaning
            files_in_directory = os.listdir(directory_path)
            if len(files_in_directory) == 0:
                if Reporter_options.verbose == True:
                    split_path = directory_path.split(os.sep)
                    print("Directory '" + split_path[-1] + "' missing fna/gff files")
                continue
            #if fname.endswith('.gff') or file_name.endswith('.fna'):
            #
            for fname in files_in_directory:
                if '_StORF-Reporter_Extended' in fname and Reporter_options.overwrite == False:
                    parser.error(
                        'Prokka/Bakta directory not clean and already contains a StORF-Reporter output. Please delete or use "-overwrite True" and try again.')
                elif '_StORF-Reporter_Extended' in fname and Reporter_options.overwrite == True:
                    file_path = os.path.join(directory_path, fname)
                    os.remove(file_path)
                    if Reporter_options.verbose == True:
                        print('StORF-Reporter output ' + fname + ' will be overwritten.')
            ####
            Reporter_options.path = directory_path
            output_file = get_outfile_name(Reporter_options)
            Reporter_options.output_file = output_file
            Reporter_options.gene_ident = "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"
            Contigs, Reporter_options = run_UR_Extractor_Directory(Reporter_options)
            ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
            Reporter_StORFs = StORF_Reported(Reporter_options, Contigs)
            StORF_Filler(Reporter_options, Reporter_StORFs)
            if Reporter_options.verbose == True:
                print("Finished: " + Reporter_options.gff.split(os.sep)[-1])

    ############### Setup for directory of a single Prokka/Bakta GFF or a directory of multiple Prokka/Bakta GFFs
    elif Reporter_options.annotation_type[0] in ('Prokka', 'Bakta') and Reporter_options.annotation_type[1] in ('Single_GFF', 'Multiple_GFFs'):
        if Reporter_options.annotation_type[1] == "Single_GFF":
            gff_list = [Reporter_options.path]
        elif Reporter_options.annotation_type[1] == "Multiple_GFFs":
            gff_list = list(pathlib.Path(Reporter_options.path).glob('*.gff'))
            gff_list.extend(pathlib.Path(Reporter_options.path).glob('*.gff3'))
        #### Checking and cleaning
        gff_list = list(map(str, gff_list))
        gffs_to_filter = []
        for gff in gff_list:
            if '_StORF-Reporter_Extended.gff' in gff and os.path.isfile(gff) and Reporter_options.overwrite == False:
                gffs_to_filter.append(gff)
                gffs_to_filter.append(gff.replace('_StORF-Reporter_Extended.gff','.gff'))
                print('Prokka/Bakta GFF has already been processed and a StORF-Reporter output exists for ' + gff.split(os.sep)[-1] + '. Please delete or use "-overwrite True" and try again.')
            elif '_StORF-Reporter_Extended.gff' in gff and os.path.isfile(gff) and Reporter_options.overwrite == True:
                os.remove(gff)
                gffs_to_filter.append(gff)
                if Reporter_options.verbose == True:
                    print('StORF-Reporter output '  + gff.split(os.sep)[-1] +  ' will be overwritten.')
        gff_list = [x for x in gff_list if x not in gffs_to_filter]
        ####
        Reporter_options.gene_ident = "misc_RNA,gene,mRNA,CDS,rRNA,tRNA,tmRNA,CRISPR,ncRNA,regulatory_region,oriC,pseudo"
        file_counter = 0
        for gff in gff_list:
            # Finalising output_file name
            if Reporter_options.annotation_type[1]  == "Multiple_GFFs" and Reporter_options.o_name != None:
                Reporter_options.output_file = output_file + '_' + str(file_counter)
                file_counter += 1
            elif Reporter_options.annotation_type[1]  == "Multiple_GFFs":
                tmp_filename = gff.split(os.sep)[-1].split('.gff')[0]  # could be .fa/.fasta etc
                Reporter_options.output_file = output_file.replace('_StORF-Reporter',tmp_filename + '_StORF-Reporter')
            else:
                Reporter_options.output_file = output_file
            if Reporter_options.verbose == True:
                print("Starting: " + str(gff.split(os.sep)[-1]))
            Contigs, Reporter_options = run_UR_Extractor_GFF(Reporter_options,gff) # used to be extended
            ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
            Reporter_StORFs = StORF_Reported(Reporter_options, Contigs)
            StORF_Filler(Reporter_options, Reporter_StORFs)
            if Reporter_options.verbose == True:
                print("Finished: " + gff.split(os.sep)[-1]) # Will add number of additional StORFs here

    ################ Run StORF-Reporter on standardGFFs
    elif Reporter_options.annotation_type[0] in ('Ensembl', 'Feature_Types'):
        if Reporter_options.annotation_type[0] == 'Ensembl':
            Reporter_options.gene_ident = "ID=gene"
        elif Reporter_options.annotation_type[0] == 'Feature_Types':
            Reporter_options.gene_ident = Reporter_options.gene_ident.split(',')
        ####
        if Reporter_options.annotation_type[1] in ("Single_Genome", "Single_Combined_GFF"):
            gff_list = [Reporter_options.path]
        elif Reporter_options.annotation_type[1] in ("Multiple_Genomes", "Multiple_Combined_GFFs"):
            gff_list = list(pathlib.Path(Reporter_options.path).glob('*.gff'))
            gff_list.extend(pathlib.Path(Reporter_options.path).glob('*.gff3'))
        else:
            parser.error('Please select two compatible options for required argument -anno.')
        #### Checking and cleaning
        gff_list = list(map(str, gff_list))
        gffs_to_filter = []
        for gff in gff_list:
            if '_StORF-Reporter_Extended.gff' in gff and os.path.isfile(gff) and Reporter_options.overwrite == False:
                gffs_to_filter.append(gff)
                gffs_to_filter.append(gff.replace('_StORF-Reporter_Extended.gff', '.gff'))
                print('Prokka/Bakta GFF has already been processed and a StORF-Reporter output exists for ' +
                      gff.split(os.sep)[-1] + '. Please delete or use "-overwrite True" and try again.')
            elif '_StORF-Reporter_Extended.gff' in gff and os.path.isfile(gff) and Reporter_options.overwrite == True:
                os.remove(gff)
                gffs_to_filter.append(gff)
                if Reporter_options.verbose == True:
                    print('StORF-Reporter output ' + gff.split(os.sep)[-1] + ' will be overwritten.')
        gff_list = [x for x in gff_list if x not in gffs_to_filter]
        ####
        file_counter = 0
        for gff in gff_list:
            # Finalising output_file name
            if (Reporter_options.annotation_type[1]  == "Multiple_Combined_GFFs" or Reporter_options.annotation_type[1] == "Multiple_Genomes") and Reporter_options.o_name != None:
                Reporter_options.output_file = output_file + '_' + str(file_counter)
                file_counter += 1
            elif (Reporter_options.annotation_type[1]  == "Multiple_Combined_GFFs" or Reporter_options.annotation_type[1] == "Multiple_Genomes"):
                tmp_filename = gff.split(os.sep)[-1].split('.gff')[0]  # could be .fa/.fasta etc
                Reporter_options.output_file = output_file.replace('_StORF-Reporter',tmp_filename + '_StORF-Reporter')
            else:
                Reporter_options.output_file = output_file

            if Reporter_options.verbose == True:
                print("Starting: " + str(gff.split(os.sep)[-1]))
            if Reporter_options.annotation_type[1] in ('Single_Combined_GFF', 'Multiple_Combined_GFFs'):
                Contigs, Reporter_options = run_UR_Extractor_Extended_GFFs(Reporter_options, gff)
            else:
                Contigs, Reporter_options = run_UR_Extractor_Matched(Reporter_options, gff)
            ################## Find StORFs in URs - Setup StORF_Reporter-Finder Run
            Reporter_StORFs = StORF_Reported(Reporter_options, Contigs)
            StORF_Filler(Reporter_options, Reporter_StORFs)
            if Reporter_options.verbose == True:
                print("Finished: " + gff.split(os.sep)[-1])  # Will add number of additional StORFs here`

    ####### Run Pyrodigal and then StORF-Reporter on either a single or directory of FASTA files
    elif Reporter_options.annotation_type[0] == 'Pyrodigal' and Reporter_options.annotation_type[1] in ("Single_FASTA", "Multiple_FASTA"):  # needs cleaning
        Reporter_options.gene_ident = "gene"
        if os.path.isdir(Reporter_options.path) and Reporter_options.annotation_type[1]  == "Multiple_FASTA":
            fasta_list = list(pathlib.Path(Reporter_options.path).glob('*.fasta'))
            fasta_list.extend(list(pathlib.Path(Reporter_options.path).glob('*.fna')))
            fasta_list.extend(list(pathlib.Path(Reporter_options.path).glob('*.fa')))
            fasta_list = list(map(str, fasta_list))
        elif os.path.isfile(Reporter_options.path) and Reporter_options.annotation_type[1]  == "Single_FASTA":
            fasta_list = [Reporter_options.path]
        else:
            parser.error('Please provide a directory or file path to a genomic FASTA file.')
        file_counter = 0
        for fasta in fasta_list:
            # Finalising output_file name
            if Reporter_options.annotation_type[1]  == "Multiple_FASTA" and Reporter_options.o_name != None:
                Reporter_options.output_file = output_file + '_' + str(file_counter)
                file_counter += 1
            elif Reporter_options.annotation_type[1]  == "Multiple_FASTA":
                file_name = fasta.split(os.sep)[-1]
                Reporter_options.output_file = output_file.replace('_Pyrodigal',file_name.split('.')[0] + '_Pyrodigal')
            else:
                Reporter_options.output_file = output_file
            if Reporter_options.py_unstorfed == True:  # User wants to output initial Pyrodigal GFF
                Reporter_options.pyrodigal_out_file_fasta = open(Reporter_options.output_file + '_Pyrodigal_Predictions.fasta','w')
                Reporter_options.pyrodigal_out_file_gff = open(Reporter_options.output_file + '_Pyrodigal_Predictions.gff', 'w')
            if Reporter_options.verbose == True:
                print("Starting: " + str(fasta.split(os.sep)[-1]))
            ###### Run Pyrodigal
            Contigs, Reporter_options = pyrodigal_predict(str(fasta),Reporter_options)
            ###### Find StORFs in URs - Setup StORF_Reporter-Finder Run
            Reporter_StORFs = StORF_Reported(Reporter_options, Contigs)
            ###### Append StORFs to current GFF file from provided genome annotation
            StORF_Filler(Reporter_options, Reporter_StORFs)
            split_path = fasta.split(os.sep)
            file_name = split_path[-1]
            print("Finished predicting Pyrodigal CDS genes and StORFs for -  " + file_name) # Can add number of additional StORFs here
    else:
        parser.error('Incompatible user-arguments selected.')




if __name__ == "__main__":
    main()
    print("Complete")


