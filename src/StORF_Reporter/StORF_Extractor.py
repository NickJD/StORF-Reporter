import argparse
import collections
from datetime import date
import gzip
import sys
import os
import pathlib


try:
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *

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

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate


def reverseCorrectLoci(seq_length,first,second,third): # here for the negative loci correction

    if second == None:
        corrected_start = max(seq_length - int(third),1)
        corrected_stop = max(seq_length - int(first-1),1)
        return corrected_start, corrected_stop
    else: # Needs to be checked
        corrected_start = max(seq_length - int(third),1)
        corrected_mid = max(seq_length - int(second-3),1)
        corrected_stop = max(seq_length - int(first),1)
        return corrected_start, corrected_mid, corrected_stop
#################################
def revCompIterative(watson): #Gets Reverse Complement

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev: # make dict to catch bad nts - if more than 1 output them to std error
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt
            #ns_nt[nt] +=1
    return crick


def write_fasta(dna_regions, fasta_outfile):
    for dna_region, dna_region_ur in dna_regions.items():
        fasta_outfile.write('##sequence-region\t' + dna_region + ' 1 ' + str(len(dna_region_ur[0])) + '\n')
        if dna_region_ur[3]:
            for storf, seq in dna_region_ur[3].items():
                fasta_outfile.write('>Start=' + storf.split('_')[0] + ';Stop=' + storf.split('_')[1] + ';Frame=' +
                                storf.split('_')[2] + storf.split(';')[1] + '\n' + seq + '\n')
    fasta_outfile.close()

def write_gff(dna_regions,options,gff_outfile, gff):
    gff_outfile.write("##gff-version\t3\n#\tStORF-Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    gff_outfile.write('##StORF-Reporter ' + StORF_Reporter_Version + '\n')
    for seq_reg in dna_regions:
        gff_outfile.write('##sequence-region\t' + seq_reg + ' 1 ' + str(dna_regions[seq_reg][1]) + '\n')
    gff_outfile.write("##Original File: " + gff.split('/')[-1] + '\n\n')
    for dna_region, storf_data in dna_regions.items():
        #ur_ident = dna_region + options.ident
        if storf_data[3]:
            for storf, seq in storf_data[3].items():

                entry = (dna_region + '\tStORF-Reporter\tCDS\t' + storf.split('_')[0] + '\t' + storf.split('_')[1]+ '\t.\t' +
                         storf.split('_')[2].split(';')[0]  + '\t.\t' + storf.split(';',1)[1] + '\n')
                gff_outfile.write(entry)
    gff_outfile.close()

def fasta_load(fasta_in):
    dna_regions = collections.OrderedDict()
    first = True
    if '>' in fasta_in.readline().rstrip():
        fasta_in.seek(0)
        #### Default for when presented with standard fasta file
        for line in fasta_in:
            line = line.strip()
            if line.startswith('>') and first == False:  # Check if first seq in file
                dna_region_length = len(seq)
                dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            elif line.startswith('>'):
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            else:
                seq += str(line)
                first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
    elif '##' in  fasta_in.readline().rstrip(): # Clunky and may fall over
        fasta_in.seek(0)
        #### Called when presented with Prokka GFF file so must get fasta from inside it
        ### Get to genome seq
        at_FASTA = False
        for line in fasta_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
            if line.startswith('##FASTA'):  # Not to crash on empty lines in GFF
                at_FASTA = True
            elif at_FASTA == True:
                line = line.strip()
                if line.startswith('>') and first == False:  # Check if first seq in file
                    dna_region_length = len(seq)
                    dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                elif line.startswith('>'):
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                else:
                    seq += str(line)
                    first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})

    return dna_regions


def gff_load(options,gff_in,dna_regions):
    has_storfs = False
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n') or line.startswith('#'):  # Not to crash on empty lines in GFF
            continue
        else:
            try:
                if line_data[0] in dna_regions:
                    if any(tool in line_data[1] for tool in options.tool_ident):
                        has_storfs = True
                        pos = line_data[3] + '_' + line_data[4] + '_' + line_data[6] + ';' + line_data[8]
                        if pos not in dna_regions[line_data[0]][2]:
                            dna_regions[line_data[0]][2].append(pos) # This will add to list
            except IndexError:
                continue
    if has_storfs == False:
        dna_regions = None
    return dna_regions


def run_StORF_Extractor_Combined_GFFs(options,gff): # When given a directory with multiple GFFs but without accompianing .fna
    gff = str(gff)
    options.gff = gff
    options.fasta = gff
    storf_extractor(options, gff)
    #return StORFs, options

def run_StORF_Extractor_Matched(options,gff): # When given a directory with multiple GFFs but with accompianing .fna
    gff = str(gff)
    options.gff = gff
    fasta = gff.replace('.gff','.fasta')
    fna = gff.replace('.gff', '.fna')
    fa = gff.replace('.gff', '.fa')
    if os.path.isfile(fasta):
        options.fasta = fasta
    elif os.path.isfile(fna):
        options.fasta = fna
    elif os.path.isfile(fa):
        options.fasta = fa
    else:
        sys.exit('No matching FASTA file found for ' + gff)
    storf_extractor(options, gff)
    #return StORFs, options


def storf_extractor(options, gff):
    try:
        try: # Detect whether fasta/gff files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta,'rt')
            dna_regions = fasta_load(fasta_in)
        except:
            fasta_in = open(options.fasta,'r',encoding='unicode_escape')
            dna_regions = fasta_load(fasta_in)
        try:
            gff_in = gzip.open(options.gff,'rt')
            dna_regions = gff_load(options,gff_in,dna_regions)
        except:
            gff_in = open(options.gff,'r',encoding='unicode_escape')
            dna_regions = gff_load(options,gff_in,dna_regions)
    except AttributeError:
        sys.exit("Attribute Error:\nStORF'ed GFF probably already exists - Must be deleted before running")
    if dna_regions == None and options.verbose == True:
        print("No StORFs to extract from " + gff)
        os.remove(options.fasta_outfile.name)
        os.remove(options.gff_outfile.name)
        return
    elif dna_regions == None:
        os.remove(options.fasta_outfile.name)
        if options.gff_out == True:
            os.remove(options.gff_outfile.name)
        return


    for (key,(seq,seq_length,posns,tmp))  in dna_regions.items(): #Extract URs from 1 dna_region at a time
        Extracted_StORFs = collections.OrderedDict()
        seq_rev = revCompIterative(seq)
        if posns: # If UR has a pos
            for pos in posns: # Iterate over GFF loci and measure flanking regions for potential URs
                start = int(pos.split('_')[0])
                stop = int(pos.split('_')[1])
                frame = pos.split('_')[2].split(';')[0]
                ###### This hack is to get over GFF errors where genome-long annotations
                if stop-start >= 100000:
                    if options.verbose == True:
                        print("UR " + pos + " is more than 100,000 kbs - Please Check Annotation")
                    continue
                if frame == '+':
                    StORF_seq = seq[start:stop]
                elif frame == '-':
                    rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(seq_length, start, None, stop)
                    StORF_seq = seq_rev[rev_corrected_start:rev_corrected_stop]

                Extracted_StORFs[pos] = StORF_seq
        dna_regions.update({key: (seq, seq_length, posns, Extracted_StORFs)})



    write_fasta(dna_regions, options.fasta_outfile)
    if options.gff_out != False:
        write_gff(dna_regions, options, options.gff_outfile, gff)

def main():

    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': StORF-Extractor Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-storf_input', action='store', dest='storf_input', required=False,
                        choices=['Combined', 'Separate'],
                        help='Are StORFs to be extracted from Combined GFF/FASTA or Separate GFF/FASTA files?\n')
    required.add_argument('-p', action='store', dest='path', default='', required=False,
                        help='Provide input file or directory path')

    ### Not implemented yet
    # optional = parser.add_argument_group('Optional Arguments')
    # optional.add_argument('-tool_ident', action='store', dest='tool_ident', default='StORF-Reporter',
    #                         help='Default "StORF-Reporter": Tool-name Identifier used for extraction of StORFs or other genomic elements'
    #                              ' "StORF-Reporter, Prodigal, Pyrodigal" (second GFF column)')


    output = parser.add_argument_group('Output')
    output.add_argument('-gff_out', action='store', dest='gff_out', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output StORFs in GFF format')
    output.add_argument('-oname', action="store", dest='o_name', required=False,
                        help='Default - Appends \'_Extracted_StORFs\' to end of input GFF filename')
    output.add_argument('-odir', action="store", dest='o_dir', required=False,
                        help='Default -  Same directory as input FASTA')
    output.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')


    misc = parser.add_argument_group('Misc')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages')
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Default - False: Print out version number and exit')



    options = parser.parse_args()
    options.tool_ident = ['StORF-Reporter']#options.tool_ident.split(',')


    if options.storf_input == None or options.path == None:
        if options.version:
            sys.exit(StORF_Reporter_Version)
        else:
            exit('StORF-Extractor: error: the following arguments are required: -storf_input, -p')

    print("Thank you for using StORF-Reporter -- A detailed user manual can be found at https://github.com/NickJD/StORF-Reporter\n"
          "Please report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")


    #### Output Directory and Filename handling
    if options.o_dir == None and options.o_name == None:
        tmp_extension = options.path.split('.')[-1]  # could be .fa/.fasta etc
        output_file = options.path.replace('.' + tmp_extension, '')
        output_file = output_file + '_Extracted_StORFs'
    elif options.o_dir != None and options.o_name != None:
        output_file = options.o_dir
        output_file = output_file + '/' if not output_file.endswith('/') else output_file
        output_file = output_file + options.o_name
    elif options.o_dir != None:
        tmp_extension = options.path.split('.')[-1]  # could be .fa/.fasta etc
        output_file = options.path.replace('.' + tmp_extension, '').split('/')[-1]
        output_file = options.o_dir + output_file + '_Extracted_StORFs'
    elif options.o_name != None:
        tmp_filename = options.path.split('/')[-1]  # could be .fa/.fasta etc
        output_file = options.path.replace(tmp_filename, '')
        output_file = output_file + options.o_name

######################################################################
    if os.path.isdir(options.path):
        options.path = options.path + '/' if not options.path.endswith('/') else options.path
        gff_list = list(pathlib.Path(options.path).glob('*.gff'))
        gff_list.extend(pathlib.Path(options.path).glob('*.gff3'))
    else:
        gff_list = [options.path]
    gff_list = list(map(str, gff_list))
    ####
    file_counter = 0
    for gff in gff_list:
        ### Checking if file has already been StORF-Extracted
        with open(gff) as f:
            second_line = f.readlines()[1]
            if 'StORF-Extractor' in second_line:
                if options.verbose == True:
                    print(gff + " is an extracted file and will be ignored.")
                continue
            #f.seek(0) # needed?
        # Finalising output_file name
        if os.path.isdir(options.path) and options.o_name != None:
            tmp_output_file = output_file + '_' + str(file_counter)
            file_counter += 1
        elif os.path.isdir(options.path):
            tmp_filename = gff.split('/')[-1].split('.gff')[0]  # could be .fa/.fasta etc
            tmp_output_file = output_file.replace('_Extracted_StORFs', tmp_filename + '_Extracted_StORFs')
        else:
            tmp_output_file = output_file

        if options.gz == False:
            options.fasta_outfile = open(tmp_output_file + '.fasta', 'w', newline='\n', encoding='utf-8')
            if options.gff_out == True:
                options.gff_outfile = open(tmp_output_file + '.gff', 'w', newline='\n', encoding='utf-8')
        elif options.gz == True:
            options.fasta_outfile = open(tmp_output_file + '.fasta.gz', 'wt', newline='\n', encoding='utf-8')
            if options.gff_out == True:
                options.gff_outfile = open(tmp_output_file + '.gff.gz', 'wt', newline='\n', encoding='utf-8')

        if options.verbose == True:
            print("Starting: " + str(gff))
        if  options.storf_input == 'Separate':
            run_StORF_Extractor_Matched(options, gff)
        else:
            run_StORF_Extractor_Combined_GFFs(options, gff)


        ##################

        if options.verbose == True:
            print("Finished: " + gff.split('/')[-1])  # Will add number of additional StORFs here`




if __name__ == "__main__":
    main()
    print("Complete")


























