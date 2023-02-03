import argparse
import collections
from datetime import date
import gzip
import sys


try:
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *

#Output FASTA and GFF separately using the same filename but with respective extensions - gz output optional
def write_fasta(dna_regions, options, fasta_out):
    fasta_out.write("##\tUR-Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    fasta_out.write('##StORF-Reporter ' + StORF_Reporter_Version + '\n')
    fasta_out.write("##Original Files: " + options.fasta.split('/')[-1] + ' | ' + options.gff.split('/')[-1] + '\n')
    for dna_region, dna_region_ur in dna_regions.items():
        fasta_out.write('\n##sequence-region\t' + dna_region + ' 1 ' + str(len(dna_region_ur[0])) + '\n')
        ur_ident = dna_region + options.ident # Add user ident onto name of dna regions
        if dna_region_ur[3]:
            for ex_ur, data in dna_region_ur[3].items():
                original_ur = data[0]
                ur_seq = data[1]
                fasta_out.write('>' + ur_ident + '_' + ex_ur + '\n' + ur_seq + '\n')
    fasta_out.close()

def write_gff(dna_regions,options,gff_out):
    gff_out.write("##gff-version\t3\n#\tUR-Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    gff_out.write('##StORF-Reporter ' + StORF_Reporter_Version + '\n')
    for seq_reg in dna_regions:
        gff_out.write('##sequence-region\t' + seq_reg + ' 1 ' + str(dna_regions[seq_reg][1]) + '\n')
    gff_out.write("##Original Files: " + options.fasta.split('/')[-1] + ' | ' + options.gff.split('/')[-1] + '\n\n')
    for dna_region, dna_region_ur in dna_regions.items():
        ur_ident = dna_region + options.ident
        if dna_region_ur[3]:
            for ex_ur, data in dna_region_ur[3].items():
                length = len(data[1])
                ex_ur_pos = ex_ur.replace('_','\t')
                entry = (dna_region + '\tUR_Extractor\tunannotated_region\t' + ex_ur_pos + '\t.\t.\t.\tID='+ur_ident+'_'+ex_ur+';' + 'Original_UR=' + str(data[0]) + ';Note=UR_Length(Extended):' + str(length) + '\n')
                gff_out.write(entry)
    gff_out.close()

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

def pyrodigal_virtual_gff_load(gff_in,dna_regions):
    with open(gff_in.name, "r") as gff:
        tmp = gff.read()
        for line in tmp.splitlines():
            line_data = line.split()
            if line.startswith('\n') or line.startswith('#'):  # Not to crash on empty lines in GFF
                continue
            else:
                try:
                    if line_data[0] in dna_regions:
                        if 'CDS' in line_data[2]: # line[2] for normal run
                            pos = line_data[3] + '_' + line_data[4]
                            if pos not in dna_regions[line_data[0]][2]:
                                dna_regions[line_data[0]][2].append(pos) # This will add to list
                except IndexError:
                    continue
    return dna_regions

def gff_load(options,gff_in,dna_regions):
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n') or line.startswith('#'):  # Not to crash on empty lines in GFF
            continue
        elif options.gene_ident[0] == 'ID=gene':
            if line_data[0] in dna_regions and options.gene_ident[0] in line_data[8]:
                pos = line_data[3] + '_' + line_data[4]
                dna_regions[line_data[0]][2].append(pos) # This will add to list
        else:
            try:
                if line_data[0] in dna_regions:
                    if any(gene_type in line_data[2] for gene_type in options.gene_ident): # line[2] for normal run
                        pos = line_data[3] + '_' + line_data[4]
                        if pos not in dna_regions[line_data[0]][2]:
                            dna_regions[line_data[0]][2].append(pos) # This will add to list
            except IndexError:
                continue
    return dna_regions

def extractor(options):
    if options.annotation_type[0] == 'Pyrodigal':
        try:  # Detect whether fasta/gff files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta, 'rt')
            dna_regions = fasta_load(fasta_in)
        except:
            fasta_in = open(options.fasta, 'r', encoding='unicode_escape')
            dna_regions = fasta_load(fasta_in)
        dna_regions = pyrodigal_virtual_gff_load(options.gff, dna_regions)
    else:
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


    for (key,(seq,seq_length,posns,URs))  in dna_regions.items(): #Extract URs from 1 dna_region at a time
        unannotated_regions = collections.OrderedDict()
        unannotated_start = 0
        if posns: # If UR has a pos
            for pos in posns: # Iterate over GFF loci and measure flanking regions for potential URs
                start = int(pos.split('_')[0])
                stop = int(pos.split('_')[1])
                ###### This hack is to account for GFF errors which contain genome-long annotations
                if stop-start >= 100000:
                    if options.verbose == True:
                        print("UR " + pos + " is more than 100,000 kbs - Please Check Annotation")
                    continue
                if start > unannotated_start:
                    length = start - unannotated_start
                    if length >= options.minlen and length <= options.maxlen: #default between 30 - 100,000
                        original_UR = str(max(unannotated_start,1)) + '_' + str(start)
                        unannotated_start = max(unannotated_start - options.exlen, 0) # 0 here for pythons base-o here is for the GFF base-1 reporting format
                        start += options.exlen
                        unannotated_seq = seq[unannotated_start:start]
                        unannotated_loci = str(max(unannotated_start,1)) + '_' + str(start)
                        unannotated_regions.update({unannotated_loci: [original_UR,unannotated_seq]})
                if stop > unannotated_start:
                    unannotated_start = stop
            try:
                if (seq_length - stop) >= options.minlen and (seq_length - stop) <= options.maxlen: #default between 30 - 100,000: # Get UR at end of dna_region if longer than minlen
                    original_UR = str(max(unannotated_start, 1)) + '_' + str(start)
                    stop -= options.exlen
                    unannotated_seq = seq[stop:seq_length]
                    unannotated_loci = str(stop) + '_' + str(seq_length)
                    unannotated_regions.update({unannotated_loci: [original_UR,unannotated_seq]})
            except UnboundLocalError:
                pass
            dna_regions.update({key: (seq, seq_length, posns, unannotated_regions)})


    if options.nout == False:
        write_fasta(dna_regions, options, options.fasta_out)
        write_gff(dna_regions, options, options.gff_out)
    else:
        return dna_regions


def main():

    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': UR-Extractor Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-f', action='store', dest='fasta', required=False,
                        help='FASTA file for Unannotated Region seq extraction')
    required.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=False)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-ident', action='store', dest='ident', default='_UR',
                        help='Identifier given for Unannotated Region output sequences - Do not modify if output is '
                             'to be used by StORF-Finder: Default "Sequence-ID"_UR')
    optional.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                        help='Minimum UR Length: Default 30')
    optional.add_argument('-max_len', action='store', dest='maxlen', default='100000', type=int,
                        help='Maximum UR Length: Default 100,000')
    optional.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                        help='UR Extension Length on 5\' and 3\': Default 50')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene',
                        help='Identifier used for extraction of Unannotated Regions "CDS,rRNA,tRNA": Default for Ensembl_Bacteria = '
                             '"ID=gene" or "-gene_ident CDS" for "most" genome annotations')

    output = parser.add_argument_group('Output')
    output.add_argument('-oname', action="store", dest='o_name', required=False,
                        help='Default - Appends \'_UR\' to end of input GFF filename')
    output.add_argument('-odir', action="store", dest='o_dir', required=False,
                        help='Default -  Same directory as input GFF')
    output.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')


    misc = parser.add_argument_group('Misc')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages')
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Default - False: Print out version number and exit')
    misc.add_argument('-nout', action='store', dest='nout', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)
    misc.add_argument('-pyrodigal', action='store', dest='pyrodigal', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)
    misc.add_argument('-nout_pyrodigal', action='store', dest='nout_pyrodigal', default='True', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)

    options = parser.parse_args()
    if options.fasta == None or options.gff == None:
        if options.version:
            sys.exit(StORF_Reporter_Version)
        else:
            exit('UR-Extractor: error: the following arguments are required: -f, -gff')

    print("Thank you for using StORF-Reporter -- A detailed user manual can be found at https://github.com/NickJD/StORF-Reporter\n"
          "Please report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")

    options.annotation_type = [None,None]

    #### Output Directory and Filename handling
    if options.o_dir == None and options.o_name == None:
        tmp_extension = options.gff.split('.')[-1]
        output_file = options.gff.replace('.' + tmp_extension, '')
        output_file = output_file + '_UR'
    elif options.o_dir != None and options.o_name != None:
        output_file = options.o_dir
        output_file = output_file + '/' if not output_file.endswith('/') else output_file
        output_file = output_file + options.o_name
    elif options.o_dir != None:
        tmp_extension = options.gff.split('.')[-1]
        output_file = options.gff.replace('.' + tmp_extension, '').split('/')[-1]
        output_file = options.o_dir + output_file + '_UR'
    elif options.o_name != None:
        tmp_filename = options.gff.split('/')[-1]
        output_file = options.gff.replace(tmp_filename, '')
        output_file = output_file + options.o_name

    if options.gz == False:
        options.fasta_out = open(output_file + '.fasta','w', newline='\n', encoding='utf-8')
        options.gff_out =  open(output_file + '.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        options.fasta_out = gzip.open(output_file + '.fasta.gz','wt', newline='\n', encoding='utf-8')
        options.gff_out =  gzip.open(output_file + '.gff.gz','wt', newline='\n', encoding='utf-8')

    options.gene_ident = options.gene_ident.split(',')
    extractor(options)


if __name__ == "__main__":
    main()
    print("Complete")




