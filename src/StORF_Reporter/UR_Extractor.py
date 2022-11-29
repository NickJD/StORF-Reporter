import argparse
import collections
from datetime import date
import gzip
import sys


try:
    # Calling from ORForise via pip
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *

#Output FASTA and GFF separately using the same out_filename but with respective extensions - gz output optional
def write_fasta(dna_regions, options):
    if options.out_file == None:
        options.out_file = options.gff.split('.gff')[0]
    if options.gz == False:
        out =  open(options.out_file + '_UR.fasta','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UR.fasta.gz', 'wt', newline='\n', encoding='utf-8')

    out.write("##\tUR Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    out.write("##Original Files: " + options.fasta + ' | ' + options.gff + '\n')
    for dna_region, dna_region_ur in dna_regions.items():
        out.write('\n##sequence-region\t' + dna_region + ' 1 ' + str(len(dna_region_ur[0])) + '\n')
        ur_ident = dna_region + options.ident # Add user ident onto name of dna regions
        if dna_region_ur[3]:
            for ex_ur, data in dna_region_ur[3].items():
                original_ur = data[0]
                ur_seq = data[1]
                out.write('>' + ur_ident + '_' + ex_ur + '\n' + ur_seq + '\n')
    out.close()

def write_gff(dna_regions,options):
    if options.out_file == None:
        options.out_file = options.gff.split('.gff')[0]
    if options.gz == False:
        out =  open(options.out_file + '_UR.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UR.gff.gz', 'wt', newline='\n', encoding='utf-8')
    out.write("##gff-version\t3\n#\tUR Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    out.write('##StORF-Reporter ' + StORF_Reporter_Version + '\n')
    for seq_reg in dna_regions:
        out.write('##sequence-region\t' + seq_reg + ' 1 ' + str(dna_regions[seq_reg][1]) + '\n')
    out.write("##Original Files: " + options.fasta + ' | ' + options.gff + '\n\n')
    for dna_region, dna_region_ur in dna_regions.items():
        ur_ident = dna_region + options.ident
        if dna_region_ur[3]:
            for ex_ur, data in dna_region_ur[3].items():
                length = len(data[1])
                ex_ur_pos = ex_ur.replace('_','\t')
                entry = (dna_region + '\tUR_Extractor\tunannotated_region\t' + ex_ur_pos + '\t.\t.\t.\tID='+ur_ident+'_'+ex_ur+';' + 'Original_UR=' + str(data[0]) + ';Note=UR_Length(Extended):' + str(length) + '\n')
                out.write(entry)
    out.close()

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
        #### Called when presented with PROKKA GFF file so must get fasta from inside it
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
    #Will code in different versions for different types of GFF3 files (Prodigal,Ensembl etc)
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n') or line.startswith('##'):  # Not to crash on empty lines in GFF
            continue
        elif options.gene_ident == 'ID=gene':
            if line_data[0] in dna_regions and options.gene_ident in line_data[8]:
                pos = line_data[3] + '_' + line_data[4]
                dna_regions[line_data[0]][2].append(pos) # This will add to list
        else:
            gene_types = options.gene_ident.split(',')
            try:
                if line_data[0] in dna_regions:
                    if any(gene_type in line_data[2] for gene_type in gene_types): # line[2] for normalrun
                        if options.verbose == True:
                            print(line_data[2])
                        pos = line_data[3] + '_' + line_data[4]
                        if pos not in dna_regions[line_data[0]][2]:
                            dna_regions[line_data[0]][2].append(pos) # This will add to list
            except IndexError:
                continue

    return dna_regions

def extractor(options):
    try:
        try: # Detect whether fasta/gff files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta,'rt')
            dna_regions = fasta_load(fasta_in)
        except:
            fasta_in = open(options.fasta,'r',encoding='unicode_escape') # Not sure if needed in long term
            dna_regions = fasta_load(fasta_in)
        try:
            gff_in = gzip.open(options.gff,'rt')
            dna_regions = gff_load(options,gff_in,dna_regions)
        except:
            gff_in = open(options.gff,'r',encoding='unicode_escape') # Not sure if needed in long term
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
                ###### This hack is to get over GFF errors where genome-long annotations
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
        write_fasta(dna_regions, options)
        write_gff(dna_regions, options)
    else:
        return dna_regions


def main():
    print("Thank you for using StORF-Reporter\nPlease report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")

    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': UR-Extractor Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-f', action='store', dest='fasta', required=True,
                        help='FASTA file for Unannotated Region seq extraction')
    required.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=True)

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
    output.add_argument('-o', action='store', dest='out_file', required=False,
                        help='Output file - Without filetype - default appends "_UR" to end of input gff filename (replaces \'.gff\')')
    output.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-nout', action='store', dest='nout', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-v', action='store', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')


    options = parser.parse_args()
    extractor(options)

    # Contig name could have a ';' which will mess up later on in StORF_Reporter-R
    # UR output should state original non extended

if __name__ == "__main__":
    main()
    print("Complete")
