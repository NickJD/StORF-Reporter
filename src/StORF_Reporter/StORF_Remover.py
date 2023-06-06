import argparse
import collections
from datetime import date
import gzip
import sys


try:
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from Constants import *

def gff_load_and_write(options,gff_in,blast_hits):
    for line in gff_in: 
        line_data = line.split()
        if line.startswith('\n') or line.startswith('#'):  # Not to crash on empty lines in GFF
            options.gff_out.write(line)
        else:
            if 'StORF-Reporter' in line:
                StORF = line_data[8].split('ID=')[1].split(';')[0]
                if StORF in blast_hits:
                    options.gff_out.write(line)
            else:
                options.gff_out.write(line)


def load_blast_6(options, blast_in):
    blast_hits = collections.OrderedDict()
    for line in blast_in:
        line_data = line.split('\t')
        if 'StORF' in line_data[0]:
            if int(line_data[11].split('.')[0]) >= options.minscore:
                blast_hits.update({line_data[0]:[]})

    return blast_hits

def remover(options):
    try:
        try:
            blast_in = gzip.open(options.blast,'rt')
            blast_hits = load_blast_6(options,blast_in)
        except:
            blast_in = open(options.blast,'r')#,encoding='unicode_escape')
            blast_hits = load_blast_6(options,blast_in)
        try:
            gff_in = gzip.open(options.gff,'rt')
            gff_load_and_write(options,gff_in,blast_hits)
        except:
            gff_in = open(options.gff,'r')#,encoding='unicode_escape')
            gff_load_and_write(options,gff_in,blast_hits)
    except AttributeError:
        sys.exit("Attribute Error")


def main():

    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': UR-Extractor Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=False)
    required.add_argument('-blast', action='store', dest='blast', help='BLAST format 6 annotation file',
                        required=False)

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-min_score', action='store', dest='minscore', default='30', type=int,
                        help='Minimum BitScore to keep StORF: Default 30')

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


    options = parser.parse_args()
    if options.gff == None or options.blast == None:
        if options.version:
            sys.exit(StORF_Reporter_Version)
        else:
            exit('StORF-Remover: error: the following arguments are required: -f, -gff, -blast')

    print("Thank you for using StORF-Reporter -- A detailed user manual can be found at https://github.com/NickJD/StORF-Reporter\n"
          "Please report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")

    options.annotation_type = [None,None]

    #### Output Directory and Filename handling
    if options.o_dir == None and options.o_name == None:
        tmp_extension = options.gff.split('.')[-1]
        output_file = options.gff.replace('.' + tmp_extension, '')
        output_file = output_file + '_NoHit_StORFs_Removed'
    elif options.o_dir != None and options.o_name != None:
        output_file = options.o_dir
        output_file = output_file + '/' if not output_file.endswith('/') else output_file
        output_file = output_file + options.o_name
    elif options.o_dir != None:
        tmp_extension = options.gff.split('.')[-1]
        output_file = options.gff.replace('.' + tmp_extension, '').split('/')[-1]
        output_file = options.o_dir + output_file + '_NoHit_StORFs_Removed'
    elif options.o_name != None:
        tmp_filename = options.gff.split('/')[-1]
        output_file = options.gff.replace(tmp_filename, '')
        output_file = output_file + options.o_name

    if options.gz == False:
        options.gff_out =  open(output_file + '.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        options.gff_out =  gzip.open(output_file + '.gff.gz','wt', newline='\n', encoding='utf-8')

    remover(options)


if __name__ == "__main__":
    main()
    print("Complete")



