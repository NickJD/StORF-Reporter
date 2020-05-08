import collections
import argparse
from datetime import date

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta',
                    help='FASTA file for Intergenic Region seq extraction')
parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA')
parser.add_argument('-ident', action='store', dest='ident', default='_IR',
                    help='Identifier given for Intergenic Region output sequences: Default "Input"_IR')
parser.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                    help='Minimum IR Length: Default 30')
parser.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                    help='IR Extension Length: Default 50')
parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene:',
                    help='Identifier used for extraction: Default = ID=gene:')
parser.add_argument('-o', '--output_file', action='store', dest='out_filename',
                    help='Output file name for Intergenic Regions to GFF and FASTA')

options = parser.parse_args()
fasta = options.fasta
gff = options.gff
ident = options.ident
minlen = options.minlen
exlen = options.exlen
gene_ident = options.gene_ident
out_file = options.out_filename

#Output FASTA and GFF separately using the same out_filename but with respective extensions
def write_fasta(dna_regions, ident, out_filename):
    with open(out_filename + '.fasta', 'w') as out_file:
        for dna_region, dna_region_ir in dna_regions.items():
            ident = dna_region + ident # Add user ident onto name of dna regions
            for ir, ir_seq in dna_region_ir[3].items():
                out_file.write('>' + ident + '|' + ir + '\n' + ir_seq + '\n')

def write_gff(dna_regions, ident, out_filename):
    with open(out_filename + '.gff', 'w') as out:
        out.write("##gff-version\t3\n#\tIR Extractor \n#\tRun Date:" + str(date.today()) + '\n')
        for dna_region, dna_region_ir in dna_regions.items():
            ident = dna_region + ident
            for ir, ir_seq in dna_region_ir[3].items():
                seq_id = ident + '|' + ir
                length = len(ir_seq)
                entry = (seq_id + '\tIR_Extractor\tDNA\t' + 'IR_Length(Extended):' + str(length) + '\n')
                out.write(entry)

def comparator(fasta, gff, ident, minlen, exlen, gene_ident, out_filename):
    dna_regions = collections.OrderedDict()
    first = True
    with open(fasta, 'r') as dna_seq: # Get dna_region sequences
        for line in dna_seq:
            line = line.strip()
            if ">" in line.strip() and first == False: # Check if first seq in file
                dna_region_length = len(seq)
                dna_regions.update({dna_region_id: (seq, str(dna_region_length), list())})
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            elif '>' in line:
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            else:
                seq += str(line)
                first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, str(dna_region_length), list())})

    with open(gff, 'r') as dna_gff: # Get gene locis from GFF - ID=Gene will also classify Pseudogenes as genes
        for line in dna_gff:
            line_data = line.split()
            if line_data[0] in dna_regions and gene_ident in line:
                pos = line_data[3] + '_' + line_data[4]
                dna_regions[line_data[0]][-1].append(pos)

    for key, value in dna_regions.items(): #Extract IRs from 1 dna_region at a time
        intergenic_regions = collections.OrderedDict()
        inter_start = 0
        for pos in value[2]: # Iterate over GFF loci and measure flanking regions for potential IRs
            start = int(pos.split('_')[0])
            stop = int(pos.split('_')[1])
            seq = value[0]
            if start > inter_start:
                length = start - inter_start
                if length >= minlen:
                    inter_start = max(inter_start - exlen, 0)
                    start += exlen
                    inter_seq = seq[inter_start:start]
                    inter_loci = str(inter_start) + '_' + str(start)
                    intergenic_regions.update({inter_loci: inter_seq})
            if stop > inter_start:
                inter_start = stop
        if int(value[1]) - stop >= minlen: # Get IR at end of dna_region if longer than minlen
            stop -= exlen
            inter_seq = seq[stop:int(value[1])]
            inter_loci = str(stop) + '_' + str(value[1])
            intergenic_regions.update({inter_loci: inter_seq})

        dna_regions.update({key: (value[0], value[1], value[2], intergenic_regions)})

    write_fasta(dna_regions, ident, out_filename)
    write_gff(dna_regions, ident, out_filename)

if __name__ == "__main__":
    comparator(**vars(options))
