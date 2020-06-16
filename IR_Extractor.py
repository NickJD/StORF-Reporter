import collections
import argparse
from datetime import date
import gzip


#Output FASTA and GFF separately using the same out_filename but with respective extensions - gz output optional
def write_fasta(dna_regions, options):
    if options.gz == False:
        out =  open(options.out_prefix + '.fasta','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_prefix + '.fasta.gz', 'wt', newline='\n', encoding='utf-8')
    for dna_region, dna_region_ir in dna_regions.items():
        ir_ident = dna_region + options.ident # Add user ident onto name of dna regions
        for ir, ir_seq in dna_region_ir[3].items():
            out.write('>' + ir_ident + '|' + ir + '\n' + ir_seq + '\n')
    out.close()

def write_gff(dna_regions,options):
    if options.gz == False:
        out =  open(options.out_prefix + '.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_prefix + '.gff.gz', 'wt', newline='\n', encoding='utf-8')
    out.write("##gff-version\t3\n#\tIR Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    out.write("##Original File: " + options.fasta + '\n')
    for dna_region, dna_region_ir in dna_regions.items():
        ir_ident = dna_region + options.ident
        for ir, ir_seq in dna_region_ir[3].items():
            length = len(ir_seq)
            ir_pos = ir.replace('_','\t')
            entry = (dna_region + '\tIR_Extractor\tintergenic_region\t' + ir_pos + '\t.\t.\t.\tID='+ir_ident+'_'+ir+';Note=IR_Length(Extended):' + str(length) + '\n')
            out.write(entry)
    out.close()

def fasta_load(fasta_in):
    first = True
    dna_regions = collections.OrderedDict()
    for line in fasta_in:
        line = line.strip()
        if ">" in line and first == False:  # Check if first seq in file
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
    return dna_regions

def gff_load(gff_in,dna_regions):
    for line in gff_in:         # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line_data[0] in dna_regions and options.gene_ident in line:
            pos = line_data[3] + '_' + line_data[4]
            dna_regions[line_data[0]][-1].append(pos)
    return dna_regions
def comparator(options):
    try: # Detect whether fasta/gff files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta,'rt')
        dna_regions = fasta_load(fasta_in)
    except:
        fasta_in = open(options.fasta,'r')
        dna_regions = fasta_load(fasta_in)
    try:
        gff_in = gzip.open(options.gff,'rt')
        dna_regions = gff_load(gff_in,dna_regions)
    except:
        gff_in = open(options.gff,'r')
        dna_regions = gff_load(gff_in,dna_regions)

    for key, value in dna_regions.items(): #Extract IRs from 1 dna_region at a time
        intergenic_regions = collections.OrderedDict()
        inter_start = 0
        for pos in value[2]: # Iterate over GFF loci and measure flanking regions for potential IRs
            start = int(pos.split('_')[0])
            stop = int(pos.split('_')[1])
            seq = value[0]
            if start > inter_start:
                length = start - inter_start
                if length >= options.minlen:
                    inter_start = max(inter_start - options.exlen, 0)
                    start += options.exlen
                    inter_seq = seq[inter_start:start]
                    inter_loci = str(inter_start) + '_' + str(start)
                    intergenic_regions.update({inter_loci: inter_seq})
            if stop > inter_start:
                inter_start = stop
        if (int(value[1]) - stop) >= options.minlen: # Get IR at end of dna_region if longer than minlen
            stop -= options.exlen
            inter_seq = seq[stop:int(value[1])]
            inter_loci = str(stop) + '_' + str(value[1])
            intergenic_regions.update({inter_loci: inter_seq})

        dna_regions.update({key: (value[0], value[1], value[2], intergenic_regions)})

    write_fasta(dna_regions, options)
    write_gff(dna_regions, options)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta', required=True,
                        help='FASTA file for Intergenic Region seq extraction')
    parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=True, )
    parser.add_argument('-ident', action='store', dest='ident', default='_IR',
                        help='Identifier given for Intergenic Region output sequences: Default "Input"_IR')
    parser.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                        help='Minimum IR Length: Default 30')
    parser.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                        help='IR Extension Length: Default 50')
    parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene:',
                        help='Identifier used for extraction: Default = ID=gene:')
    parser.add_argument('-o', '--output_prefix', action='store', dest='out_prefix', required=True,
                        help='Output file prefix - Without filetype')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')

    options = parser.parse_args()
    comparator(options)

