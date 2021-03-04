import argparse
import collections
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
        if dna_region_ir[3]:
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
        if dna_region_ir[3]:
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

def gff_load(gff_in,dna_regions):
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        #temp
        line = line.replace('NZ_','')
        line_data = line.split()
        if 'ID=gene' in options.gene_ident:
            if line_data[0] in dna_regions and options.gene_ident in line_data[8]:
                pos = line_data[3] + '_' + line_data[4]
                dna_regions[line_data[0]][2].append(pos) # This will add to list
        else:
            gene_types = options.gene_ident.split(',')
            if line_data[0] in dna_regions:
                if any(gene_type in line_data[2] for gene_type in gene_types):
                    print(line_data[2])
                    pos = line_data[3] + '_' + line_data[4]
                    dna_regions[line_data[0]][2].append(pos) # This will add to list

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

    for (key,(seq,seq_length,posns,irs))  in dna_regions.items(): #Extract IRs from 1 dna_region at a time
        intergenic_regions = collections.OrderedDict()
        inter_start = 0
        for pos in posns: # Iterate over GFF loci and measure flanking regions for potential IRs
            start = int(pos.split('_')[0])
            stop = int(pos.split('_')[1])
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
        try:
            if (seq_length - stop) >= options.minlen: # Get IR at end of dna_region if longer than minlen
                stop -= options.exlen
                inter_seq = seq[stop:seq_length]
                inter_loci = str(stop) + '_' + str(seq_length)
                intergenic_regions.update({inter_loci: inter_seq})
        except UnboundLocalError:
            pass
        dna_regions.update({key: (seq, seq_length, posns, intergenic_regions)})
        # if intergenic_regions:
        #     dna_regions.update({key: (seq, seq_length, posns, intergenic_regions)})
        # else:
        #     del dna_regions[key]



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
    parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene',
                        help='Identifier used for extraction of "genic" regions ("CDS,rRNA,tRNA"): Default for Ensembl_Bacteria = "ID=gene"')
    parser.add_argument('-o', '--output_prefix', action='store', dest='out_prefix', required=True,
                        help='Output file prefix - Without filetype')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')

    options = parser.parse_args()
    comparator(options)

    # Contig name could have a ';' which will mess up later on in StORF-R