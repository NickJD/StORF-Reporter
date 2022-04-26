import argparse
import collections
from datetime import date
import gzip


#Output FASTA and GFF separately using the same out_filename but with respective extensions - gz output optional
def write_fasta(dna_regions, options):
    if options.out_file == None:
        options.out_file = options.fasta.split('.')[0]
    if options.gz == False:
        out =  open(options.out_file + '_UR.fasta','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UR.fasta.gz', 'wt', newline='\n', encoding='utf-8')
    for dna_region, dna_region_ur in dna_regions.items():
        ur_ident = dna_region + options.ident # Add user ident onto name of dna regions
        if dna_region_ur[3]:
            for ir, ur_seq in dna_region_ur[3].items():
                out.write('>' + ur_ident + '|' + ir + '\n' + ur_seq + '\n')
    out.close()

def write_gff(dna_regions,options):
    if options.out_file == None:
        options.out_file = options.fasta.split('.')[0]
    if options.gz == False:
        out =  open(options.out_file + '_UR.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UR.gff.gz', 'wt', newline='\n', encoding='utf-8')
    out.write("##gff-version\t3\n#\tUR Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    out.write("##Original File: " + options.fasta + '\n')
    for dna_region, dna_region_ur in dna_regions.items():
        ur_ident = dna_region + options.ident
        if dna_region_ur[3]:
            for ir, ur_seq in dna_region_ur[3].items():
                length = len(ur_seq)
                ur_pos = ir.replace('_','\t')
                entry = (dna_region + '\tUR_Extractor\tunannotated_region\t' + ur_pos + '\t.\t.\t.\tID='+ur_ident+'_'+ir+';Note=UR_Length(Extended):' + str(length) + '\n')
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

def gff_load(options,gff_in,dna_regions):
    #Will code in different versions for different types of GFF3 files (Prodigal,Ensembl etc)
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n'): # Not to crash on empty lines in GFF
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
    try: # Detect whether fasta/gff files are .gz or text and read accordingly
        fasta_in = gzip.open(options.fasta,'rt')
        dna_regions = fasta_load(fasta_in)
    except:
        fasta_in = open(options.fasta,'r')
        dna_regions = fasta_load(fasta_in)
    try:
        gff_in = gzip.open(options.gff,'rt')
        dna_regions = gff_load(options,gff_in,dna_regions)
    except:
        gff_in = open(options.gff,'r')
        dna_regions = gff_load(options,gff_in,dna_regions)

    for (key,(seq,seq_length,posns,irs))  in dna_regions.items(): #Extract IRs from 1 dna_region at a time
        intergenic_regions = collections.OrderedDict()
        inter_start = 0
        if posns: # If IR has a pos
            for pos in posns: # Iterate over GFF loci and measure flanking regions for potential IRs
                start = int(pos.split('_')[0])
                stop = int(pos.split('_')[1])
                ###### This hack is to get over GFF errors where genome-long annotations
                if stop-start >= 100000:
                    continue
                if start > inter_start:
                    length = start - inter_start
                    if length >= options.minlen and length <= options.maxlen: #default between 30 - 100,000
                        inter_start = max(inter_start - options.exlen, 0)
                        start += options.exlen
                        inter_seq = seq[inter_start:start]
                        inter_loci = str(inter_start) + '_' + str(start)
                        intergenic_regions.update({inter_loci: inter_seq})
                if stop > inter_start:
                    inter_start = stop
            try:
                if (seq_length - stop) >= options.minlen and (seq_length - stop) <= options.maxlen: #default between 30 - 100,000: # Get IR at end of dna_region if longer than minlen
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


    if options.nout == False:
        write_fasta(dna_regions, options)
        write_gff(dna_regions, options)
    else:
        return dna_regions
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta', required=True,
                        help='FASTA file for Unannotated Region seq extraction')
    parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=True)
    parser.add_argument('-ident', action='store', dest='ident', default='_UR',
                        help='Identifier given for Unannotated Region output sequences: Default "Input"_UR')
    parser.add_argument('-min_len', action='store', dest='minlen', default='30', type=int,
                        help='Minimum UR Length: Default 30')
    parser.add_argument('-max_len', action='store', dest='maxlen', default='100000', type=int,
                        help='Maximum UR Length: Default 100,000')
    parser.add_argument('-ex_len', action='store', dest='exlen', default='50', type=int,
                        help='UR Extension Length: Default 50')
    parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='ID=gene',
                        help='Identifier used for extraction of "unannotated" regions "CDS,rRNA,tRNA": Default for Ensembl_Bacteria = "ID=gene"')
    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=False,
                        help='Output file - Without filetype - default appends "_UR" to end of input fasta')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-v', action='store', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')
    parser.add_argument('-nout', action='store', dest='nout', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)

    options = parser.parse_args()
    extractor(options)

    # Contig name could have a ';' which will mess up later on in StORF-Reporter-R
    # UR output should state original non extended