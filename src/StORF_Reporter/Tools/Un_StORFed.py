import argparse
import collections
from datetime import date
import gzip
import sys



###################
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
    return crick
############################

#Output FASTA and GFF separately using the same out_filename but with respective extensions - gz output optional
def write_fasta(UnStORFed_Regions, options):
    if options.out_file == None:
        options.out_file = options.gff.split('.gff')[0]
    if options.gz == False:
        out = open(options.out_file + '_UnStORFed.fasta','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UnStORFed.fasta.gz', 'wt', newline='\n', encoding='utf-8')

    for pos, seq in UnStORFed_Regions.items():
        out.write('>UnStORFed_' + pos + '\n' + seq + '\n')
    out.close()

def write_gff(UnStORFed_Regions,options):
    if options.out_file == None:
        options.out_file = options.gff.split('.gff')[0]
    if options.gz == False:
        out =  open(options.out_file + '_UnStORFed.gff','w', newline='\n', encoding='utf-8')
    elif options.gz == True:
        out = gzip.open(options.out_file + '_UnStORFed.gff.gz', 'wt', newline='\n', encoding='utf-8')
    out.write("##gff-version\t3\n#\tUnStORFed Regions \n#\tRun Date:" + str(date.today()) + '\n')
    out.write("##Original Files: " + options.fasta + ' | ' + options.gff + '\n')
    for pos, seq in UnStORFed_Regions.items():
        length = len(seq)
        pos = pos.replace('_',' ')
        entry = ('Genome\tUnStORFed\tunannotated_region\t' + pos + '\t.\t.\t.\tID=UnStORFed_Region_'+pos+';Note=Length:' + str(length) + '\n')
        out.write(entry)
    out.close()

def fasta_load(fasta_in):
    dna_regions = collections.OrderedDict()
    first = True
    #### Default for when presented with standard fasta file
    for line in fasta_in:
        line = line.strip()
        if line.startswith('#'):
            continue
        elif line.startswith('>') and first == False:  # Check if first seq in file
            dna_regions.update({dna_region_id: [seq, []]})
            seq = ''
            dna_region_id = line.split()[0].replace('>', '').split('_UR_')[1]
        elif line.startswith('>'):
            seq = ''
            dna_region_id = line.split()[0].replace('>', '').split('_UR_')[1]
        else:
            seq += str(line)
            first = False
    dna_regions.update({dna_region_id: [seq, []]})
    return dna_regions

def gff_load(options,gff_in,dna_regions):
    #Will code in different versions for different types of GFF3 files (Prodigal,Ensembl etc)
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n') or line.startswith('#'): # Not to crash on empty lines in GFF
            continue
        elif line_data[1] == 'StORF_Reporter':
            if 'Con-StORF' not in line:
                current_UR = line_data[8].split('_UR_')[1].split('_StORF')[0]
            else:
                current_UR = line_data[8].split('_UR_')[1].split('_Con-StORF')[0]
            current_UR_StORF = line_data[3] + '_' + line_data[4]
            current_UR_StORF_Frame = line_data[6]
            dna_regions[current_UR][1].append([current_UR_StORF,current_UR_StORF_Frame])
## append storf info to this dict and then go to another func where if theres one storf we just extract normally, if more than one, we likely have to use numpy set info
    return dna_regions

def extractor(options):
    try:
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
    except AttributeError:
        sys.exit("Something Happened")

    UnStORFed_Regions = collections.defaultdict(str)

    for (key,(seq,StORFs))  in dna_regions.items(): #Extract URs from 1 dna_region at a time
        UR_start = int(key.split('_')[0])
        UR_stop = int(key.split('_')[1])
        UR_length = len(seq)
        if len(StORFs) == 1:
            StORF_start = int(StORFs[0][0].split('_')[0])
            StORF_stop = int(StORFs[0][0].split('_')[1])
            relative_StORF_distance_from_start = StORF_start - UR_start
            relative_StORF_distance_from_end = UR_stop - StORF_stop
            if relative_StORF_distance_from_start >= options.minlen:
                UnStORFed_pos = str(UR_start) + '_' + str(StORF_start)
                UnStORFed_start_seq = seq[0:relative_StORF_distance_from_start-1] # -1 is to remove the first bp of the first stop codon of the StORF
                UnStORFed_Regions[UnStORFed_pos] = UnStORFed_start_seq
            if relative_StORF_distance_from_end >= options.minlen:
                UnStORFed_pos = str(StORF_stop) + '_' + str(UR_stop)
                UnStORFed_stop_seq = seq[UR_length-relative_StORF_distance_from_end:]
                UnStORFed_Regions[UnStORFed_pos] = UnStORFed_stop_seq

        elif len(StORFs) > 2:
            prev_StORF_end = 0
            for idx, StORF in enumerate(StORFs):
                StORF_start = int(StORF[0].split('_')[0])
                StORF_stop = int(StORF[0].split('_')[1])

                try:  # Either use UR boundries or the next StORFs boundries
                    next_StORF_start = int(StORFs[idx + 1][0].split('_')[0])
                    relative_StORF_distance_from_start = max(StORF_start - prev_StORF_end, 0)
                    relative_StORF_distance_from_end = max(next_StORF_start - StORF_stop, 0)
                    segment_end = next_StORF_start
                    seq_segment_end = next_StORF_start - UR_start
                except IndexError:  # Next StORF is last StORF
                    relative_StORF_distance_from_start = max(StORF_start - prev_StORF_end, 0)
                    relative_StORF_distance_from_end = UR_stop - StORF_stop
                    segment_end = UR_stop
                    seq_segment_end = -1
                if idx == 0: # sloppy and overwriting from above but its needed for now
                    relative_StORF_distance_from_start = StORF_start - UR_start
                    relative_StORF_distance_from_end = max(next_StORF_start - StORF_stop, 0)
                    segment_start = UR_start
                    seq_segment_start = 0

                else:
                    segment_start = prev_StORF_end
                    seq_segment_start = prev_StORF_end - UR_start


                if relative_StORF_distance_from_start >= options.minlen:
                    UnStORFed_pos = str(segment_start) + '_' + str(StORF_start)
                    UnStORFed_start_seq = seq[seq_segment_start:seq_segment_start+relative_StORF_distance_from_start - 1]  # -1 is to remove the first bp of the first stop codon of the StORF
                    UnStORFed_Regions[UnStORFed_pos] = UnStORFed_start_seq
                if relative_StORF_distance_from_end >= options.minlen:
                    UnStORFed_pos = str(StORF_stop) + '_' + str(segment_end)
                    UnStORFed_stop_seq = seq[seq_segment_end - relative_StORF_distance_from_end:seq_segment_end]
                    UnStORFed_Regions[UnStORFed_pos] = UnStORFed_stop_seq

                prev_StORF_end = StORF_stop

    write_fasta(UnStORFed_Regions, options)
    write_gff(UnStORFed_Regions, options)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF_Reporter v0.5.4: UR_StORFed Run Parameters.')
    parser.add_argument('-f', '--fasta_seq', action='store', dest='fasta', required=True,
                        help='FASTA file for Unannotated Region seq extraction')
    parser.add_argument('-gff', action='store', dest='gff', help='GFF annotation file for the FASTA',
                        required=True)
    parser.add_argument('-min_len', action='store', dest='minlen', default='50', type=int,
                        help='Minimum Un-StORFed UR Length: Default 50')
    parser.add_argument('-o', '--output_file', action='store', dest='out_file', required=False,
                        help='Output file - Without filetype - default appends "_UnStORFed.fasta/.gff" to end of input gff filename')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-v', action='store', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')

    options = parser.parse_args()
    extractor(options)

    # Contig name could have a ';' which will mess up later on in StORF_Reporter-R
    # UR output should state original non extended

    ## The current version reports UnStORFed regions according to the extracted URs - So the