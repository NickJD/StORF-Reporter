import argparse
from datetime import date
import re

parser = argparse.ArgumentParser(description='StORF Run Parameters.')
parser.add_argument('-seq', action="store", dest='fasta',required=True,
                    help='Input Sequence File')
parser.add_argument('-ir',dest='intergenic',action='store',default='True', type=eval, choices=[True,False],
                    help='Default - Treat input as Intergenic: Use "-ir False" for standard fasta')
parser.add_argument('-wc', action="store", dest='whole_contig',default='False', type=eval, choices=[True,False],
                    help='Default - False: StORFs reported across entire sequence')
parser.add_argument('-ps', action="store", dest='partial_storf',default='False', type=eval, choices=[True,False],
                    help='Default - False: Partial StORFs reported')
parser.add_argument('-filt',action='store',dest='filtering',default='hard', const='hard', nargs='?', choices=['none', 'soft', 'hard'],
                     help='Default - Hard: Filtering level none is not recommended, soft for single strand filtering '
                          'and hard for both-strand longest-first tiling')
parser.add_argument('-aa', action="store", dest='translate',default='False',type=eval, choices=[True,False],
                    help='Default - False: Report StORFs as amino acid sequences')
parser.add_argument('-minorf', action="store", dest='min_orf', default=100, type=int,
                    help='Default - 100: Minimum StORF size in nt')
parser.add_argument('-maxorf', action="store", dest='max_orf', default=99999, type=int,
                    help='Default - 99999: Maximum StORF size in nt')
parser.add_argument('-codons', action="store", dest='stop_codons', default="TAG,TGA,TAA",
                    help='Default - ("TAG,TGA,TAA"): List Stop Codons to use')
parser.add_argument('-olap', action="store", dest='overlap_nt', default=20, type=int,
                    help='Default - 20: Maximum number of nt of a StORF which can overlap another StORF.')
parser.add_argument('-gff',action='store',dest='gff', default='True', type=eval, choices=[True,False],
                    help='Default - True: StORF Output a GFF file')
parser.add_argument('-o', action="store", dest='output',required=True,
                    help='Output file name - Without filetype')
options = parser.parse_args()
fasta = options.fasta
intergenic = options.intergenic
whole_contig = options.whole_contig
partial_storf = options.partial_storf
filtering = options.filtering
translate = options.translate
min_orf_size = options.min_orf
max_orf_size = options.max_orf
stop_codons = options.stop_codons
overlap_nt= options.overlap_nt
gff = options.gff
output = options.output
import collections
STORF_Num = 0
#################################
def revCompIterative(watson): #Gets Reverse Complement
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev:
        crick += complements[nt]
    return crick

def tile_filtering(storfs): #Hard filtering
    storfs = collections.OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
    ################ - Order largest first filtering
    if filtering == 'hard':
        storfs = sorted(storfs.items(), key=lambda storfs:storfs[1][3],reverse=True)
        ordered_by_length = collections.OrderedDict()
        print("StORFS Ordered by Length: " + str(len(ordered_by_length)))
        import time
        start = time.time()
        for tup in storfs:
            ordered_by_length.update({tup[0]:tup[1]})
        ############## - For each STORF, remove all smaller overlapping STORFs according to filtering rules
        length = len(ordered_by_length)
        i = 0
        ordered_by_length = list(ordered_by_length.items())
        while i < length:
            pos_x, data_x = ordered_by_length[i]
            start_x = int(pos_x.split(',')[0])
            stop_x = int(pos_x.split(',')[1])
            j = i+1
            while j < length:
                pos_y, data_y = ordered_by_length[j]
                start_y = int(pos_y.split(',')[0])
                stop_y = int(pos_y.split(',')[1])
                if start_y >= stop_x or  stop_y <= start_x:
                    j+=1
                    continue  # Not caught up yet / too far
                elif start_y >= start_x and stop_y <= stop_x:
                    ordered_by_length.pop(j)
                    length = len(ordered_by_length)
                else:
                    x = set(range(start_x,stop_x+1))
                    y = set(range(start_y,stop_y+1))
                    overlap = len(x.intersection(y))
                    print(overlap)
                    if overlap >= 50: # Change ehere
                        ordered_by_length.pop(j)
                        length = len(ordered_by_length)
                    else:
                        j += 1
            length = len(ordered_by_length)
            i+=1
        print("StORFs After Filtering: " + str(len(ordered_by_length)))
        end = time.time()
        print('Time taken for filtering: ', end - start)
        storfs = collections.OrderedDict(ordered_by_length)
    return storfs

def translate_frameshifted(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

def write_gff(storfs,seq_id):
    with open(output+'.gff','a') as out:
        for k, v in storfs.items():
            sequence = v[0]
            start_stop = sequence[0:3]
            end_stop = sequence[-3:]
            pos = k.split(',')
            start = int(pos[0])
            stop = int(pos[1])
            frame = int(v[1])
            seq_id = seq_id.split()[0].replace('>','')
            storf_name = seq_id + '_' + str(list(storfs.keys()).index(k))
            length = stop - start
            if intergenic == True:
                ir_start = start + int(storf_name.split('|')[1].split('_')[0])
                ir_stop = stop + int(storf_name.split('|')[1].split('_')[1])
                entry = (seq_id + '\tSTORF\tCDS\t' + str(ir_start) + '\t' + str(ir_stop) + '\t.\t' + v[
                    2] + '\t.\tID=' + storf_name + ':IR_Position=' + str(start) + '_' + str(stop) + ':Length=' + str(
                    length) + ':Frame=' + str(frame)+  ':Start_Stop='+start_stop+':End_Stop='+end_stop+ '\n')
            elif intergenic == False:
                entry = (seq_id + '\tSTORF\tCDS\t' + str(start) + '\t' + str(stop) + '\t.\t' + v[
                    2] + '\t.\tID=' + storf_name + ':Length=' + str(
                    length) + ':Frame=' + str(frame) +':Start_Stop='+start_stop+':End_Stop='+end_stop+'\n')
            out.write(entry)

def write_fasta(storfs,seq_id,genome_size):
    storf_num = 0
    with open(output+'.fasta','a') as out:
        for k, v in storfs.items():
            strand = v[2]
            sequence = v[0]
            frame = int(v[1])
            start = int(k.split(',')[0])
            stop = int(k.split(',')[1])
            seq_id = seq_id.split()[0].replace('>','')
            storf_name = seq_id+ '_' + str(list(storfs.keys()).index(k))
            if intergenic == True:
                ir_start = start + int(storf_name.split('|')[1].split('_')[0])
                ir_stop = stop + int(storf_name.split('|')[1].split('_')[1])
                out.write(">"+str(seq_id)+'_'+str(storf_num)+"|"+str(ir_start) + strand + str(ir_stop) + "|Frame:"+str(frame)+"\n")
            elif intergenic == False:
                out.write(">"+str(seq_id)+'_'+str(storf_num)+"|"+str(start) + strand + str(stop) + "|Frame:"+str(frame)+"\n")
            if translate == True:
                if "+" in strand:
                    amino = translate_frameshifted(sequence[0:])
                    out.write(amino + '\n')
                if "-" in strand:
                    amino = translate_frameshifted(sequence[0:])
                    out.write(amino + '\n')
            elif translate == False:
                out.write(sequence+'\n')
            storf_num += 1
#            sys.stdout.flush()
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
      'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
      'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W'}
############################

def find_storfs(stops,sequence,storfs,frames_covered,counter,lengths,strand):
    first = True
    next_stops = []
    start_stops = []
    seen_stops = []
    for stop in stops:  # Finds Stop-Stop#
        seen_stops.append(stop)
        if strand == '+':
            frame = (stop % 3) + 1  # check
        elif strand == '-':
            frame = (stop % 3) + 4  # check

        frames_covered.update({frame: 1})


        for next_stop in stops[counter + 1:]:
            length = abs(next_stop - stop)
            if length % 3 == 0 and next_stop not in seen_stops:

                    if length >= min_orf_size and length <= max_orf_size:
                        if filtering == 'none': # This could be made more efficient
                            seq = sequence[stop:next_stop + 3]
                            storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length]})
                            lengths.append(length)

                        elif first == False:
                            if stop > prev_stop and next_stop < prev_next_stop:
                                break

                            storf = set(range(stop, next_stop + 4))
                            storf_overlap = len(prev_storf.intersection(storf))

                            if storf_overlap <= overlap_nt:
                                seq = sequence[stop:next_stop + 3]

                                storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length]})
                                seen_stops.append(next_stop)
                                prev_storf = storf
                                prev_stop = stop
                                prev_next_stop = next_stop
                                next_stops.append(next_stop)
                                start_stops.append(stop)
                                prevlength = prev_next_stop - prev_stop

                                break

                            elif storf_overlap >= overlap_nt:  # and length > prevlength:

                                if length > prevlength:
                                   
                                    storfs.popitem()
                                    seq = sequence[stop:next_stop + 3]
                                    storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length]})
                                    seen_stops.append(next_stop)
                                    prev_storf = storf
                                    prev_stop = stop
                                    prev_next_stop = next_stop
                                    next_stops.append(next_stop)
                                    start_stops.append(stop)
                                    prevlength = prev_next_stop - prev_stop

                                break

                            else:

                                break

                        elif first == True:
                            if partial_storf == True: # upstream partial StORF
                                if stop > min_orf_size:
                                    seq = sequence[0:stop + 3] #Start of seq to first stop identified
                                    storfs.update({",".join([str(0), str(stop)]): [seq, str(frame), strand, stop]})
                            seq = sequence[stop:next_stop + 3]
                            length = next_stop - stop
                            storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length]})
                            seen_stops.append(next_stop)
                            prev_storf = set(range(stop, next_stop + 4))
                            prev_stop = stop
                            prev_next_stop = next_stop
                            next_stops.append(next_stop)
                            start_stops.append(stop)
                            prevlength = prev_next_stop - prev_stop
                            first = False
                            break
                    else:
                        break
        counter +=1
    if partial_storf == True:  # downstream partial StORF - Last Stop to end of sequence
        if (len(sequence) - stop) > min_orf_size:
            seq = sequence[stop:len(sequence)]  # Start of seq to first stop identified
            storfs.update({",".join([str(stop), str(len(sequence))]): [seq, str(frame), strand, stop]})


    return storfs,frames_covered,counter,lengths

def STORF(sequence): #Main Function
    global storf_num
    global genome_size
    genome_size = len(sequence)
    storf_num = 0
    stops = []
    frames_covered = collections.OrderedDict()
    for x in range (1,7):
        frames_covered.update({x: 0})
    for stop_codon in stop_codons.split(','):#Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence)]
    stops.sort()
    storfs = collections.OrderedDict()
    counter = 0
    lengths = []
    storfs,frames_covered,counter,lengths = find_storfs(stops,sequence,storfs,frames_covered,counter,lengths,'+')
    ###### Reversed
    sequence_rev = revCompIterative(sequence)
    stops = []
    for stop_codon in stop_codons.split(','):#Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence_rev)]
    stops.sort()
    counter = 0
    storfs,frames_covered,counter,lengths = find_storfs(stops,sequence_rev,storfs,frames_covered,counter,lengths,'-')
    #Potential run-through StORFs
    if whole_contig == True:
        for strand,frame in frames_covered.values():
            if frame == 0:
                if strand <4:
                   storfs.update({",".join([str(0), str(len(sequence))]): [sequence, str(frame), '+', 0]})
                else:
                    storfs.update({",".join([str(0), str(len(sequence))]): [sequence, str(frame), '-', 0]})
    #Check if there are StORFs to report
    if bool(storfs):
        storfs = tile_filtering(storfs)
        # Reorder by start position
        storfs = collections.OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
        write_fasta(storfs, sequence_id, genome_size)
        write_gff(storfs, sequence_id)
    else:
        print("No StOFS Found")



if __name__ == "__main__":
    with open(output + '.gff', 'w') as out:
        out.write("##gff-version\t3\n#\tSTORF Stop - Stop CDS Predictions\n#\tRun Date:" + str(date.today()) + '\n')
    sequences = collections.OrderedDict()
    First = True
    file = open(fasta, "r")
    print(file.name)
    for line in file:
        line = line.strip()
        if '#' in line:
            continue
        elif ">" in line and First == False:
            sequences.update({sequence_name: seq})
            seq = ''
            sequence_name = line
        elif '>' in line:
            seq = ''
            sequence_name = line
        else:
            seq += str(line)
            First = False
    sequences.update({sequence_name: seq})
    for sequence_id, sequence in sequences.items():
        if len(sequence) >= min_orf_size:
            STORF(sequence)





