import sys
import argparse
from datetime import date
import re

parser = argparse.ArgumentParser(description='STORF Run Parameters.')

parser.add_argument('-seq', action="store", dest='sequence',required=True,
                    help='Input Sequence File')
parser.add_argument('-ir',dest='Intergenic',action='store_true',default=False,
                    help='Treat input as a Intergenic Region file')
parser.add_argument('-std',dest='Standard',action='store_true',default=False,
                    help='Treat input as a standard FASTA')
parser.add_argument('-WC', action="store", dest='WholeContig',default=False, type=bool,
                    help='Should STORFs covering entire contig be reported')
parser.add_argument('-PS', action="store", dest='PartialSTORF',default=False, type=bool,
                    help='Should partial STORFs be reported')
parser.add_argument('-filtering', const='hard',dest='filtering',nargs='?',default='hard',choices=['none','soft','hard'],
                    help='Choose filtering level: none is not recommended, soft for strand-only filtering and hard for longest-first tiling filtering')
parser.add_argument('-AA', action="store", dest='translate',default=False, type=bool,
                    help='Should STORFs beconverted into Amino Acid sequences')
parser.add_argument('-MinORF', action="store", dest='MinORF', default=100, type=int,
                    help='Minimum STORF Size')
parser.add_argument('-MaxORF', action="store", dest='MaxORF', default=99999, type=int,
                    help='Maximum STORF Size')
parser.add_argument('-SCodons', action="store", dest='StopCodons', default="TAG,TGA,TAA",
                    help='List Stop Codons to use ("TAG,TGA,TAA")')

parser.add_argument('-OverlapNT', action="store", dest='OverlapNT', default=20, type=int,
                    help='Maximum Number of NTs of a STORF which can overlap another STORF.')
parser.add_argument('-GFF_Out',action='store',dest='GFF', default=True, type=bool,
                    help='Should STORF Output a GFF file - Default True')
parser.add_argument('-o', action="store", dest='Output',
                    help='Output File Name')


options = parser.parse_args()


sequence = options.sequence
intergenic = options.Intergenic
standard = options.Standard
whole_contig = options.WholeContig
partial_storf = options.PartialSTORF
filtering = options.filtering
translate = options.translate
min_orf_size = options.MinORF
max_orf_size = options.MaxORF
stop_codons = options.StopCodons
overlap_nt= options.OverlapNT
GFF = options.GFF
output = options.Output
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

    #both_frames = collections.OrderedDict(list(frames.items()) + list(Frames_rev.items()))
    storfs = collections.OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))

    ################ - Order largest first filtering

    if filtering == 'hard':
        storfs = sorted(storfs.items(), key=lambda storfs:storfs[1][3],reverse=True)
        ordered_by_length = collections.OrderedDict()
        print("Num Ordered by Length " + str(len(ordered_by_length)))

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

        print("Num After Filtering " + str(len(ordered_by_length)))
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
            #What STORF should do is give original location and IR specific location.- What about the strand for original and IR specific?
            length = stop - start
            if intergenic == True:
                ir_start = start + int(storf_name.split('_')[2].split('|')[1])
                ir_stop = stop + int(storf_name.split('_')[2].split('|')[1])
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
                ir_start = start + int(storf_name.split('_')[2].split('|')[1])
                ir_stop = stop + int(storf_name.split('_')[2].split('|')[1])
                out.write(">"+str(seq_id)+'_'+str(storf_num)+"|"+str(ir_start) + "_" + strand + "_" + str(ir_stop) + "|Frame:"+str(frame)+"\n")
            elif intergenic == False:
                out.write(">"+str(seq_id)+'_'+str(storf_num)+"|"+str(start) + "_" + strand + "_" + str(stop) + "|Frame:"+str(frame)+"\n")



            if translate == True:
                start = int(k.split(",")[0])
                stop = int(k.split(",")[1])
                if frame <= 3:
                    if start != 0 and stop != genome_size:
                        frame = 0
                    elif genome_size == stop:
                        frame = 0
                    else:
                        frame = abs(stop % 3) # Need to check this!!!!
                elif frame > 3:
                    if start != 0 and stop != genome_size:
                        frame = 0
                    elif genome_size == stop:
                        frame = 0
                    else:
                        frame = abs(stop % 3) ######

                if "+" in strand:
                    #Fr -= 1#select frame

                    amino = translate_frameshifted(sequence[0:])
                    out.write(amino + '\n')
                if "-" in strand:
                    #Fr -= 4
                    amino = translate_frameshifted(sequence[0:])
                    out.write(amino + '\n')
            elif translate == False:
                out.write(sequence+'\n')
            storf_num += 1

            sys.stdout.flush()

# def writeOutPreA(Frames,Contig_ID,Di, End):
#     global STORF_Num
#     with open(output,'a') as out:
#         for k, v in Frames.iteritems():
#             Contig_ID = Contig_ID.replace('>', '')
#             out.write(">"+str(STORF_Num)+":"+str(Contig_ID)+','+str(k.split(',')[0]) + "\t" + Di + "\t" + str(k.split(',')[1]) + ",Frame: "+str(v[1])+"\n")
#             out.write(v[0]+'\n')
#             STORF_Num += 1
#             print STORF_Num
#             sys.stdout.flush()


#
# def run_Off(stops,nextStops,startStops,Contig,Frames,Genome_Size,Frames_Covered,po,Di): # Run off - STORF which start or stop n the middle of the Cotnig - Currnently has to be smaller than max size.
#     frames = (0,1,2)
#     if "start" in po:
#         runOff = []
#         runOff_Longest = 0
#         runOff_Frame = 0
#         runOff_Run = False
#
#         for stop in stops:
#             Frame = stop % 3
#             if stop not in nextStops and Frame not in runOff:
#                 runOff.append(Frame)
#                 if stop > runOff_Longest:
#                     runOff_Longest = stop
#                     runOff_Frame = Frame
#                     runOff_Run = True
#                 if all(elem in runOff for elem in frames):
#                     break
#         if runOff_Run == True and runOff_Longest <= max_orf_size and runOff_Longest >= 30:
#             seq = Contig[0:runOff_Longest + 3]
#
#             if '-' in Di:
#                 Frame +=4
#                 frameCheck_rev(Frames, Frames_Covered, Frame)
#                 runOff_Frame +=4
#                 Frames.update({",".join(['0', str(runOff_Longest)]): [seq, str(runOff_Frame), Di,runOff_Longest]})
#             elif '+' in Di:
#                 Frame +=1
#                 frameCheck(Frames, Frames_Covered, Frame)
#                 runOff_Frame +=1
#                 Frames.update({",".join(['0', str(runOff_Longest)]): [seq, str(runOff_Frame), Di,runOff_Longest]})
#
#
#
#
#     elif "end" in po:
#         runOff = []
#         runOff_Longest = Genome_Size
#         runOff_Frame = 0
#         runOff_Run = False
#
#         for stop in stops:
#             Frame = stop % 3
#             if stop not in startStops and Frame not in runOff:
#                 runOff.append(Frame)
#
#                 if stop < runOff_Longest:
#                     runOff_Longest = stop
#                     runOff_Frame = Frame
#                     runOff_Run = True
#
#                 if all(elem in runOff for elem in frames):
#                     break
#         length = abs(runOff_Longest - Genome_Size)
#         if runOff_Run == True and length <= max_orf_size and length >= 30:
#             seq = Contig[runOff_Longest: Genome_Size]
#
#             if '-' in Di:
#                 Frame +=4
#                 frameCheck_rev(Frames, Frames_Covered, Frame)
#                 runOff_Frame +=4
#                 Frames.update({",".join([str(runOff_Longest), str(Genome_Size)]): [seq, str(runOff_Frame), Di,length]})
#             elif '+' in Di:
#                 Frame +=1
#                 frameCheck(Frames, Frames_Covered, Frame)
#                 runOff_Frame +=1
#                 Frames.update({",".join([str(runOff_Longest), str(Genome_Size)]): [seq, str(runOff_Frame), Di,length]})
#




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
##Change storf into 1 method
#remove framecheck function
#none, soft and hard filtering

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


    # if partial_storf == True:
    #     run_Off(stops,nextStops,startStops,Contig_rev,Frames_rev,genome_size,Frames_Covered,"start","-")
    #     stops.reverse()
    #     nextStops.reverse()
    #     run_Off(stops,nextStops,startStops,Contig_rev,Frames_rev,genome_size,Frames_Covered,"end","-")

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
    for stop_codon in stop_codons.split(','):#Find all Stops in Contig
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence)]


    stops.sort()

    storfs = collections.OrderedDict()


    counter = 0

    lengths = []

    storfs,frames_covered,counter,lengths = find_storfs(stops,sequence,storfs,frames_covered,counter,lengths,'+')


    #########################################
    ###### Reversed
    sequence_rev = revCompIterative(sequence)
    stops = []

    for stop_codon in stop_codons.split(','):#Find all Stops in Contig
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence_rev)]
    stops.sort()

    counter = 0


    storfs,frames_covered,counter,lengths = find_storfs(stops,sequence_rev,storfs,frames_covered,counter,lengths,'-')




#########################################################

#########################################################
    # if wholeContig == True: # This is broken - It is just changing the AA seq every time and not actually updating.
    #     if len(Contig) >= minOrfSize  and any(v == 0 for v in Frames_Covered.values()):  # all(v == 0 for v in Frames_Num.values()) and wholeContig == "True": # Mix up where it prints what ever the current contig is rather than the pa\rticular frame is or something
    #
    #         if Translate == True:  # PreA - Change name to something ekse  - supposed to state that something elke
    #             translated = translate_frameshifted(Contig[0:])
    #             Frames.update({",".join(["0", str(len(Contig))]): [translated, "1"]})
    #             #writeOutPreA(Frames, Contig_ID, "+", End)
    #             translated = translate_frameshifted(Contig[1:])
    #             Frames.update({",".join(["0", str(len(Contig))]): [translated, "2"]})
    #             #writeOutPreA(Frames, Contig_ID, "+", End)
    #             translated = translate_frameshifted(Contig[2:])
    #             Frames.update({",".join(["0", str(len(Contig))]): [translated, "3"]})
    #             #writeOutPreA(Frames, Contig_ID, "+", End)
    #             translated = translate_frameshifted(Contig_rev[0:])
    #             Frames_rev.update({",".join(["0", str(len(Contig))]): [translated, "4"]})
    #             #writeOutPreA(Frames, Contig_ID, "-", End)
    #             translated = translate_frameshifted(Contig_rev[1:])
    #             Frames_rev.update({",".join(["0", str(len(Contig))]): [translated, "5"]})
    #             #writeOutPreA(Frames, Contig_ID, "-", End)
    #             translated = translate_frameshifted(Contig_rev[2:])
    #             Frames_rev.update({",".join(["0", str(len(Contig))]): [translated, "6"]})
    #             #writeOutPreA(Frames, Contig_ID, "-", End)
    #         elif Translate == False:
    #             for x in range(0, 3):
    #                 Frames.update({",".join(["0", str(len(Contig))]): [Contig, "0"]})
    #                 #writeOut(Frames, Contig_ID, "+", End)
    #
    #             for x in range(0, 3):
    #                 Frames_rev.update({",".join(["0", str(len(Contig))]): [Contig_rev, "0"]})
    #                 #writeOut(Frames, Contig_ID, "-", End)

    ###########################################


    ##################################################




    storfs = tile_filtering(storfs)
        #Reorder by start position
    storfs = collections.OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))

    write_fasta(storfs,sequence_id,genome_size)
    write_gff(storfs,sequence_id)


        # storf_nuc_count = np.zeros((genome_size), dtype=np.int)
        # for pos in both_frames.keys():
        #
        #     storf_start = int(pos.split(',')[0])
        #     storf_stop = int(pos.split(',')[1])
        #     storf_nuc_count[storf_start:storf_stop] = [1] * (storf_stop - storf_start)  # Changing all between the two positions to 1's
        # gene_coverage_genome = 100 * float(np.count_nonzero(storf_nuc_count)) / float(genome_size)
        # print(gene_coverage_genome)





with open(output+'.gff','w') as out:
    out.write("##gff-version\t3\n#\tSTORF Stop - Stop CDS Predictions\n#\tRun Date:"+str(date.today())+'\n')
# with open(output+'.fasta','w') as out:
#     out.write("##STORF Stop - Stop CDS Predictions \t Run Date:"+str(date.today())+'\n')
#
#

seq=''
Contigs = collections.OrderedDict()
Contig_Name = ''
First = True
file = open(sequence, "r")
print(file)
for line in file:
    line=line.strip()
    if '#' in line:
        continue
    elif ">" in line and First == False:
        Contigs.update({Contig_Name:seq})
        seq=''
        Contig_Name = line
    elif '>' in line:
        seq = ''
        Contig_Name = line
    else:
        seq+=str(line)
        First = False
Contigs.update({Contig_Name:seq})

for sequence_id, sequence in Contigs.items():

    STORF(sequence)



