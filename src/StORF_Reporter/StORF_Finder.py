import argparse
import collections
import re
from collections import defaultdict, OrderedDict
from datetime import date
import textwrap
import gzip
import os
import sys
#from line_profiler_pycharm import profile


try:
    from .utils import sortORFs  # Calling from ORForise via pip
    from .constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from utils import sortORFs
    from constants import *




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
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate



############################
def check_non_standard_chars(options,storf_seq):
    non_standard_set = set(['A', 'T', 'G', 'C'])
    total_chars = len(storf_seq)
    non_standard_count = sum(1 for char in storf_seq if char not in non_standard_set)
    non_standard_percentage = non_standard_count / total_chars
    if non_standard_percentage > float(options.non_standard):
        return False
    else:
        return True


############################

def reverseCorrectLoci(options,sequence_length,sequence_id,first,second,third): # here for the negative loci correction
    if second == None:
        corrected_start = max(sequence_length - int(third-3),1)
        corrected_stop = max(sequence_length - int(first-3),1)
        return corrected_start, corrected_stop
    else:
        corrected_start = max(sequence_length - int(third-3),1)
        corrected_mid = max(sequence_length - int(second-3),1)
        corrected_stop = max(sequence_length - int(first-3),1)
        return corrected_start, corrected_mid, corrected_stop


#################################
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
            #ns_nt[nt] +=1
    return crick

######## Might only be the stop which is the first constorfs end.
def prev_con_StORF_CHECKER(prev_con_StORF,sequence,options):
    last_stop_str = int(prev_con_StORF.split(',')[-1])
    seq = sequence[last_stop_str:last_stop_str + 3]
    if not any(seq in s for s in options.stop_codons.split(',')):
        last_stop = int(last_stop_str) - 3
        last_stop = str(last_stop)
        prev_con_StORF = prev_con_StORF.replace(str(last_stop_str), last_stop)
    return prev_con_StORF

#########
def cut_seq(wc_seq,end):
    while len(wc_seq) % 3 != 0:
        if '+' in end:
            wc_seq = wc_seq[:-1] # keep removing char
        elif '-' in end:
            wc_seq = wc_seq[1:]  # keep removing char
    return wc_seq

def start_filtering(storfs):
    keep_storfs = {}
    for pos, data in storfs.items():
        starts = []
        for start_codon in ['ATG', 'GTG', 'TTG']:
            starts += [match.start() for match in re.finditer(re.escape(start_codon), data[0])]
        #### This must be slow?
        for start in starts:
            if start % 3 == 0: # maybe add a clause to report only those storfs with starts in first n of storfs
                keep_storfs.update({pos:data})
                break
    return keep_storfs

#@profile
def tile_filtering(storfs,options): #both-strand filtering
    ################ - Order largest first filtering
    storfs = sorted(storfs.items(), key=lambda storfs: storfs[1][3], reverse=True)
    ordered_by_length = OrderedDict()

    for tup in storfs:
        ordered_by_length.update({tup[0]:tup[1]})

    ############## - For each StORF_Reporter, remove all smaller overlapping STORFs according to filtering rules
    # Convert to a list for easier index-based manipulation
    ordered_by_length = list(ordered_by_length.items())
    num_storfs = len(ordered_by_length)
    i = 0

    # Step 2: Filter out nested and overlapping StORFs
    while i < num_storfs:
        pos_x, data_x = ordered_by_length[i]
        start_x = int(pos_x.split(',')[0])
        stop_x = int(pos_x.split(',')[-1])

        j = i + 1
        while j < num_storfs:
            pos_y, data_y = ordered_by_length[j]
            start_y = int(pos_y.split(',')[0])
            stop_y = int(pos_y.split(',')[-1])

            # Check if StORF Y overlaps with or is contained within StORF X
            if start_y >= stop_x or stop_y <= start_x:
                # No overlap, continue checking
                j += 1
            elif start_y >= start_x and stop_y <= stop_x:
                # StORF Y is completely contained within StORF X, remove Y
                ordered_by_length.pop(j)
                num_storfs -= 1
            else:
                # Check for partial overlap
                overlap_start = max(start_x, start_y)
                overlap_end = min(stop_x, stop_y)
                overlap = max(0, overlap_end - overlap_start + 1)  # Calculate the overlap length

                if overlap >= options.overlap_nt:
                    # Overlap exceeds the threshold, remove the shorter StORF Y
                    ordered_by_length.pop(j)
                    num_storfs -= 1
                else:
                    # Overlap is within the allowed limit, keep both
                    j += 1

        # Move to the next StORF
        i += 1

    # Step 3: Reorganize the remaining StORFs after filtering
    filtered_storfs = OrderedDict(ordered_by_length)

    # Final ordering of filtered StORFs
    if options.storf_order == 'start_pos':
        final_filtered_storfs = sortORFs(filtered_storfs)  # Sort by start position
    elif options.storf_order == 'strand':
        final_filtered_storfs = OrderedDict()
        storf_nums = sorted([item[-1] for item in filtered_storfs.values()])  # Sort by StORF number
        for num in storf_nums:
            for key, value in filtered_storfs.items():
                if value[-1] == num:
                    final_filtered_storfs.update({key: value})
                    break

    return final_filtered_storfs


def prepare_out(options, storfs, seq_id):
    gff_entries = []
    fasta_entries = {}
    for pos, data in storfs.items():
        sequence = data[0]
        strand = data[2]
        idx = data[5]
        pos_ = pos.split(',')
        start_stop = sequence[0:3]
        if len(pos_) == 3:
            if strand == '-':
                mid = int(pos_[2]) - int(pos_[1])
                mid_stop = sequence[mid:mid + 3]
            else:
                mid = int(pos_[1]) - int(pos_[0])
                mid_stop = sequence[mid:mid + 3]
        else:
            mid_stop = 'N/A'
        end_stop = sequence[-3:]
        start = int(pos_[0])
        stop = int(pos_[-1])
        ur_frame = int(data[1])
        storf_Type = data[4]
        native_seq = seq_id.replace('>','')
        ur_name = seq_id.replace('|', ':')
        length = len(sequence)
        if options.unannotated == True:
            if strand == '+':
                gff_start = str(start + 1 + int(ur_name.split('_')[-2])) # + 1 to adjust the first stop codon loci
                gff_stop = str(stop + int(ur_name.split('_')[-2]))
                if options.stop_inclusive == False:  # To remove the start and stop codon positions.
                    gff_start = str(int(gff_start) + 3)
                frame = (int(gff_stop) % 3) + 1
            elif strand == '-':
                gff_start = str(start - 2 + int(ur_name.split('_')[-2])) # -2 / -3 to adjust the first stop codon loci
                gff_stop = str(stop - 3 + int(ur_name.split('_')[-2]))
                if options.stop_inclusive == False:  # To remove the start and stop codon positions.
                    gff_stop = str(int(gff_stop) - 3)
                frame = (int(gff_stop) % 3) + 4

            storf_name = native_seq + '_' + storf_Type + '_' + str(idx) + ':' + gff_start + '-' + gff_stop

            gff_entries.append(native_seq.split('_UR')[0] + '\tSingle_Genome\t' + options.feature_type + '\t' + gff_start + '\t' + gff_stop + '\t.\t' + data[2] +
                '\t.\tID=' + storf_name + ';UR=' + ur_name.replace('>','')  + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' + str(
                    length) + ';Strand=' + data[2] +
                ';Frame=' + str(frame) + ';UR_Frame=' + str(ur_frame) +
                ';Start_Stop=' + start_stop + ';Mid_Stop=' + mid_stop  + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n')

            ### need to add - if con-storf then add middle stops
            fasta_entries.update({'>' + storf_name + ';UR=' + ur_name.replace('>','')  + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' +
                                str(length) + ';Strand=' + data[2] + ';Frame=' + str(frame) + ';UR_Frame=' + str(ur_frame) +
                ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n':data[0]})

#########################################################################
        elif options.unannotated == False: # Not done yet
            if strand == '+':
                frame = (int(stop) % 3) + 4
            elif strand == '-':
                frame = (int(stop) % 3) + 1
            storf_name = native_seq + '_' + storf_Type + '_' + str(idx) + ':' + str(start) + '-' + str(stop)

            gff_entries.append(native_seq + '\tStORF_Reporter\t' + options.feature_type + '\t' + str(start) + '\t' + str(stop) + '\t.\t' + data[2] +
                '\t.\tID=' + storf_name + ';=' + native_seq + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' + str(length) +
                               ';Frame=' + str(frame) + ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n')

            fasta_entries.update({'>' + native_seq + ';Length=' + str(length) + ';Strand=' + data[2] + ';Frame=' + str(frame) +
                ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n':data[0]})

    return gff_entries, fasta_entries


def write_gff(gff_entries,gff_out):
    ###GFF Out
    for entry in gff_entries:
            gff_out.write(entry)

def write_fasta(options, fasta_entries, fasta_out,aa_fasta_out):
    ###FASTA Prepare
    storf_num = 0 # This requires a much more elegant solution.
    ###FASTA Out
    for fasta_id, sequence in fasta_entries.items(): # could be made more efficient
        if options.stop_inclusive == False: # Remove first stop codon.
            sequence = sequence[3:]
        if options.aa_only == False:# and options.translate == False:
            fasta_out.write(fasta_id)
            if options.line_wrap:
                wrapped = textwrap.wrap(sequence, width=60)
                for wrap in wrapped:
                    fasta_out.write(wrap + '\n')
            else:
                fasta_out.write(sequence + '\n')
        if options.translate == True or options.aa_only == True:
            aa_fasta_out.write(fasta_id)
            amino = translate_frame(sequence[0:])
            if options.stop_ident == False:
                amino = amino.replace('*', '') # Remove * from sequences
            if options.line_wrap:
                amino = textwrap.wrap(amino, width=60)
                for wrap in amino:
                    aa_fasta_out.write(wrap + '\n')
            else:
                aa_fasta_out.write(amino + '\n')
        storf_num += 1

#@profile
def find_storfs(working_frame,sequence_id,stops,sequence,storfs,short_storfs,con_StORFs,frames_covered,counter,lengths,strand,StORF_idx,short_StORF_idx,Con_StORF_idx,options):
    first = True
    con_StORF_tracker = ''
    next_stops = []
    start_stops = []
    seen_stops = []
    for stop in stops:  # Finds Stop-Stop#
        seen_stops.append(stop)
        if strand == '+':
            frame = (stop % 3) + 1
        elif strand == '-':
            frame = (stop % 3) + 4
        frames_covered.update({frame: 1})
        for next_stop in stops[counter + 1:]:
            length = abs((next_stop + 3) - stop)
            if length % 3 == 0 and next_stop not in seen_stops: # In frame and not already seen?
                    if length >= options.min_orf and length <= options.max_orf:
                        if not first and stop == prev_next_stop:
                            if prev_next_stop != con_StORF_tracker:
                                con_StORF_tracker = next_stop
                                seq = sequence[prev_stop:next_stop + 3]
                                length = next_stop - prev_stop
                                ##### Needed to correct for negative frame loci
                                if working_frame == 'negative':
                                    rev_corrected_start, rev_corrected_mid, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,prev_stop,stop,next_stop  + 3)
                                    con_StORF_Pos = ",".join([str(rev_corrected_start), str(rev_corrected_mid), str(rev_corrected_stop)])
                                else:
                                    con_StORF_Pos = ",".join([str(prev_stop), str(stop), str(next_stop  + 3)])
                                #####
                                con_length = abs((next_stop + 3) - prev_stop)
                                con_StORFs.update({con_StORF_Pos: [seq, str(frame), strand, con_length,'Con-StORF',Con_StORF_idx]})
                                Con_StORF_idx +=1
                            elif stop == con_StORF_tracker:
                                con_StORF_tracker = next_stop
                                prev_con_StORF = next(reversed(con_StORFs.keys())) #Get last key
                                seq_start = int(prev_con_StORF.split(',')[0])
                                #seq = sequence[seq_start:next_stop  + 3]
                                seq = sequence[int(prev_con_StORF.split(',')[0]):next_stop + 3]
                                length = next_stop - seq_start # Check
                                prev_con_StORF = prev_con_StORF_CHECKER(prev_con_StORF,sequence,options)
                                ##### Needed to correct for negative frame loci
                                if working_frame == 'negative':
                                    rev_corrected_start, rev_corrected_mid, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,prev_stop,stop,next_stop  + 3)
                                    con_StORF_Pos = ",".join([str(rev_corrected_start), str(rev_corrected_mid), str(rev_corrected_stop)])
                                else:
                                    con_StORF_Pos = ",".join([str(prev_stop), str(stop), str(next_stop  + 3)])
                                #####
                                if options.olap_filtering == 'both-strand':
                                    con_StORFs.popitem()
                                con_length = abs((next_stop + 3) - prev_stop)
                                con_StORFs.update({con_StORF_Pos: [seq, str(frame), strand, con_length,'Con-StORF',Con_StORF_idx]})
                                Con_StORF_idx +=1

                        if options.olap_filtering == 'none': # This could be made more efficient
                            seq = sequence[stop:next_stop  + 3]
                            ##### Needed to correct for negative frame loci
                            if working_frame == 'negative':
                                rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                storfs.update({",".join([str(rev_corrected_start), str(rev_corrected_stop)]): [seq, str(frame), strand, length, 'StORF', StORF_idx]})
                            else:
                                storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame), strand, length,'StORF',StORF_idx]})
                            #####
                            StORF_idx +=1
                            lengths.append(length)
                        elif not first:
                            if stop > prev_stop and next_stop < prev_next_stop: # Check
                                break
                            storf = set(range(stop, next_stop + 4)) # + 4 to account for set use
                            storf_overlap = len(prev_storf.intersection(storf))
                            if storf_overlap <= options.overlap_nt:
                                seq = sequence[stop:next_stop  + 3]
                                ##### Needed to correct for negative frame loci
                                if working_frame == 'negative':
                                    rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                    storfs.update({",".join([str(rev_corrected_start), str(rev_corrected_stop)]): [
                                        seq, str(frame), strand, length, 'StORF', StORF_idx]})
                                else:
                                    storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame), strand,
                                                                                               length, 'StORF',
                                                                                               StORF_idx]})
                                #####
                                StORF_idx +=1
                                seen_stops.append(next_stop)#  + 3)
                                prev_storf = storf
                                prev_stop = stop
                                prev_next_stop = next_stop # Check
                                next_stops.append(next_stop)
                                start_stops.append(stop)
                                prevlength = prev_next_stop - prev_stop
                            elif storf_overlap >= options.overlap_nt and options.olap_filtering == 'both-strand':  # and length > prevlength:
                                if length > prevlength:
                                    storfs.popitem()
                                    seq = sequence[stop:next_stop  + 3]
                                    ##### Needed to correct for negative frame loci
                                    if working_frame == 'negative':
                                        rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                        storfs.update({",".join(
                                            [str(rev_corrected_start), str(rev_corrected_stop)]): [seq, str(frame),
                                                                                                       strand, length,
                                                                                                       'StORF',
                                                                                                       StORF_idx]})
                                    else:
                                        storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame),
                                                                                                   strand, length,
                                                                                                   'StORF', StORF_idx]})
                                    #####
                                    StORF_idx +=1
                                    seen_stops.append(next_stop)
                                    prev_storf = storf
                                    prev_stop = stop
                                    prev_next_stop = next_stop
                                    next_stops.append(next_stop)
                                    start_stops.append(stop)
                                    prevlength = prev_next_stop - prev_stop
                            else: # If filtering is none or single-strand, we do not remove overlapping StORFs on the same strand
                                seq = sequence[stop:next_stop  + 3]
                                ##### Needed to correct for negative frame loci
                                if working_frame == 'negative':
                                    rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                    storfs.update({",".join([str(rev_corrected_start), str(rev_corrected_stop)]): [
                                        seq, str(frame), strand, length, 'StORF', StORF_idx]})
                                else:
                                    storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame), strand,
                                                                                               length, 'StORF',
                                                                                               StORF_idx]})
                                #####
                                StORF_idx += 1
                                seen_stops.append(next_stop)
                                prev_storf = storf
                                prev_stop = stop
                                prev_next_stop = next_stop
                                next_stops.append(next_stop)
                                start_stops.append(stop)
                                prevlength = prev_next_stop - prev_stop
                        elif first:
                            if options.partial_storf: # upstream partial StORF_Reporter
                                if stop > options.min_orf and frames_covered[frame] != 1:
                                    seq = sequence[0:stop + 3] #Start of seq to first stop identified
                                    #ps_seq = cut_seq(seq, '-')
                                    ##### Needed to correct for negative frame loci
                                    if working_frame == 'negative':
                                        rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                        storfs.update({",".join(
                                            [str(rev_corrected_start), str(rev_corrected_stop)]): [seq, str(frame),
                                                                                                       strand, length,
                                                                                                       'StORF',
                                                                                                       StORF_idx]})
                                    else:
                                        storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame),
                                                                                                   strand, length,
                                                                                                   'StORF', StORF_idx]})
                                    #####
                                    StORF_idx +=1
                            seq = sequence[stop:next_stop  + 3]
                            length = next_stop - stop
                            ##### Needed to correct for negative frame loci
                            if working_frame == 'negative':
                                rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(options,len(sequence),sequence_id,stop,None,next_stop  + 3)
                                storfs.update({",".join([str(rev_corrected_start), str(rev_corrected_stop)]): [seq, str(frame), strand, length, 'StORF', StORF_idx]})
                            else:
                                storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame), strand, length,'StORF',StORF_idx]})
                            #####
                            StORF_idx +=1
                            seen_stops.append(next_stop)
                            prev_storf = set(range(stop, next_stop + 4)) # + 4 to account for set use
                            prev_stop = stop
                            prev_next_stop = next_stop
                            next_stops.append(next_stop)
                            start_stops.append(stop)
                            prevlength = prev_next_stop - prev_stop
                            first = False
                        break

                    if options.short_storfs != False and length >= 30: # Report short (<= 120) StORFs
                        seq = sequence[stop:next_stop  + 3]
                        length = next_stop - stop # Check
                        short_storfs.update({",".join([str(stop), str(next_stop  + 3)]): [seq, str(frame), strand, length,'Short-StORF', short_StORF_idx]})
                        short_StORF_idx +=1
                    else:
                        break
        counter +=1
    if options.partial_storf:  # downstream partial StORF_Reporter - Last Stop to end of sequence
        try:
            if (len(sequence) - stop) > options.min_orf:
                seq = sequence[stop:len(sequence)]  # Start of seq to first stop identified - not working
                ps_seq = cut_seq(seq, '+')
                storfs.update({",".join([str(stop), str(len(sequence))]): [ps_seq, str(frame), strand, stop,'Partial-StORF',StORF_idx]})
                StORF_idx +=1
        except UnboundLocalError:
            pass

    return storfs, short_storfs, con_StORFs, frames_covered, counter, lengths, StORF_idx, Con_StORF_idx

def STORF_Finder(options, sequence_info, sequence_id, fasta_out, aa_fasta_out, gff_out, split_index): #Main Function
    ## If UR is the start of a sequence the 0/1 base position throws off the start of the StORF
    if sequence_id.split('_')[split_index] == '1':
        start_of_seq = True
    else:
        start_of_seq = False


    sequence = sequence_info[1]
    stops = []
    frames_covered = OrderedDict()
    for x in range (1,7):
        frames_covered.update({x: 0})
    for stop_codon in options.stop_codons.split(','): #Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence)]
    stops.sort()
    storfs = OrderedDict()
    short_storfs = OrderedDict()
    con_StORFs = OrderedDict()
    counter = 0
    lengths = []
    StORF_idx = 0
    short_StORF_idx = 0
    Con_StORF_idx = 0
    storfs, short_storfs, con_StORFs,frames_covered,counter,lengths,StORF_idx,Con_StORF_idx = find_storfs("positive",sequence_id,stops,sequence,storfs,short_storfs,con_StORFs,frames_covered,counter,lengths,'+',StORF_idx,short_StORF_idx,Con_StORF_idx,options)
    ###### Reversed Comppliment
    sequence_rev = revCompIterative(sequence)
    stops = []
    for stop_codon in options.stop_codons.split(','): #Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence_rev)]
    stops.sort()
    counter = 0
    storfs, short_storfs, con_StORFs,frames_covered,counter,lengths,StORF_idx,Con_StORF_idx = find_storfs("negative",sequence_id,stops,sequence_rev,storfs,short_storfs,con_StORFs,frames_covered,counter,lengths,'-',StORF_idx,short_StORF_idx,Con_StORF_idx,options)

    ## The correction for base 0/1 position in the UR
    if start_of_seq == True:
        for key in list(storfs.keys()):
            if storfs[key][2] == '-':  # Check if the strand is '-'
                positions = key.split(',')
                new_positions = [str(int(pos) - 1) for pos in positions]
                new_key = ','.join(new_positions)
                storfs[new_key] = storfs.pop(key)

    #Potential run-through StORFs
    if options.whole_contig:
        for frame,present in frames_covered.items():
            if present == 0:
                if frame <4:
                    wc_seq = sequence[frame-1:]
                    wc_seq = cut_seq(wc_seq,'+')
                    storfs.update({",".join([str(0), str(len(sequence))]): [wc_seq, str(frame), '+', len(sequence),'Run-Through-StORF',StORF_idx]})
                    StORF_idx +=1
                else:
                    wc_seq = sequence_rev[frame-4:]
                    wc_seq = cut_seq(wc_seq,'+')
                    storfs.update({",".join([str(0), str(len(sequence))]): [wc_seq, str(frame), '-', len(sequence_rev),'Run-Through-StORF',StORF_idx]})
                    StORF_idx +=1

####################################### Writing output
    ######## Only StORFs
    #Check if there are StORFs to report
    if options.con_storfs == False and options.con_only == False and options.short_storfs == False:
        if bool(storfs):
            if options.start_filtering == True:
                storfs = start_filtering(storfs)
            if options.olap_filtering == 'both-strand':
                storfs = tile_filtering(storfs, options)  # Filtering
            storfs = OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))  # Reorder by start position
            if options.reporter == True:
                return storfs
            ###Data Prepare
            gff_entries, fasta_entries = prepare_out(options, storfs, sequence_id)
            write_fasta(options, fasta_entries, fasta_out, aa_fasta_out)
            if not options.aa_only:
                write_gff(gff_entries, gff_out)

    ###### Only Short-StORFs
    # elif options.short_storfs != False and options.short_storfs_only == True: # Short-StORFs ONLY
    #     if bool(short_storfs):
    #         if options.start_filtering == True:
    #             short_storfs = start_filtering(short_storfs)
    #         if options.olap_filtering == 'both-strand':
    #             short_storfs = tile_filtering(short_storfs, options)  # Filtering
    #         short_storfs = OrderedDict(sorted(short_storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))  # Reorder by start position
    #         ###Data Prepare
    #         gff_entries, fasta_entries = prepare_out(options, short_storfs, sequence_id)
    #         write_fasta(options, fasta_entries, fasta_out, aa_fasta_out)
    #         if not options.aa_only:
    #             write_gff(gff_entries, gff_out)

    ### Short-StORFs
    elif options.short_storfs != False: # and options.short_storfs_only == False:
        ### Don't allow short-storfs to overlap with storfs
        if options.short_storfs == 'Nolap':
            all_StORFs = {**storfs, **short_storfs}
            if options.olap_filtering == 'both-strand':
                all_StORFs = tile_filtering(all_StORFs,options) # Filtering
            filtered_StORFs = OrderedDict(sorted(all_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(","))))) # Reorder by start position

        ### short-storfs can onverlap with storfs
        elif options.short_storfs == 'Olap':
            if options.olap_filtering == 'both-strand': # Filter individually
                storfs = tile_filtering(storfs,options) # Filtering
                short_storfs = tile_filtering(short_storfs,options) # Filtering
            all_StORFs = {**storfs, **short_storfs}
            filtered_StORFs = OrderedDict(sorted(all_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
        if options.short_storfs_only == True: # Checking what short_storfs survived filtering and extracting them
            final_StORFs = dict(short_storfs.items() & filtered_StORFs.items())
            if len(final_StORFs) != 0:
                print("we have one: " + options.fasta)
        else:
            final_StORFs = filtered_StORFs

        if options.reporter == True:
            return final_StORFs
            ###Data Prepare
        gff_entries, fasta_entries = prepare_out(options, final_StORFs, sequence_id)
        write_fasta(options, fasta_entries, fasta_out, aa_fasta_out)
        if not options.aa_only:
            write_gff(gff_entries, gff_out)

    ####### StORFs and Con-StORFs
    elif options.con_storfs == True and options.con_only == False:
        all_StORFs = {**storfs, **con_StORFs}
        if bool(storfs):
            if options.olap_filtering == 'both-strand':
                all_StORFs = tile_filtering(all_StORFs, options) # Filtering
            all_StORFs = OrderedDict(sorted(all_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(","))))) # Reorder by start position
            if options.reporter == True:
                return all_StORFs
            ###Data Prepare
            gff_entries, fasta_entries = prepare_out(options, all_StORFs, sequence_id)
            write_fasta(options, fasta_entries, fasta_out, aa_fasta_out)
            if not options.aa_only:
                write_gff(gff_entries, gff_out)

    ###### Con-StORFs only
    elif options.con_only == True:
        if options.olap_filtering == 'both-strand':
            con_StORFs = tile_filtering(con_StORFs, options)
            con_StORFs = OrderedDict(sorted(con_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(","))))) # Reorder by start position
        else:
            con_StORFs = OrderedDict(sorted(con_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(","))))) # Reorder by start position
        if options.reporter == True:
            return con_StORFs
        ###Data Prepare
        gff_entries, fasta_entries = prepare_out(options, con_StORFs, sequence_id)
        write_fasta(options, fasta_entries, fasta_out, aa_fasta_out)
        if not options.aa_only:
            write_gff(gff_entries, gff_out)
    ###### Below won't work..?
    elif options.verbose == True:
        print("No StOFS Found")

def fasta_load(fasta_in, sequence_regions, sequences):
    first = True
    sequence_region_length = 0
    for line in fasta_in:
        line = line.strip()
        if '##sequence-region' in line:
            sequence_region_length = int(line.split(' ')[-1]) # bug and wont work on non-UR runs
            sequence_regions.append(line)
        elif line.startswith((';','\n','#')):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: [sequence_region_length,seq.strip()]})
            seq = ''
            sequence_name = line.strip()
        elif line.startswith('>'):
            seq = ''
            sequence_name = line.strip()
        elif line:
            seq += str(line)
            first = False
    sequences.update({sequence_name: [sequence_region_length,seq.strip()]})
    return sequence_regions, sequences

## Function to control how StORF-Finder handles Single_Genome output
def StORF_Reported(options, Contigs):
    options.unannotated = True
    Reporter_StORFs = collections.OrderedDict()
    for Contig_ID, Contig_URs in Contigs.items():
        Reporter_StORFs.update({Contig_ID:[]})
        URs = Contig_URs[3]
        try:
            for UR in URs:
                if len(URs[UR][1]) >= options.min_orf:
                    contig_length = Contigs[Contig_ID][1]
                    sequence_info = [contig_length,URs[UR][1]]
                    if UR.split('_')[0] == '0': # This is to account for the GFF base-1 system
                        sequence_id = "1_" + UR.split('_')[1]
                    else:
                        sequence_id = UR # mockup sequence_id in correct format for later
                    StORFs = STORF_Finder(options, sequence_info, sequence_id, None, None, None,0) # need to pass seq and original seq length
                    if StORFs: #  Left out for now to allow for tracking of non-StORF URs
                        for StORF in StORFs.values():
                            StORF.append(URs[UR][0]) # True UR
                            StORF.append(sequence_id) # Extended UR
                        Reporter_StORFs[Contig_ID].append(StORFs)
        except TypeError:
        #    if options.verbose == True:
            print("No URs in seq")
    return Reporter_StORFs


def main():
    parser = argparse.ArgumentParser(description='StORF-Reporter ' + StORF_Reporter_Version + ': StORF-Finder Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-f', action="store", dest='fasta', required=True,
                        help='Input FASTA File - (UR_Extractor output)')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-reporter', action="store", dest='reporter', default=False, required=False,
                        help=argparse.SUPPRESS)
    optional.add_argument('-ua', dest='unannotated', action='store', default=True, type=eval, choices=[True, False],
                        help='Default - Treat input as Unannotated: Use "-ua False" for standard fasta')
    optional.add_argument('-wc', action="store", dest='whole_contig', default=False, type=eval, choices=[True, False],
                        help='Default - False: StORFs reported across entire sequence')
    optional.add_argument('-ps', action="store", dest='partial_storf', default=False, type=eval, choices=[True, False],
                        help='Default - False: Partial StORFs reported')
    optional.add_argument('-olap_filt', action='store', dest='olap_filtering', default='both-strand', const='both-strand', nargs='?',
                        choices=['none', 'single-strand', 'both-strand'],
                        help='Default - "both-strand": Filtering level "none" is not recommended, "single-strand" for single strand filtering '
                             'and both-strand for both-strand longest-first tiling')
    optional.add_argument('-start_filt', action="store", dest='start_filtering', default=False, type=eval, choices=[True, False],
                        help='Default - False: Filter out StORFs without at least one of the 3 common start codons (best used for short-storfs).')
    optional.add_argument('-con_storfs', action="store", dest='con_storfs', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output Consecutive StORFs')
    optional.add_argument('-con_only', action="store", dest='con_only', default=False, type=eval, choices=[True, False],
                        help='Default - False: Only output Consecutive StORFs')
    optional.add_argument('-short_storfs', action="store", dest='short_storfs', default=False, type=str, choices=[False, 'Nolap', 'Olap'],
                        help='Default - False: Run StORF-Finder in "Short-StORF" mode. Will only return StORFs between 30 and 120 nt '
                             'that do not overlap longer StORFs - Only works with StORFs for now. "Nolap" will filter Short-StORFs which are'
                             'overlapped by StORFs and Olap will report Short-StORFs which do overlap StORFs. Overlap is defined by "-olap".')
    optional.add_argument('-short_storfs_only', action="store", dest='short_storfs_only', default=False, type=eval, choices=[True, False],
                        help='Default - True. Only report Short-StORFs?')
    optional.add_argument('-f_type', action='store', dest='feature_type', default='CDS', const='StORF', nargs='?',
                        choices=['StORF', 'CDS', 'ORF'],
                        help='Default - "StORF": Which GFF feature type for StORFs to be reported as in GFF')
    optional.add_argument('-minorf', action="store", dest='min_orf', default=99, type=int,
                        help='Default - 99: Minimum StORF size in nt')
    optional.add_argument('-maxorf', action="store", dest='max_orf', default=60000, type=int,
                        help='Default - 60kb: Maximum StORF size in nt')
    optional.add_argument('-codons', action="store", dest='stop_codons', default="TAG,TGA,TAA",
                        help='Default - (\'TAG,TGA,TAA\'): List Stop Codons to use')
    optional.add_argument('-non_standard', action="store", dest='non_standard', default="0.20",
                          help='Default - 0.20: Reject StORFs with >=20%% non-standard nucleotides (A,T,G,C) - Provide %% as decimal')
    optional.add_argument('-olap', action="store", dest='overlap_nt', default=50, type=int,
                        help='Default - 50: Maximum number of nt of a StORF which can overlap another StORF.')
    optional.add_argument('-s', action="store", dest='suffix', required=False,
                        help='Default - Do not append suffix to genome ID')
    optional.add_argument('-so', action="store", dest='storf_order', default='start_pos', nargs='?', choices=['start_pos','strand'],
                        required=False,
                        help='Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.')



    output = parser.add_argument_group('Output')
    output.add_argument('-oname', action="store", dest='o_name', required=False,
                        help='Default - Appends \'_StORF-Finder\' to end of input FASTA filename')
    output.add_argument('-odir', action="store", dest='o_dir', required=False,
                        help='Default -  Same directory as input FASTA')
    # output.add_argument('-af', action="store", dest='affix', required=False,
    #                     help='Default - \'StORF-R\': can be used to append an identifier to output filename '
    #                          'to distinguish Con-StORF from StORF runs)')
    output.add_argument('-gff', action='store', dest='gff', default=True, type=eval, choices=[True, False],
                        help='Default - True: Output a GFF file')
    output.add_argument('-aa', action="store", dest='translate', default=False, type=eval, choices=[True, False],
                        help='Default - False: Report StORFs as amino acid sequences')
    output.add_argument('-aa_only', action="store", dest='aa_only', default=False, type=eval, choices=[True, False],
                        help='Default - False: Only output Amino Acid Fasta')
    output.add_argument('-lw', action="store", dest='line_wrap', default=True, type=eval, choices=[True, False],
                        help='Default - True: Line wrap FASTA sequence output at 60 chars')
    output.add_argument('-spos', action="store", dest='stop_inclusive', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output StORF sequences and GFF positions inclusive of first stop codon -'
                             'This can break some downstream tools if changed to True.')
    output.add_argument('-stop_ident', action="store", dest='stop_ident', default=False, choices=[True, False],
                        help='Default - True: Identify Stop Codon positions with \'*\'')
    output.add_argument('-gff_fasta', action="store", dest='gff_fasta', default=False, type=eval, choices=[True, False],
                        help='Default - False: Report all gene sequences (nt) at the bottom of GFF files in Prokka output mode')
    output.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    optional.add_argument('-nout', action='store', dest='nout', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)

    misc = parser.add_argument_group('Misc')
    misc.add_argument('-verbose', action='store', dest='verbose', default=False, type=eval, choices=[True, False],
                        help='Default - False: Print out runtime messages')
    misc.add_argument('-v', action='store_true', dest='version',
                        help='Default - False: Print out version number and exit')

    options = parser.parse_args()

    if options.fasta == None:
        if options.version:
            sys.exit(StORF_Reporter_Version)
        else:
            exit('StORF-Finder: error: the following arguments are required: -f')

    print("Thank you for using StORF-Finder -- A detailed user manual can be found at https://github.com/NickJD/StORF-Reporter\n"
          "Please report any issues to: https://github.com/NickJD/StORF-Reporter/issues\n#####")

    #############
    options.fasta = os.path.normpath(options.fasta)
    options.fasta = os.path.realpath(options.fasta)
    #############


    #ns_nt = defaultdict  # Used to Record non-standard nucleotides - not implemented yet
    ##### Load in fasta file
    sequences = OrderedDict()
    sequence_regions = []
    if options.reporter == False:
        try: # Detect whether fasta files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta,'rt')
            sequence_regions, sequences = fasta_load(fasta_in, sequence_regions, sequences)
        except:
            fasta_in = open(options.fasta,'r')
            sequence_regions, sequences =  fasta_load(fasta_in, sequence_regions, sequences)
        if options.verbose == True:
            print(fasta_in.name)
    #### Output Directory and Filename handling
    if options.o_dir == None and options.o_name == None:
        tmp_extension = options.fasta.split('.')[-1]  # could be .fa/.fasta etc
        output_file = options.fasta.replace('.' + tmp_extension, '')
        output_file = output_file + '_StORF-Finder'
    elif options.o_dir != None and options.o_name != None:
        output_file = options.o_dir
        output_file = output_file + options.o_name
    elif options.o_dir != None:
        output_file = options.fasta.split(os.sep)[-1].split('.')[0]
        output_file = options.o_dir + output_file + '_StORF-Finder'
    elif options.o_name != None:
        tmp_filename = options.fasta.split(os.sep)[-1]  # could be .fa/.fasta etc
        output_file = options.fasta.replace(tmp_filename, '')
        output_file = output_file + options.o_name

    if not options.gz: # Clear fasta and gff files if not empty - Needs an elegant solution
        if not options.aa_only:
            gff_out = open(output_file + '.gff', 'w', newline='\n', encoding='utf-8')
            gff_out.write("##gff-version\t3\n#\tSingle_Genome - Stop ORF Predictions\n#\tRun Date:" + str(date.today()) + '\n')
            gff_out.write('##Single_Genome ' + StORF_Reporter_Version + '\n')
            for seq_reg in sequence_regions:
                gff_out.write(seq_reg + '\n')
            gff_out.write("##Original File: " + options.fasta.split(os.sep)[-1] + '\n\n')
            fasta_out = open(output_file + '.fasta', 'w', newline='\n', encoding='utf-8')
            if options.translate:
                aa_fasta_out = open(output_file + '_aa.fasta', 'w', newline='\n', encoding='utf-8')
            else:
                aa_fasta_out = None
        elif options.aa_only:
            aa_fasta_out = open(output_file + '_aa.fasta', 'w', newline='\n', encoding='utf-8')
    elif options.gz:
        if not options.aa_only:
            gff_out = gzip.open(output_file + '.gff.gz', 'wt', newline='\n', encoding='utf-8')
            gff_out.write("##gff-version\t3\n#\tSingle_Genome - Stop ORF Predictions\n#\tRun Date:" + str(date.today()) + '\n')
            gff_out.write('##Single_Genome ' + StORF_Reporter_Version + '\n')
            for seq_reg in sequence_regions:
                gff_out.write(seq_reg + '\n')
            gff_out.write("##Original File: " + options.gff + '\n\n')
            fasta_out = gzip.open(output_file + '.fasta.gz', 'wt', newline='\n', encoding='utf-8')
            if options.translate:
                aa_fasta_out = gzip.open(output_file + '_aa.fasta.gz', 'wt', newline='\n', encoding='utf-8')
            else:
                aa_fasta_out = None
        elif options.aa_only:
            aa_fasta_out = gzip.open(output_file + '_aa.fasta.gz', 'wt', newline='\n', encoding='utf-8')

    for sequence_id, sequence_info in sequences.items():
        if len(sequence_info[1]) >= options.min_orf:
            STORF_Finder(options, sequence_info, sequence_id, fasta_out, aa_fasta_out, gff_out,3)

if __name__ == "__main__":
    main()
    print("Complete")
