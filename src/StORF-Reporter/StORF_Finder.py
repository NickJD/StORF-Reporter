import argparse
import collections
import re
from collections import defaultdict, OrderedDict
from datetime import date
import textwrap
import gzip

try:
    from ORForise.utils import sortORFs  # Calling from ORForise via pip
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from ORForise.src.ORForise.utils import sortORFs  # Calling from ORForise locally (StORF-Reporter and ORForise in same dir)



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
############################

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
            ns_nt[nt] +=1
    return crick

######## REWRTIE THIS!! - Might only be the stop which is the first constorfs end.
def prev_con_StORF_CHECKER(prev_con_StORF,sequence,options):
    last_stop_str = int(prev_con_StORF.split(',')[-1])
    seq = sequence[last_stop_str:last_stop_str + 3]
    if not any(seq in s for s in options.stop_codons.split(',')):
        last_stop = int(last_stop_str) - 3
        last_stop = str(last_stop)
        prev_con_StORF = prev_con_StORF.replace(str(last_stop_str), last_stop)
        seq = sequence[int(last_stop):int(last_stop) + 3]

    return prev_con_StORF

#########
def cut_seq(wc_seq,end):
    while len(wc_seq) % 3 != 0:
        if '+' in end:
            wc_seq = wc_seq[:-1] # keep removing char
        elif '-' in end:
            wc_seq = wc_seq[1:]  # keep removing char
    return wc_seq

def tile_filtering(storfs,options): #Hard filtering
    storfs = OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
    ################ - Order largest first filtering
    storfs = sorted(storfs.items(), key=lambda storfs:storfs[1][3],reverse=True)
    ordered_by_length = OrderedDict()
    for tup in storfs:
        ordered_by_length.update({tup[0]:tup[1]}) # -1 for last?
    ############## - For each StORF_Reporter, remove all smaller overlapping STORFs according to filtering rules
    num_storfs = len(ordered_by_length)
    i = 0
    ordered_by_length = list(ordered_by_length.items())
    while i < num_storfs:
        pos_x, data_x = ordered_by_length[i]
        start_x = int(pos_x.split(',')[0])
        stop_x = int(pos_x.split(',')[-1])
        j = i+1
        while j < num_storfs:
            pos_y, data_y = ordered_by_length[j]
            start_y = int(pos_y.split(',')[0])
            stop_y = int(pos_y.split(',')[-1])
            if start_y >= stop_x or stop_y <= start_x:
                j+=1
                continue  # Not caught up yet / too far
            elif start_y >= start_x and stop_y <= stop_x:
                ordered_by_length.pop(j)
                num_storfs = len(ordered_by_length)
            else: # +1 needed for stop codon
                x = set(range(start_x,stop_x+1))
                y = set(range(start_y,stop_y+1))
                overlap = len(x.intersection(y))
                if overlap >= options.overlap_nt:
                    ordered_by_length.pop(j)
                    num_storfs = len(ordered_by_length)
                else:
                    j += 1
        num_storfs = len(ordered_by_length)
        i+=1
    filtered_storfs = OrderedDict(ordered_by_length)
    #### Clunky - reordering of StORFs - Through number found (pos then neg strand) or start position?
    if options.storf_order == 'start_pos':
        final_filtered_storfs = sortORFs(filtered_storfs)
    elif options.storf_order == 'strand':
        final_filtered_storfs = collections.OrderedDict()
        storf_nums = [item[-1] for item in filtered_storfs.values()]
        storf_nums = sorted(storf_nums)
        for num in storf_nums:
            for key,value in filtered_storfs.items():
                if value[-1] == num:
                    final_filtered_storfs.update({key:value})
                    break

    return final_filtered_storfs

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate

def prepare_out(options,storfs,seq_id):
    gff_entries = []
    fasta_entries = {}
    for pos, data in storfs.items():
        sequence = data[0]
        strand = data[2]
        idx = data[5]
        start_stop = sequence[0:3]
        end_stop = sequence[-3:]
        pos_ = pos.split(',')
        start = int(pos_[0])
        stop = int(pos_[-1])
        ur_frame = int(data[1])
        storf_Type = data[4]
        #if options.suffix == None:
        #    seq_id = seq_id.replace('_UR','')
        #seq_id = seq_id.split()[0].replace('>', '')
        native_seq = seq_id.replace('>','')
        ur_name = seq_id.replace('|', ':')
        length = len(sequence)
        if options.unannotated == True:
            gff_start = str(start + int(ur_name.split('_')[2]))
            gff_stop = str(stop + int(ur_name.split('_')[2]))
            if strand == '+':
                frame = (int(gff_stop) % 3) + 1
            elif strand == '-':
                frame = (int(gff_stop) % 3) + 4
            storf_name = native_seq + '_' + storf_Type + '_' + str(idx) + ':' + gff_start + '-' + gff_stop
            gff_entries.append(native_seq + '\tStORF_Reporter\t' + options.feature_type + '\t' + gff_start + '\t' + gff_stop + '\t.\t' + data[2] +
                '\t.\tID=' + storf_name + ';UR=' + ur_name.replace('>','')  + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' + str(
                    length) + ';Strand=' + data[2] +
                ';Frame=' + str(frame) + ';UR_Frame=' + str(ur_frame) +
                ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n')
            ### need to add - if con-storf then add middle stops
            fasta_entries.update({'>' + storf_name + ';UR=' + ur_name.replace('>','')  + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' +
                                str(length) + ';Strand=' + data[2] + ';Frame=' + str(frame) + ';UR_Frame=' + str(ur_frame) +
                ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n':data[0]})

#########################################################################
        elif options.intergenic == False: # Not done yet
            fa_id = (">" + str(storf_name) + "|" + str(start) + strand + str(stop) + "|Frame:" + str(
                frame) + '|Start_Stop=' + start_stop +
                     '|End_Stop=' + end_stop + '|StORF_Type:' + storf_Type + "\n")

            gff_entries.append(                native_seq + '\tStORF_Reporter\t' + options.feature_type + '\t' + gff_start + '\t' + gff_stop + '\t.\t' + data[2] +
                '\t.\tID=' + storf_name + ';UR=' + ur_name + ';UR_Stop_Locations=' + '-'.join(pos_) + ';Length=' + str(
                    length) +
                ';Frame=' + str(frame) + ';UR_Frame=' + str(ur_frame) +
                ';Start_Stop=' + start_stop + ';End_Stop=' + end_stop + ';StORF_Type=' + storf_Type + '\n')

    return gff_entries, fasta_entries


def write_gff(gff_entries,gff_out):
    ###GFF Out
    for entry in gff_entries:
            gff_out.write(entry)

def write_fasta(fasta_entries, fasta_out,aa_fasta_out):  # Some Lines commented out for BetaRun of ConStORF
    ###FASTA Prepare
    storf_num = 0 # This requires a much more elegant solution.
    ###FASTA Out
    for fasta_id, sequence in fasta_entries.items(): # could be made more efficient
        strand = fasta_id.split('Strand=')[1].split(';')
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
            if "+" in strand:
                amino = translate_frame(sequence[0:])
                if options.stop_ident == False:
                    amino = amino.replace('*', '') # Remove * from sequences
                if options.line_wrap:
                    amino = textwrap.wrap(amino, width=60)
                    for wrap in amino:
                        aa_fasta_out.write(wrap + '\n')
                else:
                    aa_fasta_out.write(amino + '\n')
            if "-" in strand:
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


def find_storfs(working_frame,stops,sequence,storfs,con_StORFs,frames_covered,counter,lengths,strand,StORF_idx,Con_StORF_idx,options):
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
            length = abs(next_stop - stop)
            if length % 3 == 0 and next_stop not in seen_stops:
                    if length >= options.min_orf and length <= options.max_orf:
                        if not first and stop == prev_next_stop:
                            if prev_next_stop != con_StORF_tracker:
                                con_StORF_tracker = next_stop
                                seq = sequence[prev_stop:next_stop + 3]
                                length = next_stop - prev_stop
                                con_StORF_Pos = ",".join([str(prev_stop), str(stop+3), str(next_stop+3)])
                                con_StORFs.update({con_StORF_Pos: [seq, str(frame), strand, length,'Con-Stop-ORF',Con_StORF_idx]})
                                Con_StORF_idx +=1
                            elif stop == con_StORF_tracker:
                                con_StORF_tracker = next_stop
                                prev_con_StORF = next(reversed(con_StORFs.keys())) #Get last key
                                seq_start = int(prev_con_StORF.split(',')[0])
                                seq = sequence[seq_start:next_stop + 3]
                                length = next_stop - seq_start
                                ############# Fix the +3 issue for last but now internal stop position
                                prev_con_StORF = prev_con_StORF_CHECKER(prev_con_StORF,sequence,options)
                                ############ Fixed
                                con_StORF_Pos = prev_con_StORF+','+str(next_stop)
                                if options.filtering == 'hard':
                                    con_StORFs.popitem()
                                con_StORFs.update({con_StORF_Pos: [seq, str(frame), strand, length,'Con-Stop-ORF',Con_StORF_idx]})
                                Con_StORF_idx +=1

                        if options.filtering == 'none': # This could be made more efficient
                            seq = sequence[stop:next_stop + 3]
                            storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length,'Stop-ORF',StORF_idx]})
                            StORF_idx +=1
                            lengths.append(length)
                        elif not first:
                            if stop > prev_stop and next_stop < prev_next_stop:
                                break
                            storf = set(range(stop, next_stop + 4))
                            storf_overlap = len(prev_storf.intersection(storf))
                            if storf_overlap <= options.overlap_nt:
                                seq = sequence[stop:next_stop + 3]
                                storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length,'Stop-ORF',StORF_idx]})
                                StORF_idx +=1
                                seen_stops.append(next_stop)
                                prev_storf = storf
                                prev_stop = stop
                                prev_next_stop = next_stop
                                next_stops.append(next_stop)
                                start_stops.append(stop)
                                prevlength = prev_next_stop - prev_stop
                            elif storf_overlap >= options.overlap_nt and options.filtering == 'hard':  # and length > prevlength:
                                if length > prevlength:
                                    storfs.popitem()
                                    seq = sequence[stop:next_stop + 3]
                                    storfs.update({",".join([str(stop), str(next_stop + 3)]): [seq, str(frame), strand,
                                                                                               length, 'Stop-ORF',
                                                                                               StORF_idx]})
                                    StORF_idx +=1
                                    seen_stops.append(next_stop)
                                    prev_storf = storf
                                    prev_stop = stop
                                    prev_next_stop = next_stop
                                    next_stops.append(next_stop)
                                    start_stops.append(stop)
                                    prevlength = prev_next_stop - prev_stop
                            else: # If filtering is none or soft, we do not remove overlapping StORFs on the same strand
                                seq = sequence[stop:next_stop + 3]
                                storfs.update({",".join([str(stop), str(next_stop + 3)]): [seq, str(frame), strand,
                                                                                           length, 'Stop-ORF',
                                                                                           StORF_idx]})
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
                                    ps_seq = cut_seq(seq, '-')
                                    storfs.update({",".join([str(0), str(stop)]): [ps_seq, str(frame), strand, stop,'Partial-StORF',StORF_idx]})
                                    StORF_idx +=1
                            seq = sequence[stop:next_stop + 3]
                            length = next_stop - stop
                            storfs.update({",".join([str(stop), str(next_stop+3)]): [seq, str(frame), strand, length,'Stop-ORF',StORF_idx]})
                            StORF_idx +=1
                            seen_stops.append(next_stop)
                            prev_storf = set(range(stop, next_stop + 4))
                            prev_stop = stop
                            prev_next_stop = next_stop
                            next_stops.append(next_stop)
                            start_stops.append(stop)
                            prevlength = prev_next_stop - prev_stop
                            first = False
                        break
                    break
        counter +=1
    if options.partial_storf:  # downstream partial StORF_Reporter - Last Stop to end of sequence
        try:
            if (len(sequence) - stop) > options.min_orf:
                seq = sequence[stop:len(sequence)]  # Start of seq to first stop identified
                ps_seq = cut_seq(seq, '+')
                storfs.update({",".join([str(stop), str(len(sequence))]): [ps_seq, str(frame), strand, stop,'Partial-StORF',StORF_idx]})
                StORF_idx +=1
        except UnboundLocalError:
            pass
    return storfs,con_StORFs,frames_covered,counter,lengths,StORF_idx,Con_StORF_idx

def STORF_Finder(sequence, sequence_id, options): #Main Function
    stops = []
    frames_covered = OrderedDict()
    for x in range (1,7):
        frames_covered.update({x: 0})
    for stop_codon in options.stop_codons.split(','):#Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence)]
    stops.sort()
    storfs = OrderedDict()
    con_StORFs = OrderedDict()
    counter = 0
    lengths = []
    StORF_idx = 0
    Con_StORF_idx = 0
    storfs,con_StORFs,frames_covered,counter,lengths,StORF_idx,Con_StORF_idx = find_storfs("positive",stops,sequence,storfs,con_StORFs,frames_covered,counter,lengths,'+',StORF_idx,Con_StORF_idx,options)
    ###### Reversed
    sequence_rev = revCompIterative(sequence)
    stops = []
    for stop_codon in options.stop_codons.split(','):#Find all Stops in seq
        stops += [match.start() for match in re.finditer(re.escape(stop_codon), sequence_rev)]
    stops.sort()
    counter = 0
    storfs,con_StORFs,frames_covered,counter,lengths,StORF_idx,Con_StORF_idx = find_storfs("negative",stops,sequence_rev,storfs,con_StORFs,frames_covered,counter,lengths,'-',StORF_idx,Con_StORF_idx,options)
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

    #Check if there are StORFs to report
    if options.con_storfs == False and options.con_only == False:
        if bool(storfs):
            if options.filtering == 'hard':
                all_StORFs = tile_filtering(storfs,options)
                all_StORFs = OrderedDict(sorted(all_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
            else:
                all_StORFs = OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
            # Reorder by start position

            if options.reporter == True:
                return all_StORFs
            ###Data Prepare
            gff_entries, fasta_entries = prepare_out(options, all_StORFs, sequence_id)
            write_fasta(fasta_entries, fasta_out, aa_fasta_out)
            if not options.aa_only:
                write_gff(gff_entries, gff_out)
    elif options.con_only == False:
        all_StORFs = {**storfs, **con_StORFs}
        if bool(storfs):
            if options.filtering == 'hard':
                all_StORFs = tile_filtering(storfs, options)
                all_StORFs = OrderedDict(sorted(all_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
            else:
                all_StORFs = OrderedDict(sorted(storfs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
            #### reporter  -
            if options.reporter == True:
                return all_StORFs
            ###Data Prepare
            gff_entries, fasta_entries = prepare_out(options, all_StORFs, sequence_id)
            write_fasta(fasta_entries, fasta_out, aa_fasta_out)
            if not options.aa_only:
                write_gff(gff_entries, gff_out)
    elif options.con_only == True and bool(con_StORFs):
        if options.filtering == 'hard':
            con_StORFs = tile_filtering(con_StORFs, options)
            con_StORFs = OrderedDict(sorted(con_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
        else:
            con_StORFs = OrderedDict(sorted(con_StORFs.items(), key=lambda e: tuple(map(int, e[0].split(",")))))
        #### reporter  -
        if options.reporter == True:
            return con_StORFs
        ###Data Prepare
        gff_entries, fasta_entries = prepare_out(options, con_StORFs, sequence_id)
        write_fasta(fasta_entries, fasta_out, aa_fasta_out)
        if not options.aa_only:
            write_gff(gff_entries, gff_out)
    elif options.verbose == True:
        print("No StOFS Found")

def fasta_load(fasta_in):
    first = True
    for line in fasta_in:
        line = line.strip()
        if line.startswith(';'):
            continue
        elif line.startswith('>') and not first:
            sequences.update({sequence_name: seq})
            seq = ''
            sequence_name = line
        elif line.startswith('>'):
            seq = ''
            sequence_name = line
        else:
            seq += str(line)
            first = False
    sequences.update({sequence_name: seq})


def StORF_Reported(Contigs,options):
    Reporter_StORFs = collections.OrderedDict()
    for Contig_ID, Contig_URs in Contigs.items():
        Reporter_StORFs.update({Contig_ID:[]})
        URs = Contig_URs[3]
        try:
            for UR in URs:
                if len(URs[UR][1]) >= options.min_orf: # Here
                    StORFs = STORF_Finder(URs[UR][1], Contig_ID, options)

                    if StORFs: #  Left out for now to allow for tracking of non-StORF URs
                        for StORF in StORFs.values():
                            StORF.append(URs[UR][0]) # True UR
                            StORF.append(UR) # Extended UR
                        Reporter_StORFs[Contig_ID].append(StORFs)
        except TypeError:
            if options.verbose == True:
                print("No URs in seq")
    return Reporter_StORFs




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF-Reporter v0.5.0: StORF_Finder Run Parameters.')
    parser.add_argument('-reporter', action="store", dest='reporter', default=False, required=False,
                        help=argparse.SUPPRESS)
    parser.add_argument('-f', action="store", dest='fasta', required=True,
                        help='Input FASTA File')
    parser.add_argument('-ua', dest='unannotated', action='store', default=True, type=eval, choices=[True, False],
                        help='Default - Treat input as Unannotated: Use "-ua False" for standard fasta')
    parser.add_argument('-wc', action="store", dest='whole_contig', default=False, type=eval, choices=[True, False],
                        help='Default - False: StORFs reported across entire sequence')
    parser.add_argument('-ps', action="store", dest='partial_storf', default=False, type=eval, choices=[True, False],
                        help='Default - False: Partial StORFs reported')
    parser.add_argument('-filt', action='store', dest='filtering', default='hard', const='hard', nargs='?',
                        choices=['none', 'soft', 'hard'],
                        help='Default - "hard": Filtering level "none" is not recommended, "soft" for single strand filtering '
                             'and hard for both-strand longest-first tiling')
    parser.add_argument('-aa', action="store", dest='translate', default=False, type=eval, choices=[True, False],
                        help='Default - False: Report StORFs as amino acid sequences')
    parser.add_argument('-con_storfs', action="store", dest='con_storfs', default=False, type=eval, choices=[True, False],
                        help='Default - False: Output Consecutive StORFs')
    parser.add_argument('-aa_only', action="store", dest='aa_only', default=False, type=eval, choices=[True, False],
                        help='Default - False: Only output Amino Acid Fasta')
    parser.add_argument('-con_only', action="store", dest='con_only', default=False, type=eval, choices=[True, False],
                        help='Default - False: Only output Consecutive StORFs')
    parser.add_argument('-stop_ident', action="store", dest='stop_ident', default=True, type=eval, choices=[True, False],
                        help='Default - True: Identify Stop Codon positions with \'*\'')
    parser.add_argument('-type', action='store', dest='feature_type', default='StORF', const='StORF', nargs='?',
                        choices=['StORF', 'CDS', 'ORF'],
                        help='Default - "StORF": Which GFF feature type for StORFs to be reported as in GFF')
    parser.add_argument('-minorf', action="store", dest='min_orf', default=100, type=int,
                        help='Default - 100: Minimum StORF size in nt')
    parser.add_argument('-maxorf', action="store", dest='max_orf', default=50000, type=int,
                        help='Default - 50kb: Maximum StORF size in nt')
    parser.add_argument('-codons', action="store", dest='stop_codons', default="TAG,TGA,TAA",
                        help='Default - (\'TAG,TGA,TAA\'): List Stop Codons to use')
    parser.add_argument('-olap', action="store", dest='overlap_nt', default=50, type=int,
                        help='Default - 50: Maximum number of nt of a StORF which can overlap another StORF.')
    parser.add_argument('-gff', action='store', dest='gff', default=True, type=eval, choices=[True, False],
                        help='Default - True: StORF Output a GFF file')
    parser.add_argument('-s', action="store", dest='suffix', required=False,
                        help='Default - Do not append suffix to genome ID')
    parser.add_argument('-so', action="store", dest='storf_order', default='start_pos', nargs='?', choices=['start_pos','strand'],
                        required=False,
                        help='Default - Start Position: How should StORFs be ordered when >1 reported in a single UR.')
    parser.add_argument('-o', action="store", dest='out_file', required=False,
                        help='Default - False: Without filetype - default appends \'_StORF-R\' to end of input gff filename (replaces \'.gff\')')
    parser.add_argument('-af', action="store", dest='affix', required=False,
                        help='Default - None: \'-af Con-StORFs\' can be used to append an identifier to output filename '
                             'to distinguish Con-StORF from StORF runs)')
    parser.add_argument('-lw', action="store", dest='line_wrap', default=True, type=eval, choices=[True, False],
                        help='Default - True: Line wrap FASTA sequence output at 60 chars')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    parser.add_argument('-v', action='store', dest='verbose', default='False', type=eval, choices=[True, False],
                        help='Default - False: Print out runtime status')
    parser.add_argument('-nout', action='store', dest='nout', default='False', type=eval, choices=[True, False],
                        help=argparse.SUPPRESS)
    options = parser.parse_args()
    ns_nt = defaultdict  # Used to Record non-standard nucleotides

    if options.out_file == None: # Clunky
        tmp_outfile = options.fasta.split('.')[-1]
        options.out_file = options.fasta.replace(tmp_outfile,'')
        options.out_file = options.out_file[:-1]

    if options.affix == None: # Also Clunky
        affix = '_StORF-R'
    else:
        affix = '_StORF-R_' + options.affix

    if not options.gz: # Clear fasta and gff files if not empty - Needs an elegant solution
        if not options.aa_only:
            gff_out = open(options.out_file + affix + '.gff', 'w', newline='\n', encoding='utf-8')
            gff_out.write("##gff-version\t3\n#\tStORF Stop - Stop ORF Predictions\n#\tRun Date:" + str(date.today()) + '\n')
            gff_out.write("##Original File: " + options.fasta + '\n')
            fasta_out = open(options.out_file + affix + '.fasta', 'w', newline='\n', encoding='utf-8')
            if options.translate:
                aa_fasta_out = open(options.out_file + affix + '_aa.fasta', 'w', newline='\n', encoding='utf-8')
        elif options.aa_only:
            aa_fasta_out = open(options.out_file + affix + '_aa.fasta', 'w', newline='\n', encoding='utf-8')
    elif options.gz:
        if not options.aa_only:
            gff_out = gzip.open(options.out_file + affix + '.gff.gz', 'wt', newline='\n', encoding='utf-8')
            gff_out.write("##gff-version\t3\n#\tStORF Stop - Stop ORF Predictions\n#\tRun Date:" + str(date.today()) + '\n')
            gff_out.write("##Original File: " + options.fasta + '\n')
            fasta_out = gzip.open(options.out_file + affix + '.fasta.gz', 'wt', newline='\n', encoding='utf-8')
            if options.translate:
                aa_fasta_out = gzip.open(options.out_file + affix + '_aa.fasta.gz', 'wt', newline='\n', encoding='utf-8')
        elif options.aa_only:
            aa_fasta_out = gzip.open(options.out_file + affix + '_aa.fasta.gz', 'wt', newline='\n', encoding='utf-8')

    sequences = OrderedDict()
    if options.reporter == False:
        try: # Detect whether fasta files are .gz or text and read accordingly
            fasta_in = gzip.open(options.fasta,'rt')
            fasta_load(fasta_in)
        except:
            fasta_in = open(options.fasta,'r')
            fasta_load(fasta_in)
        if options.verbose == True:
            print(fasta_in.name)

    for sequence_id, sequence in sequences.items():
        if len(sequence) >= options.min_orf:
            STORF_Finder(sequence, sequence_id, options)

