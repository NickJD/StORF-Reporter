import argparse

import collections


def load_stops(FASTA):
    results_in = open(FASTA,'r')
    Con_StORFs = collections.defaultdict(list)
    Con_StORFs_Counting = collections.defaultdict(int)
    codons = ["TGA","TAG","TAA"]
    multi = 0
    for line in results_in:
        if line.startswith('>'):
            id = line.replace('>','')
            id = id.replace('\n','')
            stops = line.split('UR_Stop_Locations:')[1].split('|')[0].split('-')
            first_stop = stops[0]
            stops = stops[1:]
            del stops[-1]
        else:
            count = 0
            for stop in stops:
                stop = int(stop) - int(first_stop) # Need to account for the UR postional reporting of the stops
                tmp = line[stop-3:stop]
                if len(stops) >1:
                    multi +=1
                if any(x == line[stop-3:stop] for x in codons):  # fixes the bug in con-storf and just be after on newer output
                    Con_StORFs[id].append([stop,line[stop-3:stop]])
                    Con_StORFs_Counting[line[stop-3:stop]] += 1
                else:
                    tmp = line[stop:stop+3]
                    Con_StORFs[id].append([stop,line[stop:stop+3]])
                    Con_StORFs_Counting[line[stop:stop+3]] +=1
                count +=1
    return  Con_StORFs,Con_StORFs_Counting


def get_matching_codon(codon_stop_info,stop,first_stop):
    for li in codon_stop_info:
        li[0] = li[0] + int(first_stop)
        if stop == li[0]:
            matching_codon = li[1]
            break
    return matching_codon


def load_results(Con_StORF_Stops,results_input):
    results_in = open(results_input,'r')
    Con_StORFs = collections.OrderedDict()
    for line in results_in:
        line = line.split('\t')
        try: # This is where I can connect the blast output and the DNA Con-StORF_Reporter internal codons data...
            codon_stop_info = Con_StORF_Stops[line[0]]
            stops = line[0].split('UR_Stop_Locations:')[1].split('|')[0].split('-')
            stops = [int(x) for x in stops] # Yes very stupid
            storf_start = line[6]
            storf_end = line[7]
            left = int(storf_start)*3+stops[0]
            right = int(storf_end)*3+stops[0]
            Con_StORF_Alignments = []
            for stop in stops:
                if left < stop and right > stop:
                    matching_codon = get_matching_codon(codon_stop_info,stop,stops[0])
                    Con_StORF_Alignments.append(matching_codon)
            Con_StORFs.update({line[0]:Con_StORF_Alignments})
        except KeyError:
            continue

    return  Con_StORFs

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', action='store', dest='BLAST_Output', required=True,
                        help='Diamond/BLAST output to check for TRUE Con-StORFs')
    parser.add_argument('-f', action='store', dest='FASTA', required=True,
                        help='Con-StORF FASTA')
    options = parser.parse_args()

    Con_StORF_Stops, Con_StORF_Stops_Counting = load_stops(options.FASTA)


    Con_StORFs = load_results(Con_StORF_Stops,options.BLAST_Output)

    value_list = []
    multi_counting = 0

    for key, values in Con_StORFs.items():
        print(key+" : "+str(values))
        if len(values) > 1:
            multi_counting +=1
        for value in values:
            value_list.append(value)
    print(len(Con_StORFs))
    validated = list(filter(None,Con_StORFs.values()))
    validated = [item for sublist in validated for item in sublist]
    print("Number of Validated Con-StORF: " + str(len(validated)))
    print("TGA: " + str(value_list.count("TGA")) + "\nTAG: " +
          str(value_list.count("TAG")) + "\nTAA: " + str(value_list.count("TAA")))
    print("MULTI: " + str(multi_counting))