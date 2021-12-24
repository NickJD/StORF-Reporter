import argparse

import collections




def load_results(results_input):
    results_in = open(results_input,'r')
    Con_StORFs = collections.defaultdict(list)
    Con_StORFs_Counting = collections.defaultdict(int)
    codons = ["TGA","TAG","TAA"]
    multi = 0
    for line in results_in:
        if line.startswith('>'):
            id = line
            stops = line.split('UR_Stop_Locations:')[1].split('|')[0].split('-')
            first_stop = stops[0]
            stops = stops[1:]
            del stops[-1]
        else:
            count = 0
            if len(stops) > 1:
                multi += 1
            for stop in stops:
                stop = int(stop) - int(first_stop) # Need to account for the UR postional reporting of the stops
                tmp = line[stop-3:stop]
                if any(x == line[stop-3:stop] for x in codons):  # fixes the bug in con-storf and just be after on newer output
                    Con_StORFs[id].append(line[stop-3:stop])
                    Con_StORFs_Counting[line[stop-3:stop]] += 1
                else:
                    tmp = line[stop:stop+3]
                    Con_StORFs[id].append(line[stop:stop+3])
                    Con_StORFs_Counting[line[stop:stop+3]] +=1
                count +=1
            #print(id + '\t' + str(Con_StORFs[id]))



    print("Multi: " + str(multi))
    return  Con_StORFs,Con_StORFs_Counting

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input', required=True,
                        help='Con-StORFs')
    options = parser.parse_args()

    Con_StORFs,Con_StORFs_Counting = load_results(options.input)

    print(collections.Counter(Con_StORFs_Counting))

