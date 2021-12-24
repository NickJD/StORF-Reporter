import argparse
import numpy as np
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqs', default='', help='Which seqs to analyse?')
args = parser.parse_args()

def add_stops(stop_codons,genome,stop):
    if "TGA" in stop:
        stop_codons[genome][0] +=1
    elif "TAG" in stop:
        stop_codons[genome][1] += 1
    elif "TAA" in stop:
        stop_codons[genome][2] += 1
    return stop_codons

def seq_Lengths(seqs):
    seq = ""
    first = True
    genome = ""
    stop_codons = collections.defaultdict()
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")

            if line.startswith('>') and first == False:
                start = seq[:3]
                end = seq[-3:]
                stop_codons = add_stops(stop_codons,genome,end)
                stop_codons = add_stops(stop_codons, genome, start)

                seq = ""
            if not line.startswith('>'):
                seq += str(line)
                first = False

            if line.startswith('>'):
                genome = line.split('|')[0]
                if genome not in stop_codons.keys():
                    stop_codons[genome] = [0,0,0]

    start = seq[:3]
    end = seq[-3:]
    #stop_codons[start] += 1
    stop_codons = add_stops(stop_codons,genome,end)
    stop_codons = add_stops(stop_codons, genome, start)
    print(stop_codons)

    TGAs = []
    TAGs = []
    TAAs = []
    for key, values in stop_codons.items():
        TGAs.append(values[0])
        TAGs.append(values[1])
        TAAs.append(values[2])
    #print("Number of codons: " + str(collections.Counter(stop_codons('TGA'))))
    # print("Median of Seqs: " + str(np.median(lengths)))
    # print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    # print("Std: " + format(np.std(lengths),'.2f'))

    print(TGAs)
    print(TAGs)
    print(TAAs)

if __name__ == "__main__":
    seq_Lengths(**vars(args))










