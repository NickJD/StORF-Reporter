import argparse
import numpy as np
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqs', default='', help='Which seqs to analyse?')
args = parser.parse_args()

def seq_Lengths(seqs):
    seq = ""
    first = True
    non_canonical = []
    stop_codons = collections.OrderedDict({"TGA":0,"TAG":0,"TAA":0})
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")
            if line.startswith('>') and first == False:
                start = seq[:3]
                end = seq[-3:]
                #stop_codons[start] +=1
                try:
                    stop_codons[end] +=1
                except KeyError:
                    non_canonical.append(end)
                seq = ""
            if not line.startswith('>'):
                seq += str(line)
                first = False
    start = seq[:3]
    end = seq[-3:]
    #stop_codons[start] += 1
    try:
        stop_codons[end] += 1
    except KeyError:
        non_canonical.append(end)
    print(stop_codons)
    TGA = stop_codons["TGA"]
    TAG = stop_codons["TAG"]
    TAA = stop_codons["TAA"]
    All = TGA+TAG+TAA
    print((TGA / All) * 100)
    print((TAG/ All) * 100)
    print((TAA / All) * 100)
    print("Number of non-canonial: " +str(len(non_canonical)))
    #print("Number of codons: " + str(collections.Counter(stop_codons('TGA'))))
    # print("Median of Seqs: " + str(np.median(lengths)))
    # print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    # print("Std: " + format(np.std(lengths),'.2f'))

if __name__ == "__main__":
    seq_Lengths(**vars(args))










