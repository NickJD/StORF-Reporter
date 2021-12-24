import re
import argparse
import numpy as np
import collections

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqs', default='', help='Which seqs to analyse?')
args = parser.parse_args()

##########################
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
######################


def seq_Lengths(seqs):
    count = 0
    seq = ""
    first = True
    stop_codons = collections.OrderedDict()
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")

            if line.startswith('>') and first == False:
                for idx, stop_codon in enumerate("TGA,TAG,TAA".split(',')):  # Find all Stops in seq
                    stop_codons[genome][idx] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq)])
                seq_rev = revCompIterative(seq)
                for idx, stop_codon in enumerate("TGA,TAG,TAA".split(',')):  # Find all Stops in seq
                    stop_codons[genome][idx] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq_rev)])
                seq = ""
            if not line.startswith('>'):
                seq += str(line)
                first = False
            if line.startswith('>'):
                genome = line.split('|')[0]
                if '>Escherichia_coli_bw25113.asm75055v1' in genome:
                    print("")
                if genome not in stop_codons.keys():
                    stop_codons[genome] = [0,0,0]
                    count +=1
                    print(count)
                    # if count == 2:
                    #     break
    for idx, stop_codon in enumerate("TGA,TAG,TAA".split(',')):  # Find all Stops in seq
        stop_codons[genome][idx] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq)])
    seq_rev = revCompIterative(seq)
    for idx, stop_codon in enumerate("TGA,TAG,TAA".split(',')):  # Find all Stops in seq
        stop_codons[genome][idx] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq_rev)])

    TGAs = []
    TAGs = []
    TAAs = []
    for key, values in stop_codons.items():
        TGAs.append(values[0])
        TAGs.append(values[1])
        TAAs.append(values[2])

    print(TGAs)
    print(TAGs)
    print(TAAs)

    print("Plotting")

    TGA_P = []
    TAG_P = []
    TAA_P = []

    for idx, TGA in enumerate(TGAs):
        all = TGA+TAGs[idx]+TAAs[idx]
        TGA_P.append((int(TGA/all)*100))
        TAG_P.append((int(TAGs[idx] / all) * 100))
        TAA_P.append((int(TAAs[idx] / all) * 100))

    print(TGA_P)
    print(TAG_P)
    print(TAA_P)


    #print("Number of codons: " + str(collections.Counter(stop_codons('TGA'))))
    # print("Median of Seqs: " + str(np.median(lengths)))
    # print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    # print("Std: " + format(np.std(lengths),'.2f'))

if __name__ == "__main__":
    seq_Lengths(**vars(args))










