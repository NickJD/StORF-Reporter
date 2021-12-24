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
    seq = ""
    first = True
    stop_codons = collections.OrderedDict({"TGA":0,"TAG":0,"TAA":0})
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")
            if line.startswith('>') and first == False:
                for stop_codon in "TGA,TAG,TAA".split(','):  # Find all Stops in seq
                    stop_codons[stop_codon] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq)])
                seq_rev = revCompIterative(seq)
                for stop_codon in "TGA,TAG,TAA".split(','):  # Find all Stops in seq
                    stop_codons[stop_codon] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq_rev)])

                seq = ""
            if not line.startswith('>'):
                seq += str(line)
                first = False
    for stop_codon in "TGA,TAG,TAA".split(','):  # Find all Stops in seq
        stop_codons[stop_codon] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq)])
    seq_rev = revCompIterative(seq)
    for stop_codon in "TGA,TAG,TAA".split(','):  # Find all Stops in seq
        stop_codons[stop_codon] += len([match.start() for match in re.finditer(re.escape(stop_codon), seq_rev)])

    print(stop_codons)
    TGA = stop_codons["TGA"]
    TAG = stop_codons["TAG"]
    TAA = stop_codons["TAA"]
    All = TGA+TAG+TAA
    print((TGA / All) * 100)
    print((TAG/ All) * 100)
    print((TAA / All) * 100)
    #print("Number of codons: " + str(collections.Counter(stop_codons('TGA'))))
    # print("Median of Seqs: " + str(np.median(lengths)))
    # print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    # print("Std: " + format(np.std(lengths),'.2f'))

if __name__ == "__main__":
    seq_Lengths(**vars(args))










