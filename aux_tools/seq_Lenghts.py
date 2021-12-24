import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--seqs', default='', help='Which seqs to analyse?')
args = parser.parse_args()

def seq_Lengths(seqs):
    seq = ""
    lengths = []
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")
            if line.startswith('>'):
                if len(seq) != 0:
                    lengths.append(len(seq))
                seq = ""
            if not line.startswith('>'):
                seq += str(line)

    print("Number of Seqs: " + str(len(lengths)))
    print("Median of Seqs: " + str(np.median(lengths)))
    print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    print("Std: " + format(np.std(lengths),'.2f'))
    print("Min: "+ str(np.min(lengths)))
    print("Max: "+ str(np.max(lengths)))

if __name__ == "__main__":
    seq_Lengths(**vars(args))
