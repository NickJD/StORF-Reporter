import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-ur', '--URs', required= True, help='Which URs to analyse?')
args = parser.parse_args()

def seq_Lengths(seqs):
    seq = ""
    lengths = []
    prev_ID = ""
    with open(seqs, 'r') as seqs:
        for line in seqs:
            line = line.replace("\n","")
            if line.startswith('>'):
                # if len(seq) >= 100000:
                #     print(prev_ID)
                if len(seq) != 0:
                    lengths.append(len(seq))
                    if len(seq) <=50:
                        print("")
                seq = ""
                prev_ID = line
            if not line.startswith('>'):
                seq += str(line)



    print("Number of Seqs: " + str(len(lengths)))
    print("Longest Seq: " + str(max(lengths)))
    print("Shortest Seq: " + str(min(lengths)))
    print("Median of Seqs: " + str(np.median(lengths)))
    print("Mean of Seqs: " + format(np.mean(lengths), '.2f'))
    print("Std: " + format(np.std(lengths),'.2f'))
    print("25th Percentile " + format(np.percentile(lengths, 25), '.2f'))
    print("75th Percentile " + format(np.percentile(lengths,75),'.2f'))


if __name__ == "__main__":
    seq_Lengths(**vars(args))
