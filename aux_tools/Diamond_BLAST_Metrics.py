import numpy as np
import argparse

def diamond_Parser(options):
    diamond_In = open(options.diamond,'r')
    Pident = []
    Length = []
    Bit = []

    for line in diamond_In:
        line = line.strip()
        line = line.split('\t')
        Pident.append(float(line[2]))
        Length.append(float(line[3]))
        Bit.append(float(line[11]))

    print(np.average(Pident))
    print(np.average(Length))
    print(np.average(Bit))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF Run Parameters.')
    parser.add_argument('-diamond', action="store", dest='diamond', required=True,
                        help='Input Sequence File')
    options = parser.parse_args()
    diamond_Parser(options)