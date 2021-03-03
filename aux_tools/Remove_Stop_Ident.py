import gzip
import argparse




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF Run Parameters.')
    parser.add_argument('-seq', action="store", dest='fasta', required=True,
                        help='Input Sequence File')
    parser.add_argument('-o', action="store", dest='out_prefix', required=True,
                        help='Output file prefix - Without filetype')
    parser.add_argument('-gz', action='store', dest='gz', default='False', type=eval, choices=[True, False],
                        help='Default - False: Output as .gz')
    options = parser.parse_args()