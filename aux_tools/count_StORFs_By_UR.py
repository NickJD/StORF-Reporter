import gzip
import glob
import argparse
import collections

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Overlapping CDS Gene Lengths.')
    parser.add_argument('-ur', action="store", dest='ur_in', required=True,
                        help='Input Directory containing URs')
    parser.add_argument('-storf', action="store", dest='storf_in', required=True,
                        help='Input Directory containing StORFs')

    options = parser.parse_args()

    UR_list = list(glob.glob(options.ur_in + '/*.fa_UR.gff.gz'))
    StORF_list = list(glob.glob(options.storf_in + '/*fa_UR_StORFs.gff.gz'))

    current_URs = collections.defaultdict(list)#lambda: collections.defaultdict(list))

    for ur_file in UR_list:

        current_identifier = ur_file.split('/')[-1].split('.dna.toplevel.fa')[0]
        ur_in = gzip.open(ur_file, 'rt')
        for line in ur_in:
            if not line.startswith('#') and not line.startswith('\n'):
                UR = line.split('\t')[8].split('ID=')[1].split(';')[0]
                current_URs[current_identifier+'|'+UR] = 0
        storf_in = gzip.open(ur_file.replace('.fa_UR.gff.gz','.fa_UR_StORFs.gff.gz').replace(options.ur_in,options.storf_in),'rt')
        for line in storf_in:
            if not line.startswith('#') and not line.startswith('\n'):
                UR = line.split('\t')[8].split('ID=')[1].split(';')[0].split('_StORF')[0]
                current_URs[current_identifier + '|' + UR] +=1


    res = collections.defaultdict(int)
    for key, val in current_URs.items():
        res[val] += 1

    # printing result
    print("The frequency dictionary : " + str(dict(res)))






