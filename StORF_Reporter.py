import argparse
from argparse import Namespace
import pathlib


try:
    from UR_Extractor import extractor
    from StORF_Finder import StORF_Reported
except ImportError:
    from .UR_Extractor import extractor
    from .StORF_Finder import StORF_Reported


def run_UR_Extractor(options):
    fasta = list(pathlib.Path(options.prokka_dir).glob('*.fna'))
    fasta = str(fasta[0])
    gff = list(pathlib.Path(options.prokka_dir).glob('*.gff'))
    gff = str(gff[0])

    if options.gene_ident == "PROKKA":
        gene_ident = "misc_RNA,gene,mRNA,CDS,tRNA,tmRNA,CRISPR"

        options = Namespace(fasta=fasta, gff=gff, gene_ident=gene_ident,
                            minlen=30, maxlen=100000, exlen=50 ,verbose=False, nout=True)

    URs = extractor(options)
    return URs,gff


def find_prev_StORFs(all_StORFs,track_current_start,track_prev_stop):
    StORFs_to_del = []
    StORFs = []
    is_break = False
    for UR, UR_StORFs in all_StORFs.items():
        if is_break == True:
            break
        if UR == '11736_12213':
            print("Double")
        if UR_StORFs != None:
            for UR_StORF, UR_StORF_data in UR_StORFs.items():
                StORF_Pos_In_UR = UR_StORF
                StORF_Start_In_UR = int(StORF_Pos_In_UR.split(',')[0])
                StORF_Stop_In_UR = int(StORF_Pos_In_UR.split(',')[1])
                UR_StORF_list = UR_StORF_data
                strand = UR_StORF_list[2]
                key_start = int(UR.split('_')[0])
                StORF_start = key_start + StORF_Start_In_UR
                if 7937 == StORF_start:
                    print("D")
                key_stop = int(UR.split('_')[1])
                StORF_stop = key_start + StORF_Stop_In_UR
                allow_start = StORF_start + allowed_overlap
                allow_stop = StORF_stop - allowed_overlap
                if allow_stop <= track_current_start:# and allow_start >= track_prev_stop:
                    StORFs.append([StORF_start,StORF_stop,StORF_Start_In_UR,StORF_Stop_In_UR,strand,UR_StORF_list[0][0],UR_StORF_list[0][3],UR_StORF_list[0][0]])
                    if UR not in StORFs_to_del:
                        StORFs_to_del.append(UR)
                elif allow_start > track_current_start: # is this fucking with URs with more than one StORF?
                    is_break = True
                    break

    for StORF in StORFs_to_del:
        del all_StORFs[StORF]
### Might need to force remove more ur storfs regions
    return StORFs,all_StORFs

def find_after_StORFs(all_StORFs,track_current_start,track_current_stop):
    StORFs = []
    for key, value in all_StORFs.items():
        if value != None:
            StORF_Pos_In_UR = list(value.keys())
            StORF_Start_In_UR = int(StORF_Pos_In_UR[0].split(',')[0])
            StORF_Stop_In_UR = int(StORF_Pos_In_UR[0].split(',')[1])
            value_list = value.values()
            value_list = list(value_list)
            strand = value_list[0][2]
            key_start = int(key.split('_')[0])
            StORF_start = key_start + StORF_Start_In_UR
            key_stop = int(key.split('_')[1])
            StORF_stop = key_stop + StORF_Stop_In_UR
            allow_start = StORF_start + allowed_overlap
            allow_stop = StORF_stop - allowed_overlap


            if allow_start >= track_current_stop:
                StORFs.append([StORF_start,StORF_stop,StORF_Start_In_UR,StORF_Stop_In_UR,strand,value_list[0][0],value_list[0][3],value_list[0][0]])

    return StORFs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-pd', '--PROKKA_dir', action='store', dest='prokka_dir', default='', required=False,
                        help='PROKKA output directory to be used')
    parser.add_argument('-gene_ident', action='store', dest='gene_ident', default='PROKKA', const='PROKKA', required=False,
                        choices=['PROKKA', 'Ensem', 'CDS'], nargs='?',
                        help='Identifier used for extraction of "unannotated" regions ("CDS,rRNA,tRNA"):'
                             ' Default for PROKKA = "ID=gene"')

    options = parser.parse_args()
    if options.prokka_dir != "":
        URs,gff = run_UR_Extractor(options)

        options = Namespace(reporter=True, stop_codons="TGA,TAA,TAG", partial_storf=False, whole_contig=False,
                            con_storfs=False, con_only=False, max_orf=50000, filtering='hard',
                            overlap_nt=50,
                            minlen=30, maxlen=100000, min_orf=100 ,verbose=False, nout=True)

        sequence_id, all_StORFs = StORF_Reported(URs, options)

        for k, v in list(all_StORFs.items()):
            if v is None:
                del all_StORFs[k]


        gff_in = open(gff,'r')

        outfile = open(gff+'_StORF.gff','w')

        track_prev_start, track_prev_stop = 0, 0
        track_current_loci, track_current_URs = 0, 0
        StORF_Num = 0
        allowed_overlap = 50
        end = False

        for line in gff_in:
            if not line.startswith('#') and end == False:
                data = line.split('\t')
                track_current_start = int(data[3])
                track_current_stop = int(data[4])
                if track_current_start == track_prev_start and track_current_stop == track_prev_stop: # `duplicate' entry in GFF
                    print("skip")
                else:
                    StORFs,all_StORFs = find_prev_StORFs(all_StORFs,track_current_start,track_prev_stop)
                track_prev_start = track_current_start
                track_prev_stop = track_current_stop
                if StORFs:
                    for StORF in StORFs: # ([StORF_start,StORF_stop,StORF_Start_In_UR,StORF_Stop_In_UR,strand,value_list[0][0],value_list[0][3],value_list[0][0]])
                        outfile.write(sequence_id + '\tStORF-Reporter\tCDS\t' + str(StORF[0]) + '\t' + str(StORF[1]) + '\t.\t' + StORF[4] + '\t.\tID=Additional_Annotation_StORF-Reporter;StORF_Num_' + str(StORF_Num) + '\n')
                        StORF_Num +=1
                outfile.write(line)
                StORFs = None

            elif line.startswith('##FASTA'):
                end = True
                StORFs = find_after_StORFs(all_StORFs, track_current_start, track_current_stop)
                if StORFs:
                    for StORF in StORFs: # ([StORF_start,StORF_stop,StORF_Start_In_UR,StORF_Stop_In_UR,strand,value_list[0][0],value_list[0][3],value_list[0][0]])
                        outfile.write(sequence_id + '\tStORF-Reporter\tCDS\t' + str(StORF[0]) + '\t' + str(StORF[1]) + '\t.\t' + StORF[4] + '\t.\tID=Additional_Annotation_StORF-Reporter;StORF_Num_' + str(StORF_Num) + '\n')
                        StORF_Num += 1

            else:
                outfile.write(line)







