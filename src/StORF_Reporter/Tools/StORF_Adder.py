import argparse


try:
    from ORForise.GFF_Adder import gff_adder  # Calling from ORForise via pip
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    from ORForise.src.ORForise.GFF_Adder import gff_adder  # Calling from ORForise locally (StORF_Reporter and ORForise in same dir)





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='StORF-Reporter v0.5.4: StORF-Adder Run Parameters.')
    parser.add_argument('-dna', '--genome_DNA', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    parser.add_argument('-rt', '--reference_tool', required=False,
                        help='Which tool format to use as reference? - If not provided, will default to '
                             'standard Ensembl GFF format, can be Prodigal or any of the other tools available')
    parser.add_argument('-ref', '--reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')
    parser.add_argument('-gi', '--gene_ident',  default='CDS', required=False,
                        help='Identifier used for extraction of "genic" regions from reference annotation '
                             '"CDS,rRNA,tRNA": Default for is "CDS"')
    parser.add_argument('-at', '--additional_tool', default='StORF_Reporter', required=False, # needs cleaner fix
                        help='Which format to use for additional annotation?')
    parser.add_argument('-add', '--additional_annotation', required=True,
                        help='Which annotation file to add to reference annotation?')
    parser.add_argument('-olap', '--overlap', default=50, type=int, required=False,
                        help='Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Output filename')

    options = parser.parse_args()
    gff_adder(options)
