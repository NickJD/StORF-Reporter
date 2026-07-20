nextflow.enable.dsl=2

params.samples = 'samples.tsv'
params.outdir = 'results'
params.primary_tier = 'exploratory'
params.gene_ident = 'CDS'
params.gff_coordinate_mode = 'stop-to-stop'
params.overlap_coordinate_basis = 'coding'
params.include_consecutive = false
params.include_first_stop = false
params.remove_inframe_stops = false
params.protein_coordinate_mode = 'candidate'

process STORF_AUGMENT {
    tag "$sample"
    publishDir "${params.outdir}/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(gff), path(fasta), val(genetic_code)

    output:
    path "out/${sample}/**"

    script:
    def con_option = params.include_consecutive ? '--include-consecutive' : ''
    def first_stop_option = params.include_first_stop ? '--include-first-stop' : ''
    def remove_stops_option = params.remove_inframe_stops ? '--remove-inframe-stops' : ''
    """
    storf augment --input ${gff} --fasta ${fasta} --outdir out \
      --genetic-code ${genetic_code} --gene-ident '${params.gene_ident}' \
      --primary-tier ${params.primary_tier} \
      --gff-coordinate-mode ${params.gff_coordinate_mode} \
      --protein-coordinate-mode ${params.protein_coordinate_mode} \
      --overlap-coordinate-basis ${params.overlap_coordinate_basis} ${con_option} \
      ${first_stop_option} ${remove_stops_option} --overwrite
    """
}

workflow {
    Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.sample, file(row.gff), file(row.fasta), row.genetic_code ?: 11) }
        .set { samples }

    STORF_AUGMENT(samples)
}
