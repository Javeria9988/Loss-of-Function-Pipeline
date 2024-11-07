process panaroo {

    input:
    path gff_files

    output:
    path "panaroo_output/pan_genome_reference.fa", emit: pan_genome_reference
    path "panaroo_output/gene_presence_absence.csv", emit: gene_presence_absence
    path "panaroo_output/summary_statistics.txt", emit: summary_statistics

    script:
    """
    panaroo -i ${gff_files.join(' ')} \
            -o panaroo_output \
            --clean-mode ${params.panaroo.clean_mode} \
            -a ${params.panaroo.analysis_mode} \
            --aligner ${params.panaroo.aligner} \
            -t ${params.panaroo.threads}
    """
}

