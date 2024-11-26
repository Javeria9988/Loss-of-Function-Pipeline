process dwgsim {

    input:
    path sequence_file

    output:
    tuple val(sequence_file.simpleName), path("dwgsim_output/${sequence_file.simpleName}_1.fq"), path("dwgsim_output/${sequence_file.simpleName}_2.fq")

    script:
    """
    mkdir -p dwgsim_output
    dwgsim -e ${params.dwgsim.error_rate} -E ${params.dwgsim.indel_fraction} -N ${params.dwgsim.num_reads} -1 ${params.dwgsim.read_len_1} -2 ${params.dwgsim.read_len_2} ${sequence_file} dwgsim_output/${sequence_file.simpleName}

    if [ -f dwgsim_output/${sequence_file.simpleName}.bwa.read1.fastq ]; then
        mv dwgsim_output/${sequence_file.simpleName}.bwa.read1.fastq dwgsim_output/${sequence_file.simpleName}_1.fq
    fi

    if [ -f dwgsim_output/${sequence_file.simpleName}.bwa.read2.fastq ]; then
        mv dwgsim_output/${sequence_file.simpleName}.bwa.read2.fastq dwgsim_output/${sequence_file.simpleName}_2.fq
    fi
    """
}

