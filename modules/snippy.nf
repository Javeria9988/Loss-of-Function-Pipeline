process snippy {

    input:
    path reference_file
    tuple val(sample_name), path(fq1), path(fq2)

    output:
    path "snippy_output/${sample_name}/${sample_name}_snps.consensus.fa", emit: consensus_fasta
    path "snippy_output/${sample_name}/${sample_name}.vcf", emit: variant_calling_files
    path "snippy_output/${sample_name}/${sample_name}_snps.txt", emit: snps_file

    script:
    """
    snippy --outdir snippy_output/${sample_name} \
           --ref ${reference_file} \
           --R1 ${fq1} \
           --R2 ${fq2} \
           --mincov ${params.snippy.mincov} \
           --minfrac ${params.snippy.minfrac} \
           --maxsoft ${params.snippy.maxsoft} \
           --minqual ${params.snippy.minqual} \
           --cpus ${params.snippy.cpus}

    # Rename the snps.consensus.fa file to include the sample name
    mv snippy_output/${sample_name}/snps.consensus.fa snippy_output/${sample_name}/${sample_name}_snps.consensus.fa
    
    # Rename the VCF file to include the sample name
    mv snippy_output/${sample_name}/snps.vcf snippy_output/${sample_name}/${sample_name}.vcf
    
    # Rename the snps.txt file to include the sample name
    mv snippy_output/${sample_name}/snps.txt snippy_output/${sample_name}/${sample_name}_snps.txt
    """
}

