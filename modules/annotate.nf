process annotate {
    input:
    path contigfile

    output:
    path "${contigfile.simpleName}.gff3", emit: gff_files

    script:
    """
    # Create output directory if it doesn't exist
    mkdir -p ${contigfile.simpleName}

    # Run Bakta annotation with the specific contig file and database path
    bakta --db ${params.db_path} --output ${contigfile.simpleName} --prefix ${contigfile.simpleName} ${contigfile} --force

    # Copy the generated GFF3 file to the desired location
    cp ${contigfile.simpleName}/${contigfile.simpleName}.gff3 ${contigfile.simpleName}.gff3
    """
}

