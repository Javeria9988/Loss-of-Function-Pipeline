#!/usr/bin/env nextflow

include { annotate } from './modules/annotate.nf'
include { panaroo } from './modules/panaroo.nf'
include { dwgsim } from './modules/dwgsim.nf'
include { snippy } from './modules/snippy.nf'
include { detectStopCodons } from './modules/detectStopCodons.nf'
include { genedescription } from './modules/genedescription.nf'
include { aggregateStopCodonData } from './modules/aggregateStopCodonData.nf'

workflow {
    sequence_ch = Channel.fromPath("${params.sequences}")
    
    // Annotate sequences
    annotate_ch = annotate(sequence_ch)

    // Extract the GFF files and TXT files
    gff_files_ch = annotate_ch.gff_files.collect()

    // Run Panaroo and get the reference FASTA output
    panaroo_outputs = panaroo(gff_files_ch)  

    // Generate error-free reads using dwgsim
    dwgsim_outputs = dwgsim(sequence_ch)
    
    // Run snippy on dwgsim outputs
    snippy_outputs = snippy(panaroo.out.pan_genome_reference, dwgsim_outputs)

    // Premature stop codons detection
    detectStopCodons_outputs = detectStopCodons(snippy.out.consensus_fasta)
    
    //Addition of gene description
    genedescription_output = genedescription(panaroo.out.gene_presence_absence, detectStopCodons_outputs)

    // Aggregate stop codon data across samples into a matrix 
    matrix_output = aggregateStopCodonData(genedescription.out.desc_files.collect())
    
}

