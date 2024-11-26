#!/usr/bin/env nextflow

include { annotate } from './modules/annotate.nf'
include { panaroo } from './modules/panaroo.nf'
include { dwgsim } from './modules/dwgsim.nf'
include { snippy } from './modules/snippy.nf'
include { Non_Synonymous_Nonsense_mutations_detection } from './modules/Non_Synonymous_Nonsense_mutations_detection.nf'
include { Non_Synonymous_Missense_mutations_detection } from './modules/Non_Synonymous_Missense_mutations_detection.nf'
include { Synonymous_mutations_detection } from './modules/Synonymous_mutations_detection.nf'


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

    // Premature stop codons i.e. non-sense mutations detection
    Non_Synonymous_Nonsense_mutations_detection_outputs = Non_Synonymous_Nonsense_mutations_detection(panaroo.out.pan_genome_reference, snippy.out.consensus_fasta)
    
    // Missense mutations detection
    Non_Synonymous_Missense_mutations_detection_outputs = Non_Synonymous_Missense_mutations_detection(panaroo.out.pan_genome_reference, snippy.out.consensus_fasta)
    
    //Synonymous mutations detection
    Synonymous_mutations_detection_outputs = Synonymous_mutations_detection(panaroo.out.pan_genome_reference, snippy.out.consensus_fasta)
    
}

