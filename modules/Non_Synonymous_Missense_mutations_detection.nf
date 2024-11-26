process Non_Synonymous_Missense_mutations_detection {

    input:
    path reference_file                
    path snps_consensus_fasta

    output:
    path "${snps_consensus_fasta.baseName}_Non_Synonymous_Missense_mutations.csv", emit: missense_mutations

    script:
    """
    python3 << 'EOF'
    from Bio import SeqIO
    import csv

    # Function to translate codons to amino acids
    def translate_codon(codon):
        codon_table = {
            "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M", "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
            "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K", "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
            "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L", "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
            "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q", "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
            "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
            "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E", "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
            "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S", "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
            "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*", "TGA": "*", "TGC": "C", "TGT": "C", "TGG": "W",
            "CGC": "R", "CGA": "R", "CGT": "R", "TGC": "C", "TGT": "C", "TGG": "W"
        }
        return codon_table.get(codon, "*")  # Return '*' for stop codons

    # Function to detect non-synonymous missense mutations
    def detect_missense_mutations(reference_fasta, consensus_fasta, output_file):
        reference_genome = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
        consensus_genome = SeqIO.to_dict(SeqIO.parse(consensus_fasta, "fasta"))
        
        missense_mutations = []

        # Iterate over the sequences in the reference genome
        for gene_id, reference_record in reference_genome.items():
            if gene_id in consensus_genome:
                consensus_record = consensus_genome[gene_id]
                
                # Get the coding sequence of the gene
                reference_sequence = str(reference_record.seq)
                consensus_sequence = str(consensus_record.seq)
                
                # Compare the reference and consensus sequences in triplets (codons)
                for i in range(0, len(reference_sequence) - 2, 3):
                    ref_codon = reference_sequence[i:i+3]
                    cons_codon = consensus_sequence[i:i+3]
                    
                    # If the codons differ, it might indicate a missense mutation
                    if ref_codon != cons_codon:
                        ref_amino_acid = translate_codon(ref_codon)
                        cons_amino_acid = translate_codon(cons_codon)
                        
                        if ref_amino_acid != cons_amino_acid:
                            missense_mutations.append((gene_id, i + 1, ref_codon, cons_codon, ref_amino_acid, cons_amino_acid))

        # Write the detected missense mutations to the output file
        with open(output_file, 'w') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(["Gene_ID", "Position", "Reference_Codon", "Consensus_Codon", "Reference_Amino_Acid", "Consensus_Amino_Acid"])
            for mutation in missense_mutations:
                csv_writer.writerow(mutation)

        if missense_mutations:
            print(f"Detected {len(missense_mutations)} missense mutations. Results written to {output_file}.")
        else:
            print("No missense mutations detected.")
    
    consensus_fasta = "${snps_consensus_fasta}"   
    reference = "${reference_file}"              
    output_file = "${snps_consensus_fasta.baseName}_Non_Synonymous_Missense_mutations.csv"  

    # Run the missense mutation detection
    detect_missense_mutations(reference, consensus_fasta, output_file)
    EOF
    """
}

