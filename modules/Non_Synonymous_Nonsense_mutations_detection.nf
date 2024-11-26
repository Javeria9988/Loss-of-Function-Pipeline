process Non_Synonymous_Nonsense_mutations_detection {

    input:
    path reference_file  
    path snps_consensus_fasta

    output:
    path("${snps_consensus_fasta.baseName}_Non_Synonymous_Nonsense_mutations.csv"), emit: stop_codons

    script:
    """
    python3 << 'EOF'
from Bio import SeqIO

# Function to detect premature stop codons by comparing with the reference
def detect_premature_stop_codons(reference_fasta, consensus_fasta, output_file, stop_codons=["TAA", "TAG", "TGA"]):
    premature_stop_codons = []

    # Parse the FASTA files
    reference_genome = SeqIO.to_dict(SeqIO.parse(reference_fasta, "fasta"))
    consensus_genome = SeqIO.to_dict(SeqIO.parse(consensus_fasta, "fasta"))
    
    # Compare each gene in the reference and consensus
    for gene_id, reference_record in reference_genome.items():
        if gene_id in consensus_genome:
            consensus_record = consensus_genome[gene_id]
            
            # Get sequences
            reference_sequence = str(reference_record.seq)
            consensus_sequence = str(consensus_record.seq)
            
            # Scan in codons
            for i in range(0, len(reference_sequence) - 2, 3):
                ref_codon = reference_sequence[i:i+3]
                cons_codon = consensus_sequence[i:i+3]
                
                # Check if the reference codon is not a stop codon but the consensus is
                if cons_codon in stop_codons and ref_codon not in stop_codons:
                    premature_stop_codons.append((gene_id, i + 1, ref_codon, cons_codon))
    
    # Write results to the output CSV
    with open(output_file, 'w') as out_file:
        out_file.write("Gene_ID,Stop_Codon_Position,Reference_Codon,Consensus_Codon\\n")
        for gene, pos, ref_codon, cons_codon in premature_stop_codons:
            out_file.write(f"{gene},{pos},{ref_codon},{cons_codon}\\n")
    
    if premature_stop_codons:
        print(f"Detected premature stop codons in {len(premature_stop_codons)} positions. Details written to {output_file}.")
    else:
        print("No premature stop codons detected.")

# Input and output files
consensus_fasta = "${snps_consensus_fasta}"
reference_fasta = "${reference_file}"
output_file = "${snps_consensus_fasta.baseName}_Non_Synonymous_Nonsense_mutations.csv"

# Run the detection function
detect_premature_stop_codons(reference_fasta, consensus_fasta, output_file)
EOF
    """
}

