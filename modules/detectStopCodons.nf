process detectStopCodons {

    input:
    path snps_consensus_fasta

    output:
    path("${snps_consensus_fasta.baseName}_premature_stop_codons.csv"), emit: stop_codons

    script:
    """
    python3 << 'EOF'
from Bio import SeqIO

def detect_premature_stop_codons(consensus_fasta, output_file, stop_codons=["TAA", "TAG", "TGA"]):
    premature_stop_codons = []

    # Parse the FASTA file
    alignment = SeqIO.parse(consensus_fasta, "fasta")
    
    for record in alignment:
        sequence = str(record.seq)
        gene_name = record.id
        
        stop_positions = []
        
        # Scan the sequence in steps of 3 nucleotides (codons) from the first position
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            if codon in stop_codons:  # Look for the stop codons
                stop_positions.append(i + 1)  # Record the 1-based position of the stop codon
        
        if stop_positions:
            for pos in stop_positions:
                premature_stop_codons.append((gene_name, pos, sequence))
    
    # Write the detected premature stop codons to the output file in CSV format
    with open(output_file, 'w') as out_file:
        out_file.write("Gene_Name,Stop_Codon_Position,Original_Sequence\\n")  # Change to CSV header
        for gene, pos, seq in premature_stop_codons:
            out_file.write(f"{gene},{pos},{seq}\\n")  # Change to CSV format
    
    if premature_stop_codons:
        print(f"Premature stop codons detected in {len(premature_stop_codons)} genes. Details written to {output_file}.")
    else:
        print("No premature stop codons detected.")

consensus_fasta = "${snps_consensus_fasta}"
output_file = "${snps_consensus_fasta.baseName}_premature_stop_codons.csv"
detect_premature_stop_codons(consensus_fasta, output_file)
EOF
    """
}

