process aggregateStopCodonData {
    input:
    path desc_files  // Input is a list of files passed from previous processes

    output:
    path "gene_stop_codon_matrix.csv"  // Output is a CSV file with the matrix

    script:
    """
    python3 << EOF
import os
import pandas as pd

# Use the provided desc_files directly as a space-separated string
csv_files = "${desc_files.join(' ')}".split()

# Initialize an empty dictionary to store stop codon counts for each gene per sample
stop_codon_counts = {}

# Process each CSV file to count premature stop codons per gene
for csv_file in csv_files:
    # Get the sample name from the filename (extracting the first part before '_snps')
    sample_name = os.path.basename(csv_file).split('_snps')[0]  # Extract the part before '_snps'
    
    # Load the gene description file without headers
    data = pd.read_csv(csv_file, header=None, names=['Gene_Name', 'Stop_Codon_Position', 'Product', 'Original_Sequence'])
    
    # Count premature stop codons for each gene in this sample
    gene_stop_counts = data['Gene_Name'].value_counts()
    
    # Store counts in dictionary
    for gene, count in gene_stop_counts.items():
        if gene not in stop_codon_counts:
            stop_codon_counts[gene] = {}
        stop_codon_counts[gene][sample_name] = count

# Create a DataFrame from the dictionary
stop_codon_matrix = pd.DataFrame(stop_codon_counts).T.fillna(0).astype(int)
stop_codon_matrix.index.name = 'Gene'
stop_codon_matrix.reset_index(inplace=True)

# Write the matrix to a CSV file
stop_codon_matrix.to_csv("gene_stop_codon_matrix.csv", index=False)

EOF
    """
}

