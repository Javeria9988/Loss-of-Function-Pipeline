process genedescription {
    input:
    path gene_presence_absence_csv
    path detect_stop_codon_csv

    output:
    path "${detect_stop_codon_csv.baseName}_gene_description.csv", emit: desc_files

    script:
    """
    python3 << EOF
import pandas as pd

# Load gene annotations from the gene_presence_absence.csv
gene_annotations = pd.read_csv("${gene_presence_absence_csv}")

# Replace empty annotations with 'No annotation'
gene_annotations['Annotation'] = gene_annotations['Annotation'].fillna('No annotation')

# Load stop codon file (ensure it's a CSV with the appropriate separator)
stop_codons = pd.read_csv("${detect_stop_codon_csv}", header=None, names=['Gene_Name', 'Stop_Codon_Position', 'Original_Sequence'])

# Merge the two dataframes on the gene name
merged_data = pd.merge(stop_codons, gene_annotations[['Gene', 'Annotation']], left_on='Gene_Name', right_on='Gene', how='left')

# Replace NaN values in annotations with 'No annotation'
merged_data['Annotation'] = merged_data['Annotation'].fillna('No annotation')

# Create a new DataFrame for the output in CSV format with the desired column order
output_data = merged_data[['Gene_Name', 'Stop_Codon_Position', 'Annotation', 'Original_Sequence']].copy()

# Replace 'No annotation' with 'Product' in the first row of the annotation column
if not output_data.empty:
    output_data.iloc[0, 2] = 'Product'  # Assuming Annotation is now the third column (index 2)

# Write the result to a new CSV file without headers
output_data.to_csv("${detect_stop_codon_csv.baseName}_gene_description.csv", index=False, header=False)

EOF
    """
}

