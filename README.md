# Loss_of_Function-Pipeline

## Overview
This pipeline performs a comprehensive pangenome analysis, starting from genome annotation, through pangenome identification, fastq reads generation, SNP calling, and finally, detecting synonymous and non-synonymous mutations. The pipeline leverages various bioinformatics tools, integrated within a Nextflow workflow to streamline the analysis of multiple sequences.

## Pipeline Structure
The pipeline consists of the following steps:

1. **Annotate:** Annotates the input sequences using BAKTA.
2. **Panaroo:** Identifies core and accessory genes across the annotated genomes and generates relevant outputs.
3. **DWGSIM:** Generates error-free fastq reads from the input sequences.
4. **Snippy:** Aligns sequences and calls SNPs based on the simulated reads.
5. **Non_Synonymous_Nonsense_mutations_detection:** Identifies non-sense mutations i.e. premature stop codons from the SNP data using python code.
6. **Non_Synonymous_Missense_mutations_detection:** Identifies missense mutations from the SNP data using python code.
7. **Synonymous_mutations_detection:** Identifies synonymous mutations from the SNP data using python code.

## Workflow
![Pipeline workflow](https://github.com/user-attachments/assets/4092bbdf-4c1d-4fd8-a8e0-6ec7b7430fe7)

## Installation

### Prerequisites
- [Nextflow](https://www.nextflow.io/)
- [Docker](https://www.docker.com/) for containerized execution of processes.

### Setup
Clone the repository and navigate to the directory:
```bash
git clone https://github.com/Javeria9988/Loss-of-Function-Pipeline.git
cd Loss-of-Function-Pipeline
```
### Usage
### Running the Pipeline
To run the pipeline, use the following command:
nextflow run main.nf -c nextflow.config -with-report workflow_report.html -with-trace trace.txt -with-timeline timeline.html

### Parameters
params.sequences: Path to the input sequences in .fa format. The default is set to 'sequences/*.fa'.

### Input Data
Sequences: The pipeline expects input sequences in .fa format, located in the directory specified by the params.sequences parameter.

### Output Data
Annotated Genomes: GFF files generated by the annotation step.
Panaroo Outputs: Core and accessory gene information in CSV and FASTA formats.
fastq Reads generation: Error-free reads generated by DWGSIM.
SNPs: SNP calling results from Snippy.
Non_Synonymous_Nonsense_mutations_detection_results: A csv file detailing any detected premature stop codons in each sample.
Non_Synonymous_Missense_mutations_detection_results: A csv file detailing any detected missense mutations in each sample.
Synonymous_mutations_detection_results: A csv file detailing any detected synonymous mutations in each sample.


## Docker Integration
This pipeline utilizes Docker containers to ensure reproducibility and consistency across different environments. Each step of the pipeline is associated with a specific Docker container, as outlined below:

### Containers Used:

1. **Annotate:**
   - **Container:** `staphb/bakta:1.9.4-5.1-light`
   - **Description:** This container includes BAKTA, a tool used for rapid annotation of prokaryotic genomes.

2. **Panaroo:**
   - **Container:** `staphb/panaroo:latest`
   - **Description:** This container includes Panaroo, a tool for pangenome analysis, particularly focused on bacterial genomes.

3. **DWGSIM:**
   - **Container:** `dwgsim-with-ps:latest`
   - **Description:** This customized container includes DWGSIM, a tool for getting error-free reads from genomic FASTA sequences and ps, which is required for genration of nextflow report.

4. **Snippy:**
   - **Container:** `staphb/snippy:4.6.0`
   - **Description:** This container includes Snippy, a tool for rapid bacterial SNP calling and variant detection.

5. **Non_Synonymous_Nonsense_mutations_detection:**
   - **Container:** `javeriam21/python-biopython-pandas-ps:latest`
   - **Description:** This customized container includes Python along with the ps, Biopython and Pandas libraries, used for detecting premature stop codons and processing genomic data.

6. **Non_Synonymous_Missense_mutations_detection:**
   - **Container:** `javeriam21/python-biopython-pandas-ps:latest`
   - **Description:** This customized container includes Python along with the ps, Biopython and Pandas libraries, used for detecting missense mutations and processing genomic data.
   
7. **Synonymous_mutations_detection:**
   - **Container:** `javeriam21/python-biopython-pandas-ps:latest`
   - **Description:** This customized container includes Python along with the ps, Biopython and Pandas libraries, used for detecting synonymous mutations and processing genomic data.


### Usage
To ensure that Docker is enabled, the pipeline is configured with Docker integration:

```groovy
docker {
    enabled = true
}
