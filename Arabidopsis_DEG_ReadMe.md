
# Project: Differential Expression Analysis of Arabidopsis thaliana Leaf Vasculature in Response to UV Stress

This project involves analyzing the differential gene expression in the vasculature of Arabidopsis thaliana leaf tissues treated with UV-C stress compared to a control (Water-treated) group. The goal is to identify genes that respond to UV stress and the enriched pathways in the UV-C treated group.

## Modal 00

---

## Modal 01: Data Preprocessing
We will begin by downloading the FASTQ files, performing quality control (QC) on the raw data, and trimming the reads for adapter removal and low-quality bases.

### Step 1: Download FASTQ Files: The following FASTQ files correspond to the Arabidopsis thaliana leaf vasculature:
- Control: SRR12808527, SRR12808528, and SRR12808529
- UV-C Treated: SRR12808497, SRR12808498, and SRR12808499

```bash
mkdir clean_reads
#!/bin/bash
# Download the FASTQ files using wget
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/027/SRR12808527/SRR12808527.fastq.gz
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/029/SRR12808529/SRR12808529.fastq.gz
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/097/SRR12808497/SRR12808497.fastq.gz
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/098/SRR12808498/SRR12808498.fastq.gz
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/099/SRR12808499/SRR12808499.fastq.gz
wget -P clean_data ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/028/SRR12808528/SRR12808528.fastq.gz
```

---

## Modal 02: Quality Control and Multiple Quality Control

We will use FastQC to check the quality of the raw reads.

### Run FastQC on all downloaded FASTQ files:

```bash
mkdir -p qc/raw
fastqc -o qc/raw clean_reads/*.fastq.gz
```

We will run MultiQC to aggregate the individual FastQC reports and generate a summary.

```bash
multiqc -o qc/multiqc qc/raw qc/trimmed
```

---

## Modal 03: Trimming the Reads (Fastp)

After checking the quality, we will trim the reads using fastp to remove adapters and low-quality bases.

```bash
#!/bin/bash
# Create a directory for trimmed files if it doesn't exist
mkdir -p trim

# Loop through all .fastq.gz files in clean_data/
for filename in clean_data/*.fastq.gz; do
    # Extract the base filename without extension
    base=$(basename "$filename" .fastq.gz)

    # Print the fastp command for debugging (optional)
    echo fastp -i "$filename" -o "trim/${base}_trim.fastq.gz" -h "trim/${base}_report.html"

    # Run fastp to trim the file and output to the 'trim' directory
    fastp -i "$filename" -o "trim/${base}_trim.fastq.gz" -h "trim/${base}_report.html"
done
```

---

## Modal 04: Quality Control of Trimmed Reads

Run FastQC again on the trimmed data.

```bash
mkdir -p qc/trimmed
fastqc -o qc/trimmed trimmed/*.fastq.gz
```

---

## Modal 05: Download Reference Genome

### RNA-Seq Alignment to the Reference Genome

For Arabidopsis thaliana, we will download the reference genome and GTF file, then align the trimmed reads to the genome using STAR.

```bash
mkdir genome
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa a_thaliana.fa
```

### Build STAR Genome Index

```bash
mkdir -p genomeIndex
STAR --runMode genomeGenerate      --genomeDir genomeIndex      --genomeFastaFiles genome/a_thaliana.fa
```

### Download Annotated Genome

```bash
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz
gunzip Arabidopsis_thaliana.TAIR10.62.gff3.gz
mv Arabidopsis_thaliana.TAIR10.62.gff3 a_thaliana.gff3
```

---

## Modal 06: Align the Trimmed Reads with STAR

Now align the trimmed reads to the reference genome using STAR.

```bash
#!/bin/bash
# Create directory for mapped files
mkdir -p mapped
# Align each trimmed FASTQ file using STAR
for infile in trimmed/*.fastq.gz; do
    outfile=$(basename "$infile" .fastq.gz)
    STAR --genomeDir references/STARindex          --readFilesIn "$infile"          --outFileNamePrefix mapped/$outfile          --outSAMtype BAM SortedByCoordinate          --outSAMattributes All
done
```

---

## Modal 07: Counts

```bash
mkdir -p counts
# Run featureCounts
featureCounts -O -t gene -g ID -a a_thaliana.gff3 -o counts/counts.txt mapped/*_Aligned.sortedByCoord.out.bam
# Check the summary
cat counts/counts.txt.summary
```

---

## Modal 08: Differential Expression with DESeq2

### Set Working Directory:

```r
setwd("/path/to/your/directory")  # Adjust to your working directory
```

### Libraries

```r
library(DESeq2)        # For differential expression analysis
library(pheatmap)      # For heatmap visualizations
library(dplyr)         # For data manipulation
library(ggplot2)       # For plotting
library(clusterProfiler)  # For enrichment analysis
library(org.At.tair.db)   # For Arabidopsis ID mapping
```

---

## Results and Discussions

### 1. Volcano Plot
The volcano plot illustrates the global transcriptional changes in response to UV-C treatment. A large number of genes were significantly differentially expressed, with a clear asymmetry between upregulated and downregulated genes.

### 2. Heatmap
This heatmap confirms that UV-C treatment has a strong global impact on gene expression. UV-C treated samples cluster tightly, showing reproducibility of the transcriptional response.

### 3. Interpretation of the Top DEGs and GO Enrichment
Gene Ontology (GO) enrichment analysis revealed terms related to response to oxidative stress, oxidative stress, and reactive oxygen species generation.

### 4. KEGG Pathway Enrichment
The KEGG analysis confirmed that UV-C treatment induces strong stress-responsive pathways in Arabidopsis thaliana.

---

## Links

1. [LinkedIn](https://www.linkedin.com/posts/habiba-faried-a62bb11b0_hackbio-bash-linux-activity-7366582020027703296-s0Vk?utm_source=share&utm_medium=member_android&rcm=ACoAADFMi6oBwjx2cbYE0OvQ1vQcA7FnBUzmvb)
2. [GitHub](https://github.com/habibafaried/hackbio-bash-script-habiba/blob/main/plant%20rnaseq)
