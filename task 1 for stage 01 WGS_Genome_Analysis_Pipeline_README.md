
# **WGS: Whole Genome Sequencing (WGS) Analysis Pipeline for Patient X (human)**

#======================================================================================

### #Module 00: Setting

Create a folder to do the task:
```bash
mkdir human
cd human
```

Create a new conda environment with the necessary tools:
```bash
conda create -n genome_analysis -c bioconda -c conda-forge gatk bwa samtools fastqc -y
conda activate genome_analysis
```

Check tools versions:
```bash
echo "Checking GATK version..."
gatk --version

echo "Checking BWA version..."
bwa --version

echo "Checking Samtools version..."
samtools --version

echo "Checking FastQC version..."
fastqc --version
```

==========================================================================================================================

### #Model 1: Download Raw Data, Check Quality Control and Preprocessing

```bash
mkdir raw_reads && cd raw_reads
nano rawreads.sh
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR315/030/SRR31548730/SRR31548730_1.fastq.gz
bash rawreads.sh
```

The first step in processing the raw sequencing data is quality control (QC). We will use **FastQC** to check the quality of the raw FASTQ file and ensure there are no obvious issues with the sequencing data.

Create a directory to store QC reports:
```bash
mkdir qc
```

Run FastQC to assess the quality of the raw fastq data:
```bash
fastqc /data/human_stage_1/SRR31548730_1.fastq.gz -O qc
```

Explanation:
- **fastqc**: This tool assesses the quality of the sequencing data (e.g., base quality, GC content).
- **-O qc**: This flag specifies the output directory for the results. The output will be a report in HTML format that you can use to evaluate the raw data quality.

======================================================================================

### #Model 2: Trimming Off Low Quality Reads

```bash
nano trim.sh
mkdir trimm
input_file="SRR31548730_1.fastq.gz"
output_dir="trim"
fastp -i SRR31548730_1.fastq.gz -o output_clean.fastq.gz   -h fastp_report.html   -j fastp_report.json   --cut_front   --cut_tail   --cut_mean_quality 20   --length_required 50
```

Explanation of flags:
- **-i input.fastq.gz**: The input raw FASTQ file.
- **-o output_clean.fastq.gz**: The output cleaned FASTQ file after quality control.
- **-h fastp_report.html**: Generates an HTML report of the quality control process.
- **-j fastp_report.json**: Generates a JSON report.
- **--cut_front**: Trims bases from the 5' (front) end of the read until it meets the quality threshold.
- **--cut_tail**: Trims bases from the 3' (tail) end of the read until it meets the quality threshold.
- **--cut_mean_quality 20**: Specifies the mean quality score threshold (in this case, 20). Bases below this score will be trimmed.
- **-length_required 50**: Discards reads that are shorter than 50 bases after trimming.

```bash
cp trim.sh raw_reads && cd raw_reads
bash trim.sh
cp -R raw_reads/trim
```

======================================================================================

### #Model 3: Genome Mapping

To make the assembly, we need the reference sequence to align the genome. We can download it from NCBI.

We will use **BWA** to align sequences and create the genome mapping.

1. **Create Reference Genome Index with BWA**

```bash
bwa index /data/ref/Homo_sapiens_assembly38.fasta
```

This will generate the following index files:
- **Homo_sapiens_assembly38.fasta.bwt**
- **Homo_sapiens_assembly38.fasta.pac**
- **Homo_sapiens_assembly38.fasta.sa**
- **Homo_sapiens_assembly38.fasta.amb**

2. **Genome Mapping (BWA)**

Create a directory to store the mapped reads (BAM files):

```bash
mkdir ~/habiba/humantask/mapped
```

Run **BWA** to align the reads to the reference genome:

```bash
bwa mem -R '@RG	ID:patientX	SM:patientX	PL:ILLUMINA' ~/habiba/humantask/reference/Homo_sapiens_assembly38.fasta ~/habiba/humantask/raw_reads/output_clean.fastq.gz | samtools view -b -o ~/habiba/humantask/mapped/patientX_aligned.bam
bash align.sh
```

Explanation:
- **-R '@RG	ID:patientX	SM:patientX	PL:ILLUMINA'**: Read group information.
- **~/habiba/humantask/reference/Homo_sapiens_assembly38.fasta**: The reference genome path.
- **patientX_aligned.bam**: The output BAM file containing aligned reads.

3. **Sorting the BAM File**

```bash
mkdir -p ~/habiba/humantask/sorted
gatk SortSam -I ~/habiba/humantask/mapped/patientX_aligned.bam -O ~/habiba/humantask/sorted/patientX_sorted.bam -SORT_ORDER coordinate
```

4. **Marking Duplicates**

```bash
mkdir -p ~/habiba/humantask/marked
gatk MarkDuplicates -I ~/habiba/humantask/sorted/patientX_sorted.bam -O ~/habiba/humantask/marked/patientX_marked.bam -M ~/habiba/humantask/marked/patientX_metrics.txt
```

5. **Building BAM Index**

```bash
gatk BuildBamIndex -I ~/habiba/humantask/marked/patientX_marked.bam
```

======================================================================================

### #Model 4: Base Quality Score Recalibration (BQSR)

```bash
mkdir -p ~/habiba/humantask/BQSR
gatk BaseRecalibrator -I ~/habiba/humantask/marked/patientX_marked.bam -R /data/ref/Homo_sapiens_assembly38.fasta --known-sites ~/habiba/humantask/reference/Homo_sapiens_assembly38.known_indels.vcf.gz -O ~/habiba/humantask/BQSR/patientX_recal.table
gatk ApplyBQSR -I ~/habiba/humantask/marked/patientX_marked.bam -R /data/ref/Homo_sapiens_assembly38.fasta --bqsr-recal-file ~/habiba/humantask/BQSR/patientX_recal.table -O ~/habiba/humantask/BQSR/patientX_recal.bam
```

======================================================================================

### #Model 5: Variant Calling Using HaplotypeCaller

```bash
mkdir -p ~/habiba/humantask/VCF
gatk HaplotypeCaller -I ~/habiba/humantask/BQSR/patientX_recal.bam   -R /home/sararomi/reference/Homo_sapiens_assembly38.fasta   -O ~/habiba/humantask/VCF/patientX.g.vcf.gz   -ERC GVCF
```

======================================================================================

### #Model 6: Annotating Variants with VEP

```bash
vep -i ~/habiba/humantask/VCF/patientX.g.vcf.gz --cache --dir_cache /path/to/vep/cache --output_file ~/habiba/humantask/VCF/patientX_annotated.vcf
```

Or, upload your final VCF to Variant Effect Predictor: https://www.ensembl.org/Multi/Tools/VEP?db=core

#Explanation:
#vep: Variant Effect Predictor annotates the VCF file by providing details about the functional consequences of each variant.
#Output: Annotated VCF file (patientX_annotated.vcf).
# or Your final VCF can  be uploaded to Variant Effect Predictor: https://www.ensembl.org/Multi/Tools/VEP?db=core
# searched for different genes 
#1. rs34126315
#Gene: HBB (Hemoglobin subunit beta)
#Chromosome: 11
#Location: chr11:5,225,469 (hg38)
#Mutation: A → T (in the 6th codon)
#Amino Acid Change: Glu6Val (Glutamic acid to Valine)
#Clinical Significance: This is the canonical sickle cell mutation (also known as HbS). This mutation causes the red blood cells to become sickle-shaped, leading to sickle cell disease. It is a pathogenic variant.
#2. rs34135787
#Gene: HBB (Hemoglobin subunit beta)
#Chromosome: 11
#Location: chr11:5,225,531 (hg38)
#Mutation: C → T (in the 6th codon)
#Amino Acid Change: Glu6Lys (Glutamic acid to Lysine)
#Clinical Significance: This mutation is another form of the sickle cell mutation, where the mutation causes HbS formation, leading to sickle cell disease. This is also a pathogenic variant.
#rs33941849 - Key Information
#Gene: HBB (Hemoglobin subunit beta)
#Chromosome: 11
#Location: chr11:5,225,518 (hg38)
#Mutation: A → T
#Amino Acid Change: Glu6Val (Glutamic acid to Valine)
#Clinical Significance: This variant is also known as HbS and is associated with sickle cell disease (SCD). This mutation causes a glutamic acid to valine substitution at position 6 of the hemoglobin beta chain, leading to sickle-shaped red blood cells and reduced oxygen transport.
#=====================================================================

## **Links**

- **GitHub**: [GitHub Link](https://github.com/habibafaried/hackbio-bash-script-habiba/blob/main/human%20task%20hackbio%20stage%201)
- **LinkedIn**: [LinkedIn Profile](https://www.linkedin.com/posts/habiba-faried-a62bb11b0_hackbio-bash-linux-activity-7366582020027703296-s0Vk?utm_source=share&utm_medium=member_android&rcm=ACoAADFMi6oBwjx2cbYE0OvQ1vQcA7FnBUzmvb)
