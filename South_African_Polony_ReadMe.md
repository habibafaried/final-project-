
# **South African Polony — Listeria monocytogenes Outbreak Genomics**

## **Executive Summary (TL;DR)**

I downloaded 40 paired-end Illumina datasets linked to the 2017 South African *Listeria monocytogenes* polony outbreak, performed quality control (QC) and trimming, assembled the data with SPAdes, compared the assemblies using QUAST, inspected assembly graphs with Bandage, and screened the assemblies for antimicrobial resistance (AMR) and virulence genes using ABRicate.

The organism was confirmed as *Listeria monocytogenes* through BLAST of circular components or representative contigs. AMR screening identified **fosX** (fosfomycin resistance) and **lmo0919_fam** (ABC-F ribosomal protection; lincosamides). Virulence profiling revealed universal **LIPI-1** (prfA, hly, plcA, plcB, actA) and widespread **LIPI-3** (llsA–llsP), indicating hypervirulence.

**Therapeutic Implication**: Treat invasive disease with **ampicillin** (or **penicillin G**) ± **gentamicin**, and avoid **cephalosporins** (intrinsic resistance), **lincomycin/clindamycin** (ABC-F), and **fosfomycin** (fosX). **TMP-SMX** is a common alternative for β-lactam allergies.

---

## **Methods**

### **1. Environment Setup, Data Acquisition, and QC**

1. **Create and activate a conda environment:**
    ```bash
    conda create -n finaltask -c bioconda -c conda-forge spades quast fastqc abricate -y
    conda activate finaltask
    ```

2. **Download paired-end dataset:**
    ```bash
    mkdir raw_reads && cd raw_reads
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR270/016/SRR27013316/SRR27013316_1.fastq.gz -o SRR27013316_1.fastq.gz
    curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR270/016/SRR27013316/SRR27013316_2.fastq.gz -o SRR27013316_2.fastq.gz
    # Repeat for all other files...
    ```

3. **Run FastQC for quality check:**
    ```bash
    mkdir fastqc_results
    fastqc raw_reads/*.fastq.gz -o fastqc_results/
    ```

### **2. Trimming Adapters and Low-Quality Bases**

1. **Create a trimming script (`trim.sh`):**
    ```bash
    mkdir -p trim
    SAMPLES=("SRR27013316" "SRR27013315" "SRR27013314" # List all samples)

    for SAMPLE in "${SAMPLES[@]}"; do
      fastp -i "${SAMPLE}_1.fastq.gz" -I "${SAMPLE}_2.fastq.gz" -o "trim/${SAMPLE}_1.trim.fastq.gz" -O "trim/${SAMPLE}_2.trim.fastq.gz" --html "trim/${SAMPLE}_fastp.html" --json "trim/${SAMPLE}_fastp.json"
    done
    ```

### **3. Assembly and Assessment**

1. **SPAdes Assembly:**
    ```bash
    set -u
    TRIM_DIR="raw_reads/trim"
    OUT_BASE="assembly"
    mkdir -p "$OUT_BASE"

    for R1 in "$TRIM_DIR"/*_1.trim.fastq.gz; do
      [[ -e "$R1" ]] || { echo "No *_1.trim.fastq.gz files found"; break; }
      STEM="$(basename "$R1" "_1.trim.fastq.gz")"
      R2="$TRIM_DIR/${STEM}_2.trim.fastq.gz"
      spades.py -1 "$R1" -2 "$R2" -k 41,85 -o "$OUTDIR" > "$OUTDIR/${STEM}.spades.log" 2>&1
    done
    ```

2. **QUAST Assembly Evaluation:**
    ```bash
    FILES=(assembly/*/contigs.fasta)
    LABELS=$(for f in "${FILES[@]}"; do basename "$(dirname "$f")"; done | paste -sd, -)
    quast.py "${FILES[@]}" -o quast_all -t 8 -l "$LABELS"
    ```

### **4. Organism Identification (BLAST)**

1. **Graph Inspection with Bandage:**
    - Visualize the SPAdes assembly graph using **Bandage**.
    - Export circular components and BLAST them against the NCBI database to confirm **Listeria monocytogenes**.

    ```bash
    Bandage image assembly_group/assembly_graph.fastg graph_overview.png --drawdepth --labels
    ```

### **5. AMR Screening**

1. **Run ABRicate:**
    ```bash
    abricate --db ncbifiles/ "$contigs" > AMR/"$sample"_result.tab
    ```

2. **Summarize Results:**
    ```bash
    abricate --summary AMR/*.tab > AMR/abricate_summary.tsv
    ```

### **6. Virulence Profiling**

1. **Screen against VFDB:**
    ```bash
    for d in assembly/*/; do
      contigs="$d/contigs.fasta"
      abricate --db vfdb "$contigs" > "virulence/${sample}_vfdb.tab"
    done
    ```

2. **Create a summary table for virulence genes:**
    ```bash
    awk 'BEGIN{FS=OFS="	"} !/^#/ {if($6 ~ /(^hly$|^prfA$|^plcA$|^plcB$|^actA$|^inlA$|^inlB$|^lls)/) print s,$6,$14,$12,$10 }' "$f" > virulence/summary_key_virulence.tsv
    ```

---

## **Results**

1. **Organism Identification:**
    - **BLAST hits** for the circular components confirmed **Listeria monocytogenes** with 100% identity and coverage (e.g., CP096157.1).

2. **AMR Screening:**
    - Detected **fosX** (fosfomycin resistance) and **lmo0919_fam** (ABC-F ribosomal protection; lincosamides).
    - **Prevalence**: **fosX** and **lmo0919_fam** were detected in most isolates.
    - **Clinical Implications**: Avoid **fosfomycin**, **lincosamides** (e.g., clindamycin), and **cephalosporins** (intrinsic resistance).

3. **Virulence Profiling:**
    - **LIPI-1** (prfA, hly, plcA, plcB, actA) and **LIPI-3** (llsA–llsP) are present, indicating hypervirulence.
    - **Internalins** inlA and inlB were largely intact, with minor variation in some isolates.

---

## **Directory Structure**

```
south-african-polony/
├─ README.md
├─ env.yml
├─ data/
│  └─ raw_reads/                # (empty; populated by rawreads.sh)
├─ qc/
│  ├─ fastqc_results/           # generated
│  └─ trim/                     # generated by trim.sh (fastp)
├─ assembly/
│  └─ {SAMPLE}/                 # SPAdes outputs
├─ quast_all/                   # QUAST comparison output
├─ AMR/
│  ├─ *.tab                     # per-sample abricate outputs
│  └─ abricate_summary.tsv
├─ virulence/
│  ├─ *_vfdb.tab
│  ├─ summary_key_virulence.tsv
│  └─ virulence_matrix.tsv
├─ scripts/
│  ├─ rawreads.sh
│  ├─ trim.sh
│  ├─ assemble.sh
│  ├─ quast.sh
│  ├─ amr.sh
│  └─ virulence.sh
└─ figures/
   └─ graph_overview.png
```

---

## **Interpretation / Takeaways**

1. **LIPI-1** and **LIPI-3** are universally present, supporting a virulent Listeria strain.
2. **InlA** and **InlB** internalins are mostly intact, ensuring strong invasiveness.

---

## **Therapeutic Context**

- **AMR findings**: Avoid **fosfomycin** and **lincosamides**. Use **ampicillin** (or **penicillin G**) ± **gentamicin** for invasive listeriosis.
- **β-lactam allergy**: **TMP-SMX** is a common alternative.

---

## **GitHub Repository**: [South African Polony Outbreak Project](https://github.com/your-github-link)
