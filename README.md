
# RNA-seq Practice Workflow (Paired-End)

This repository provides a **hands-on practice pipeline** for students learning RNA-seq data analysis.  
It uses **publicly available paired-end reads** (uploaded here as *control* and *treated* samples), for potato and covers every major preprocessing step — from quality control to gene-level read counting.

## Repository Contents

```

sample/
├── 50k_subsamples/              # Paired-end FASTQ files (sample data)
│   ├── control_rep1_1.50k.fastq
│   ├── control_rep1_2.50k.fastq
│   ├── control_rep2_1.50k.fastq
│   ├── control_rep2_2.50k.fastq
│   ├── control_rep3_1.50k.fastq
│   ├── control_rep3_2.50k.fastq
│   ├── treated_rep1_1.50k.fastq
│   ├── treated_rep1_2.50k.fastq
│   ├── treated_rep2_1.50k.fastq
│   ├── treated_rep2_2.50k.fastq
│   ├── treated_rep3_1.50k.fastq
│   └── treated_rep3_2.50k.fastq
├── trim_and_qc_one_sample.sh    # Automated trimming + QC script (optional)
└── README.md                    # This guide

````

## Quick Start

### Install Required Tools

Make sure these tools are installed and available in your PATH:

| Tool | Purpose | Install command (conda) |
|------|----------|--------------------------|
| **FastQC** | Quality check of reads | `conda install -c bioconda fastqc` |
| **Trim Galore** | Adapter and quality trimming | `conda install -c bioconda trim-galore cutadapt` |
| **HISAT2** | RNA-seq read alignment | `conda install -c bioconda hisat2` |
| **Samtools** | BAM/SAM processing | `conda install -c bioconda samtools` |
| **HTSeq** | Gene-level read counting | `conda install -c bioconda htseq` |

Alternatively, if Conda is not available:

```bash
sudo apt install fastqc hisat2 samtools htseq
python3 -m pip install --user cutadapt
````

---

## Example Workflow

Below is a complete example using one control sample (**control_rep1**).
Repeat for each sample pair (both control and treated replicates).

---

### **Step 1 — Quality Check (Raw Reads)**

```bash
fastqc control_rep1_1.50k.fastq control_rep1_2.50k.fastq -o fastqc_raw
```

**Explanation**

* `control_rep1_1.50k.fastq` and `control_rep1_2.50k.fastq` are forward and reverse reads.
* `-o fastqc_raw` specifies an output directory for FastQC reports.
* Open `fastqc_raw/control_rep1_1.50k_fastqc.html` in a web browser to review read quality.

---

### **Step 2 — Trim Reads**

```bash
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
  control_rep1_1.50k.fastq control_rep1_2.50k.fastq \
  -o trimmed_reads
```

**Explanation**

* `--paired` → indicates paired-end sequencing data.
* `--clip_R1/2 10` → trims 10 bases from the 5′ end of each read.
* `--three_prime_clip_R1/2 5` → trims 5 bases from the 3′ end of each read.
* `-o trimmed_reads` → directory for trimmed reads and reports.

**Output Files**

```
trimmed_reads/
 ├── control_rep1_1.50k.fastq_trimming_report.txt
 ├── control_rep1_2.50k.fastq_trimming_report.txt
 ├── control_rep1_1.50k_val_1.fq
 └── control_rep1_2.50k_val_2.fq
```

---

### **Step 3 — Re-run QC on Trimmed Reads**

```bash
fastqc trimmed_reads/control_rep1_1.50k_val_1.fq trimmed_reads/control_rep1_2.50k_val_2.fq -o fastqc_trimmed
```

This step confirms that the **“Per base sequence content”** and overall quality improved after trimming.

---

### **Step 4 — Align Reads to the Reference Genome**

Before mapping, build a **HISAT2 index** of your reference genome (only needed once):

```bash
hisat2-build reference_genome.fa genome
```

Then align reads for each sample:

```bash
hisat2 -q -x genome \
  -1 read_file_1.fastq \
  -2 read_file_2.50k.fastq \
  -S control_rep1.sam \
  --summary-file summary_control_rep1.txt
```

**Explanation**

* `-x genome` → HISAT2 index basename.
* `-1` and `-2` → paired-end input reads.
* `-S` → SAM alignment file output.
* `--summary-file` → alignment summary with mapping percentages.

---

### **Step 5 — Convert and Sort Alignments**

```bash
samtools view -S -b control_rep1.sam > control_rep1.bam
samtools sort -o control_rep1_sorted.bam control_rep1.bam
```

**Explanation**

* Converts SAM to compressed BAM format.
* Sorts reads by genomic coordinates (required for counting and visualization).

---

### **Step 6 — Generate Read Counts**

```bash
htseq-count -f bam control_rep1_sorted.bam tomato.gtf > count_control_rep1.counts
```

**Explanation**

* `-f bam` → input format is BAM.
* `tomato.gtf` → genome annotation file with gene coordinates.
* Output file (`count_control_rep1.counts`) contains per-gene read counts.

---

## Full Workflow Summary

```bash
# 1. Quality check (raw reads)
fastqc control_rep1_1.50k.fastq control_rep1_2.50k.fastq -o fastqc_raw

# 2. Trim 5′ and 3′ ends
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
  control_rep1_1.50k.fastq control_rep1_2.50k.fastq -o trimmed

# 3. QC trimmed reads
fastqc trimmed/control_rep1_1.50k_val_1.fq trimmed/control_rep1_2.50k_val_2.fq -o fastqc_trimmed

# 4. Align reads
hisat2 -q -x genome \
  -1 control_rep1_1.50k.fastq -2 control_rep1_2.50k.fastq \
  -S control_rep1.sam --summary-file summary_control_rep1.txt

# 5. Convert and sort
samtools view -S -b control_rep1.sam > control_rep1.bam
samtools sort -o control_rep1_sorted.bam control_rep1.bam

# 6. Count reads per gene
htseq-count -f bam control_rep1_sorted.bam tomato.gtf > count_control_rep1.counts
```

Repeat the above commands for:

* **control_rep2**, **control_rep3**
* **treated_rep1**, **treated_rep2**, **treated_rep3**

---

## Output Overview

| Step | Tool        | Input         | Output                     |
| ---- | ----------- | ------------- | -------------------------- |
| 1    | FastQC      | Raw FASTQ     | HTML + ZIP quality reports |
| 2    | Trim Galore | Raw FASTQ     | Trimmed FASTQ + reports    |
| 3    | FastQC      | Trimmed FASTQ | QC reports after trimming  |
| 4    | HISAT2      | Trimmed FASTQ | SAM alignment              |
| 5    | Samtools    | SAM           | Sorted BAM                 |
| 6    | HTSeq       | BAM + GTF     | Gene-level read counts     |

---

## Notes for Students

* Process **each replicate** separately (both control and treated).
* Check FastQC reports after trimming — all quality modules should pass.
* `.counts` files are used as input for **DESeq2**, **edgeR**, or **limma-voom** for differential expression.
* Ensure the **reference genome FASTA** and **annotation GTF** are from the same version.
* Run HISAT2, Samtools, and HTSeq in the same working directory.

---

## Optional: Run the Pipeline for All Samples

Automate processing of all paired-end samples:

```bash
for s in *_1.50k.fastq; do
  base=${s%_1.50k.fastq}
  fastqc ${base}_1.50k.fastq ${base}_2.50k.fastq -o fastqc_raw
  trim_galore --paired \
    --clip_R1 10 --clip_R2 10 \
    --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
    ${base}_1.50k.fastq ${base}_2.50k.fastq -o trimmed
  hisat2 -q -x genome \
    -1 ${base}_1.50k.fastq -2 ${base}_2.50k.fastq \
    -S ${base}.sam --summary-file summary_${base}.txt
  samtools view -S -b ${base}.sam > ${base}.bam
  samtools sort -o ${base}_sorted.bam ${base}.bam
  htseq-count -f bam ${base}_sorted.bam tomato.gtf > count_${base}.counts
done
```

---

## Author

Developed by **effective-robot**
for **hands-on RNA-seq training workshops** on read preprocessing, alignment, and quantification.

---

### End of README

```
