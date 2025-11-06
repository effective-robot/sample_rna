```markdown
# 🧬 RNA-seq Practice Workflow (Paired-End)

This repository provides a **hands-on practice pipeline** for students learning RNA-seq data analysis.  
It uses **publicly available sample reads** (uploaded here), and covers every major preprocessing step — from quality control to read counting.

---
```
## 📂 Repository Contents

```

sample/
├── 50k_subsamples/              # Paired-end FASTQ files (sample data)
├── trim_and_qc_one_sample.sh    # Automated trimming + QC script (optional)
└── README.md                    # This guide

````

---

## 🚀 Quick Start

### 1️⃣ Install Required Tools

Make sure these tools are installed and available in your PATH:

| Tool | Purpose | Install command (conda) |
|------|----------|--------------------------|
| **FastQC** | Quality check of reads | `conda install -c bioconda fastqc` |
| **Trim Galore** | Adapter and quality trimming | `conda install -c bioconda trim-galore cutadapt` |
| **HISAT2** | RNA-seq read alignment | `conda install -c bioconda hisat2` |
| **Samtools** | BAM/SAM processing | `conda install -c bioconda samtools` |
| **HTSeq** | Gene-level read counting | `conda install -c bioconda htseq` |

Alternatively, if Conda is not available, you can use:

```bash
sudo apt install fastqc hisat2 samtools htseq
python3 -m pip install --user cutadapt
````

---

## 🧪 Example Workflow

Below is a complete example using one sample (e.g., **SRR866245**).
Repeat for each sample pair in the dataset.

---

### **Step 1 — Quality Check (Raw Reads)**

```bash
fastqc SRR866245_1.fastq SRR866245_2.fastq -o fastqc_raw
```

**Explanation**

* `SRR866245_1.fastq` and `SRR866245_2.fastq` are forward and reverse reads
* `-o fastqc_raw` specifies an output directory for the FastQC reports
* Open `fastqc_raw/SRR866245_1_fastqc.html` in a web browser to view quality metrics

---

### **Step 2 — Trim Reads**

```bash
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
  SRR866245_1.fastq SRR866245_2.fastq \
  -o trimmed_reads
```

**Explanation**

* `--paired` → indicates paired-end data
* `--clip_R1/2 10` → trims 10 bases from the 5′ end of each read
* `--three_prime_clip_R1/2 5` → trims 5 bases from the 3′ end of each read
* `-o trimmed_reads` → output directory for trimmed reads and trimming reports

**Output Files**

```
trimmed_reads/
 ├── SRR866245_1.fastq_trimming_report.txt
 ├── SRR866245_2.fastq_trimming_report.txt
 ├── SRR866245_1_val_1.fq
 └── SRR866245_2_val_2.fq
```

---

### **Step 3 — Re-run QC on Trimmed Reads**

```bash
fastqc trimmed_reads/SRR866245_1_val_1.fq trimmed_reads/SRR866245_2_val_2.fq -o fastqc_trimmed
```

This checks whether the **“Per base sequence content”** and other metrics improved after trimming.

---

### **Step 4 — Align Reads to the Reference Genome**

Before mapping, make sure you have a **HISAT2 index** of your reference genome (e.g., `genome.*.ht2` files).
If not yet built, create one using:

```bash
hisat2-build reference_genome.fa genome
```

Then align reads:

```bash
hisat2 -q -x genome \
  -1 SRR866245_1.fastq \
  -2 SRR866245_2.fastq \
  -S SRR866245.sam \
  --summary-file summary_SRR866245.txt
```

**Explanation**

* `-x genome` → HISAT2 index basename
* `-1` and `-2` → input paired reads
* `-S` → output SAM file
* `--summary-file` → saves alignment summary statistics

---

### **Step 5 — Convert and Sort Alignments**

```bash
samtools view -S -b SRR866245.sam > SRR866245.bam
samtools sort -o SRR866245_sorted.bam SRR866245.bam
```

**Explanation**

* Converts SAM to BAM format
* Sorts BAM file by genomic coordinates (required for downstream tools)

---

### **Step 6 — Generate Read Counts**

```bash
htseq-count -f bam SRR866245_sorted.bam tomato.gtf > count_SRR866245.counts
```

**Explanation**

* `-f bam` → input format is BAM
* `tomato.gtf` → annotation file containing gene coordinates
* Output file (`count_SRR866245.counts`) contains gene-wise read counts

---

## 🧭 Full Workflow Summary

```bash
# 1. Quality check (raw reads)
fastqc SRR866245_1.fastq SRR866245_2.fastq -o fastqc_raw

# 2. Trim 5′ and 3′ ends
trim_galore --paired \
  --clip_R1 10 --clip_R2 10 \
  --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
  SRR866245_1.fastq SRR866245_2.fastq -o trimmed

# 3. QC on trimmed reads
fastqc trimmed/SRR866245_1_val_1.fq trimmed/SRR866245_2_val_2.fq -o fastqc_trimmed

# 4. Align reads
hisat2 -q -x genome \
  -1 SRR866245_1.fastq -2 SRR866245_2.fastq \
  -S SRR866245.sam --summary-file summary_SRR866245.txt

# 5. Convert and sort
samtools view -S -b SRR866245.sam > SRR866245.bam
samtools sort -o SRR866245_sorted.bam SRR866245.bam

# 6. Count reads per gene
htseq-count -f bam SRR866245_sorted.bam tomato.gtf > count_SRR866245.counts
```

---

## 📊 Output Overview

| Step | Tool        | Input         | Output                     |
| ---- | ----------- | ------------- | -------------------------- |
| 1    | FastQC      | Raw FASTQ     | HTML + ZIP quality reports |
| 2    | Trim Galore | Raw FASTQ     | Trimmed FASTQ + reports    |
| 3    | FastQC      | Trimmed FASTQ | QC reports after trimming  |
| 4    | HISAT2      | Trimmed FASTQ | SAM alignment              |
| 5    | Samtools    | SAM           | Sorted BAM                 |
| 6    | HTSeq       | BAM + GTF     | Gene-level read counts     |

---

## 💡 Notes for Students

* Each sample must be processed separately (replace `SRR866245` with your sample ID).
* Verify that all FastQC quality modules pass after trimming.
* The `.counts` files are suitable for differential expression tools such as **DESeq2**, **edgeR**, or **limma-voom**.
* Ensure that the genome FASTA and annotation GTF files come from the same genome version.
* HISAT2, Samtools, and HTSeq should all run in the same working directory.

---

## 🧰 Optional: Run the Pipeline for All Samples

You can automate all paired-end samples with this simple loop:

```bash
for s in *_1.fastq; do
  base=${s%_1.fastq}
  fastqc ${base}_1.fastq ${base}_2.fastq -o fastqc_raw
  trim_galore --paired \
    --clip_R1 10 --clip_R2 10 \
    --three_prime_clip_R1 5 --three_prime_clip_R2 5 \
    ${base}_1.fastq ${base}_2.fastq -o trimmed
  hisat2 -q -x genome \
    -1 ${base}_1.fastq -2 ${base}_2.fastq \
    -S ${base}.sam --summary-file summary_${base}.txt
  samtools view -S -b ${base}.sam > ${base}.bam
  samtools sort -o ${base}_sorted.bam ${base}.bam
  htseq-count -f bam ${base}_sorted.bam tomato.gtf > count_${base}.counts
done
```

---

## 👨‍🏫 Author

Developed by **effective-robot**
for hands-on **RNA-seq training workshops** on data preprocessing, alignment, and quantification.

---

### ✅ End of README

```
