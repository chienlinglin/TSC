## 🔬 Analysis Pipeline Overview

This pipeline processes FASTQ files through alignment, filtering, and read classification steps, producing various BAM and count outputs.


---

### 🧩 Pipeline Steps Summary

1. **alignment.sh**
   - Aligns raw FASTQ reads to a reference genome using HISAT2.
   - Outputs:
     - Aligned SAM file
     - Unaligned FASTQ
     - Mapping statistics log

2. **bamcount.sh**
   - Filters reads with mapping quality ≥ Q60 using `samtools`
   - Extracts paired-end reads
   - Splits BAM by mitochondrial info (`mt_info` vs `no_mt_info`)

3. **Splicing Analysis**
   - Classifies reads as spliced or unspliced from `no_mt_info` BAM
   - Generates count tables using `bamcount.sh`

4. **Junction Analysis**
   - Extracts splice junctions using `regtools`
   - Outputs junction BED file with functional annotations

---

### 🛠️ Software & Tool Versions

| Tool             | Version     |
|------------------|-------------|
| `Trimmomatic`    | 0.39        |
| `HISAT2`         | 2.2.1       |
| `samtools`       | 1.15.1      |
| `bedtools`       | 2.30.0      |
| `regtools`       | 0.5.2       |
| `FASTX-Toolkit`  | 0.0.14      |

---



