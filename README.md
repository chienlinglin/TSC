## 🔬 Analysis Pipeline Overview

This pipeline processes FASTQ files through alignment, filtering, and read classification steps, producing various BAM and count outputs.


---

### 🧩 Pipeline Steps Summary

0. **Trimmomatic** (located in `01_trimmomatic/`)
   - Performs quality control and adapter trimming on paired-end FASTQ files using `Trimmomatic`.
   - Removes Nextera adapter sequences and low-quality bases.

1. **alignment.sh** (located in `02_alignment/`)
   - Aligns raw FASTQ reads to a reference genome using `HISAT2`.

2. **bamcount.sh** (located in `03_bedtools_bamcount/`)
   - Filters reads with mapping quality ≥ Q60 using `samtools`
   - Extracts paired-end reads
   - Splits BAM reads by variant position (`mt_info` vs `no_mt_info`)
   - Classifies reads as spliced or unspliced from BAM files
   - Generates count tables

3. **Junction Analysis** (located in `04_regtools_juncbed/`)
   - Extracts splice junctions using `regtools`
   - Outputs junction BED file with functional annotations

3. **R Analysis Scripts** (located in `05_rScripts/`)
   - Performs downstream statistical analyses and visualization for the 3′SS and 5′SS datasets.
   - Detailed descriptions and instructions for the R scripts are provided in the README files within the corresponding 3′SS and 5′SS subdirectories.

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



