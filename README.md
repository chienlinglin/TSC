## Analysis Pipeline Overview

This pipeline processes FASTQ files through alignment, filtering, and read classification steps, producing various BAM and count outputs.

Separate shell scripts are provided for the 3′ splice-site (3′SS) and 5′ splice-site (5′SS) libraries. 
The overall analysis workflow is the same for both libraries, but the scripts use library-specific input paths, reference sequences/indexes, and downstream analysis files. 
Therefore, 3′SS and 5′SS scripts should be applied only to data generated from their corresponding libraries.

---

### Pipeline Steps Summary

1. **Trimmomatic** (located in `01_trimmomatic/`)
   - Scripts: `20191122_3ss_trimmomatic.sh` and `20191122_5ss_trimmomatic.sh`
   - Performs quality control and adapter trimming on paired-end FASTQ files using `Trimmomatic`.
   - Removes Nextera adapter sequences and low-quality bases.

2. **Alignment** (located in `02_alignment/`)
   - Scripts: `2022_angchu_3ss_alignment.sh` and `2022_angchu_5ss_alignment.sh`
   - Aligns raw FASTQ reads to a reference genome using `HISAT2`.

3. **BAM Filtering and Read Classification** (located in `03_bedtools_bamcount/`)
   - Scripts: `2022_3ss_bedtools_bamcount.sh` and `2022_5ss_bedtools_bamcount.sh`
   - Filters reads with mapping quality ≥ Q60 using `samtools`
   - Extracts paired-end reads
   - Splits BAM reads by variant position (`mt_info` vs `no_mt_info`)
   - Classifies reads as spliced or unspliced from BAM files
   - Generates count tables

4. **Junction Analysis** (located in `04_regtools_juncbed/`)
   - Scripts: `2022_3ss_regtools_juncbed.sh` and `2022_5ss_regtools_juncbed.sh`
   - Extracts splice junctions using `regtools`
   - Outputs junction BED file with functional annotations

5. **R Analysis Scripts** (located in `05_rScripts/`)
   - Performs downstream statistical analyses and visualization for the 3′SS and 5′SS datasets.
   - Detailed descriptions and instructions for the R scripts are provided in the README files within the corresponding 3′SS and 5′SS subdirectories.

---

### Running the shell scripts

The shell scripts use predefined project, input, and output directories specified at the beginning of each script. 
Before running a script, users should modify the path variables according to their local directory structure. For example:

```bash
projectPath=/path/to/project/
readPath=${projectPath}3ss_raw/
outputPath=${projectPath}3ss_qc/
```
Input filenames are not passed explicitly on the command line. Each script automatically identifies input files from the configured input directory using predefined filename patterns.
Run the scripts using:
```bash
bash <script_name>.sh
```

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



