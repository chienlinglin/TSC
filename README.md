
---

### 🧩 Pipeline Steps Summary

1. **alignment.sh**
   - Aligns raw FASTQ reads to a reference genome.
   - Outputs:
     - Aligned SAM file
     - Unaligned FASTQ
     - Mapping statistics log

2. **bamcount.sh**
   - Filters reads with mapping quality ≥ Q60
   - Extracts paired-end reads
   - Splits BAM by mitochondrial info (`mt_info` vs `no_mt_info`)

3. **Splicing Analysis**
   - Classifies reads as spliced or unspliced from `no_mt_info` BAM
   - Generates count tables using `bamcount.sh`

4. **Junction Analysis**
   - Extracts splice junctions using `juncbed.sh`
   - Outputs junction BED file with functional annotations

---

