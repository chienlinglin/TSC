# 3'SS Library Analysis Pipeline

### Phase I: Data Preprocessing

#### 1. Quality Control (`01_2022_3ss_total__hist_pairplot.R`)
- **Input**: Raw read count files
- **Process**: 
  - Merge technical replicates (3-1, 3-2, 3-3)
  - Generate count distribution histograms (2525 oligos)
  - Calculate inter-sample correlations using Spearman correlation
- **Output**: QC reports and correlation plots

#### 2. Splice Site Annotation (`02_2022_3ss_w:wo_mt_annotation_sum.R`)
- **Function**: Identify and classify splice site types
- **Process**:
  - Mark GT-AG/CT-AC as canonical
  - Label other sequences as non-canonical
  - Remove barcode interference regions (3'ss specific filtering)
  - Process variant and wild-type data separately
- **Output**: Standardized splice site annotation files

### Phase II: Splicing Efficiency Analysis

#### 3. Canonical Splicing Analysis (`03_2022_3ss_w:wo_mt_annotation_gt_ag_only_sum.R`)
- **Target**: Focus on GT-AG canonical splice sites
- **Process**: Filter to retain only GT-AG splicing events
- **Output**: GT-AG specific datasets for 3'ss analysis

#### 4. Optimal Splice Site Identification (`04_2022_3ss_build_add_AG_bed.R`)
- **Function**: Determine primary splice sites for each oligo
- **Method**: 
  - Identify highest-expressed splice site per sequence
  - Calculate AG dinucleotide positions for 3'ss
  - Account for exon/intron variant effects on splice site positions
- **Output**: Primary splice site coordinate files (AG position)

### Phase III: Statistical Testing

#### 5. Data Integration (`06_2025_3ss_gt_ag_only_integration_for_fisher.R`)
- **Function**: Comprehensive data integration for statistical analysis
- **Integration**:
  - Spliced counts
  - Unspliced counts
  - Canonical vs non-canonical classification
  - AG addition effect analysis
- **Strategies**:
  - **"Correct canonical"**: Use predefined canonical sites
  - **"Most canonical"**: Use highest-expressed sites as canonical

#### 6. Statistical Comparisons

##### Canonical vs Non-canonical (`2022_run3_run4_3ss_gtag_only_can_noncan_fisher.R`)
```r
# Statistical methods: Chi-square test & Fisher's exact test
# Comparison: Canonical splicing vs Non-canonical splicing
# Multiple testing correction: FDR control
# Sample size: 2386 variant pairs
```

##### Canonical vs (Non-canonical + Unspliced) (`2022_run3_run4_3ss_gtag_only_0525_0711_can_noncan_un_fisher.R`)
```r
# Extended analysis: Include unspliced in non-canonical category
# Purpose: More sensitive detection of splicing efficiency changes
# Sample size: 2386 variant pairs
```

##### Spliced vs Unspliced (`2022_run3_run4_3ss_gtag_only_spliced_unspliced_fisher.R`)
```r
# Core analysis: Direct splicing efficiency comparison
# Method: Fisher's exact test
# Metric: Spliced/unspliced ratio changes
# Sample size: 2386 variant pairs
```

### Phase IV: Candidate Variant Selection

#### 7. Comprehensive Filtering (`2022_3ss_candidate_info.R`)

##### Statistical Significance Criteria
- Bonferroni-corrected p < 0.05/2386
- Significant across all three replicates

##### Biological Significance Criteria
- Read count ≥ 100
- Odds ratio ≥ 2 or ≤ 0.5
- FDR ≤ 0.05

##### AG Dinucleotide Effect Analysis
- AG addition/reduction effect assessment
- AG-specific statistical significance testing

**Last Updated**: July 2025
