# 5'SS library Analysis Pipeline


### Phase I: Data Preprocessing

#### 1. Quality Control (`01_2022_5ss_total__hist_pairplot.R`)
- **Input**: Raw read count files
- **Process**: 
  - Merge technical replicates
  - Generate count distribution histograms
  - Calculate inter-sample correlations
- **Output**: QC reports and correlation plots

#### 2. Splice Site Annotation (`02_2022_5ss_w:wo_mt_annotation_sum.R`)
- **Function**: Identify and classify splice site types
- **Process**:
  - Mark GT-AG/CT-AC as canonical
  - Label other sequences as non-canonical
  - Remove barcode interference regions
  - Process variant and wild-type data separately
- **Output**: Standardized splice site annotation files

### Phase II: Splicing Efficiency Analysis

#### 3. Canonical Splicing Analysis (`03_2022_5ss_w:wo_mt_annotation_gt_ag_only.R`)
- **Target**: Focus on GT-AG canonical splice sites
- **Process**: Filter to retain only GT-AG splicing events
- **Output**: GT-AG specific datasets

#### 4. Optimal Splice Site Identification (`04_2022_5ss_build_add_GT_bed.R`)
- **Function**: Determine primary splice sites for each oligo
- **Method**: 
  - Identify highest-expressed splice site per sequence
  - Account for exon/intron variant effects on splice site positions
- **Output**: Primary splice site coordinate files

### Phase III: Statistical Testing

#### 5. Data Integration (`05_2022_count_5ss_gt_ag_only_integration_for_fisher.R`)
- **Function**: Integrate multiple data sources for statistical analysis
- **Integration**:
  - Spliced counts
  - Unspliced counts
  - Canonical vs non-canonical classification
  - GT addition effect analysis
- **Strategies**:
  - **"Correct canonical"**: Use predefined canonical sites
  - **"Most canonical"**: Use highest-expressed sites as canonical

#### 6. Statistical Comparisons

##### Canonical vs Non-canonical (`2022_run3_run4_5ss_gtag_only_can_noncan_fisher.R`)
```r
# Statistical methods: Chi-square test & Fisher's exact test
# Comparison: Canonical splicing vs Non-canonical splicing
# Multiple testing correction: FDR control
```

##### Canonical vs (Non-canonical + Unspliced) (`2022_run3_run4_5ss_gtag_only_can_noncan_un_fisher.R`)
```r
# Extended analysis: Include unspliced in non-canonical category
# Purpose: More sensitive detection of splicing efficiency changes
```

##### Spliced vs Unspliced (`2022_run3_run4_5ss_gtag_only_spliced_unspliced_fisher.R`)
```r
# Core analysis: Direct splicing efficiency comparison
# Method: Fisher's exact test
# Metric: Spliced/unspliced ratio changes
```

### Phase IV: Candidate Variant Selection

#### 7. Comprehensive Filtering (`2022_5ss_candidate_info.R`)

##### Statistical Significance Criteria
- Bonferroni-corrected p < 0.05/2644
- Significant across all three replicates

##### Biological Significance Criteria
- Read count ≥ 100
- Odds ratio ≥ 2 or ≤ 0.5
- FDR ≤ 0.05

## 📊 Output Files

### Primary Results
- `candidate_variants_info.csv`: Comprehensive candidate variant information
- `statistical_test_summary.csv`: Statistical testing results summary
- `functional_annotation_table.csv`: Integrated functional annotations
- `clinical_relevance_report.csv`: Clinical significance assessment

### Quality Control
- `sample_correlation_plots.pdf`: Inter-sample correlation analysis
- `count_distribution_histograms.pdf`: Read count distributions
- `splicing_efficiency_plots.pdf`: Splicing efficiency visualizations

### Intermediate Files
- `normalized_splice_counts/`: Processed splice junction counts
- `statistical_tests/`: Individual statistical test results
- `annotation_files/`: Splice site annotations and classifications


**Last Updated**: January 2025
