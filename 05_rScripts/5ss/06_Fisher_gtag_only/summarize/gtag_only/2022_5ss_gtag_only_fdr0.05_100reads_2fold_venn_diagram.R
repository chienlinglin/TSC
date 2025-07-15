# Load library
library(VennDiagram)
library(dplyr)

########## filter 100 read count
# input datasets 
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/")
su <- read.csv(file = "spliced_unspliced/gtag_only/2022_run3_run4_5ss_gtag_only_spliced_unspliced_sort_fvals.csv",row.names = 1)
cn <- read.csv(file = "can_noncan/gtag_only/2022_run3_run4_5ss_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu <- read.csv(file = "can_noncan_unspliced/gtag_only/2022_run3_run4_5ss_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)


qval_3 <- (su$FDR_qval<=0.05) + (su$FDR_qval.1<=0.05)+(su$FDR_qval.2<=0.05)
wt_sum_3 <- (su$wt_sum>=100) + (su$wt_sum.1>=100)+(su$wt_sum.2>=100)
mt_sum_3 <- (su$mt_sum>=100) + (su$mt_sum.1>=100)+(su$mt_sum.2>=100)
or_3 <- (su$odds_ratio>=2 | su$odds_ratio<=0.5) + (su$odds_ratio.1>=2 | su$odds_ratio.1<=0.5)+(su$odds_ratio.2>=2 | su$odds_ratio.2<=0.5)
su_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,su)


qval_3 <- (cn$FDR_qval<=0.05) + (cn$FDR_qval.1<=0.05)+(cn$FDR_qval.2<=0.05)
wt_sum_3 <- (cn$wt_sum>=100) + (cn$wt_sum.1>=100)+(cn$wt_sum.2>=100)
mt_sum_3 <- (cn$mt_sum>=100) + (cn$mt_sum.1>=100)+(cn$mt_sum.2>=100)
or_3 <- (cn$odds_ratio>=2 | cn$odds_ratio<=0.5) + (cn$odds_ratio.1>=2 | cn$odds_ratio.1<=0.5)+(cn$odds_ratio.2>=2 | cn$odds_ratio.2<=0.5)
cn_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,cn)


qval_3 <- (cnu$FDR_qval<=0.05) + (cnu$FDR_qval.1<=0.05)+(cnu$FDR_qval.2<=0.05)
wt_sum_3 <- (cnu$wt_sum>=100) + (cnu$wt_sum.1>=100)+(cnu$wt_sum.2>=100)
mt_sum_3 <- (cnu$mt_sum>=100) + (cnu$mt_sum.1>=100)+(cnu$mt_sum.2>=100)
or_3 <- (cnu$odds_ratio>=2 | cnu$odds_ratio<=0.5) + (cnu$odds_ratio.1>=2 | cnu$odds_ratio.1<=0.5)+(cnu$odds_ratio.2>=2 | cnu$odds_ratio.2<=0.5)
cnu_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,cnu)


# Generate 3 sets of 200 words
set <- su_f %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% select(mt) %>% unlist()
set1_ <- su_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set1 <- intersect(set1_,set)
set2_ <- cn_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set2 <- intersect(set2_,set)
set3_ <- cnu_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set3 <- intersect(set3_,set)

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
# setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/summarize/gtag_only")
setwd("/Users/angchu/Desktop/phD/TSC_codes/2025_plots/5ss")
venn_diagramm_5ss <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Unspliced" , "Noncanonical" , "Noncanonical+\nUnspliced"),
  filename = NULL,
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 600 , 
  width = 600 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-16, 5, 180),
  cat.dist = c(0.005, 0.003, 0.01),
  cat.fontfamily = "sans",
  rotation = 1
)
pdf(file="2025_5ss_significant_gtag_only_fdr0.05_100reads_2fold_venn_diagramm.pdf")
grid.draw(venn_diagramm_5ss)
dev.off()

#####
# set1(spliced_unspliced): 41 -- >39
# set2(can_noncan): 503 --> 488
# set3(can_noncan+un): 528 --> 463
### total: 597  --> 526
### 2644 wt-mt pair
### 597/2644= 22.58% -- > 526/2644 = 20%
