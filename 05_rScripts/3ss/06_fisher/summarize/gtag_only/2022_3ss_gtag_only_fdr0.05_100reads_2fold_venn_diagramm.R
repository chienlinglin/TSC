# Load library
library(VennDiagram)
library(dplyr)

########## filter 100 read count
# input datasets 
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
su <- read.csv(file = "spliced_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_spliced_unspliced_sort_fvals.csv",row.names = 1)
cn <- read.csv(file = "can_nocan/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu <- read.csv(file = "can_noncan_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)

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
# setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/summarize/gtag_only/")
setwd("~/Desktop/")

venn_diagramm_3ss <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("spliced_unspliced" , "canonical_noncanonical" , "canonical_noncan+un"),
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
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-20, 6, 180),
  cat.dist = c(0.025, 0.06, 0.01),
  cat.fontfamily = "sans",
  rotation = 1
)

pdf(file="2022_3ss_significant_gtag_only_fdr0.05_100reads_2fold_venn_diagramm.pdf")
grid.draw(venn_diagramm_3ss)
dev.off()

#####
# set1(spliced_unspliced): 107
# set2(can_noncan): 655
# set3(can_noncan+un): 554
### total: 107+655+554-(2+78)-(78+435)-(78+9)+78=714
### 2386 wt-mt pair
### 714/2386= 29.9%
