# Load library
library(VennDiagram)
library(dplyr)

########## filter 100 read count
# input datasets 
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
su <- read.csv(file = "spliced_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_spliced_unspliced_sort_fvals.csv",row.names = 1)
cn <- read.csv(file = "can_nocan/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu <- read.csv(file = "can_noncan_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)
## preparing for |wt-mt|>30
wt_s_e <- na.omit(cnu[,c(2,3,8,10,18,20,28,30)])
wt_se <- 100*((wt_s_e$wt_can)/(wt_s_e$wt_sum))
wt_se1 <- 100*((wt_s_e$wt_can.1)/(wt_s_e$wt_sum.1))
wt_se2 <- 100*((wt_s_e$wt_can.2)/(wt_s_e$wt_sum.2))
wt_s_e <- cbind(wt_se,wt_se1,wt_se2,wt_s_e)
colnames(wt_s_e)[1] <- "wt_splicing_efficiency"
colnames(wt_s_e)[2] <- "wt_splicing_efficiency1"
colnames(wt_s_e)[3] <- "wt_splicing_efficiency2"
wt_s_e[is.na(wt_s_e)] <- 0

mt_s_e <- na.omit(cnu[,c(3,11,13,21,23,31,33)])
mt_se <- 100*((mt_s_e$mt_can)/(mt_s_e$mt_sum))
mt_se1 <- 100*((mt_s_e$mt_can.1)/(mt_s_e$mt_sum.1))
mt_se2 <- 100*((mt_s_e$mt_can.2)/(mt_s_e$mt_sum.2))
mt_s_e <- cbind(mt_se,mt_se1,mt_se2,mt_s_e)
colnames(mt_s_e)[1] <- "mt_splicing_efficiency"
colnames(mt_s_e)[2] <- "mt_splicing_efficiency1"
colnames(mt_s_e)[3] <- "mt_splicing_efficiency2"
mt_s_e[is.na(mt_s_e)] <- 0

se_difference <- abs(wt_s_e[,c(1)]-mt_s_e[,c(1)])
se_difference1 <- abs(wt_s_e[,c(2)]-mt_s_e[,c(2)])
se_difference2 <- abs(wt_s_e[,c(3)]-mt_s_e[,c(3)])
se_difference_1 <- cbind(se_difference,se_difference1,se_difference2,mt_s_e)


se_difference_1$se_difference_3 <- (se_difference_1[,c(1)]>=30) + (se_difference_1[,c(2)]>=30) + (se_difference_1[,c(3)]>=30)



##########################
qval_3 <- (su$FDR_qval<=0.05) + (su$FDR_qval.1<=0.05)+(su$FDR_qval.2<=0.05)
wt_sum_3 <- (su$wt_sum>=100) + (su$wt_sum.1>=100)+(su$wt_sum.2>=100)
mt_sum_3 <- (su$mt_sum>=100) + (su$mt_sum.1>=100)+(su$mt_sum.2>=100)

su_f <- cbind(qval_3,wt_sum_3,mt_sum_3,su)
su_f <- merge(su_f,se_difference_1[,c(7,14)],by.x ="mt",by.y="mt")

qval_3 <- (cn$FDR_qval<=0.05) + (cn$FDR_qval.1<=0.05)+(cn$FDR_qval.2<=0.05)
wt_sum_3 <- (cn$wt_sum>=100) + (cn$wt_sum.1>=100)+(cn$wt_sum.2>=100)
mt_sum_3 <- (cn$mt_sum>=100) + (cn$mt_sum.1>=100)+(cn$mt_sum.2>=100)

cn_f <- cbind(qval_3,wt_sum_3,mt_sum_3,cn)
cn_f <- merge(cn_f,se_difference_1[,c(7,14)],by.x ="mt",by.y="mt")

qval_3 <- (cnu$FDR_qval<=0.05) + (cnu$FDR_qval.1<=0.05)+(cnu$FDR_qval.2<=0.05)
wt_sum_3 <- (cnu$wt_sum>=100) + (cnu$wt_sum.1>=100)+(cnu$wt_sum.2>=100)
mt_sum_3 <- (cnu$mt_sum>=100) + (cnu$mt_sum.1>=100)+(cnu$mt_sum.2>=100)

cnu_f <- cbind(qval_3,wt_sum_3,mt_sum_3,cnu)
cnu_f <- merge(cnu_f,se_difference_1[,c(7,14)],by.x ="mt",by.y="mt")

# Generate 3 sets of 200 words
set <- su_f %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% select(mt) %>% unlist() #1380
set1_ <- su_f %>% filter(qval_3 == 3) %>% filter(se_difference_3 == 3) %>% select(mt) %>% unlist() #71
set1 <- intersect(set1_,set) #spliced_unspliced = 46
set2_ <- cn_f %>% filter(qval_3 == 3) %>% filter(se_difference_3 == 3) %>% select(mt) %>% unlist() #286
set2 <- intersect(set2_,set) #canonical_nocanonical = 163
set3_ <- cnu_f %>% filter(qval_3 == 3) %>% filter(se_difference_3 == 3) %>% select(mt) %>% unlist() #295
set3 <- intersect(set3_,set) #cn+uns = 163

# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")

# Chart
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/summarize/gtag_only/")

venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("spliced_unspliced" , "canonical_noncanonical" , "canonical_noncan+un"),
  filename = '2022_3ss_significant_gtag_only_30_fdr0.05_100reads_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 900 , 
  width = 900 , 
  resolution = 400,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .8,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(2, 2, 180),
  cat.dist = c(-0.12, -0.01, 0.05),
  cat.fontfamily = "sans",
  rotation = 3
)

#####
# set1(spliced_unspliced): 107
# set2(can_noncan): 655
# set3(can_noncan+un): 554
### total: 107+655+554-(2+78)-(78+435)-(78+9)+78=714
### 2386 wt-mt pair
### 714/2386= 29.9%
