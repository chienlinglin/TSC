# Load library
library(dplyr)
library(stringr)
library("Biostrings")
library(readxl)
library(tidyr)

########## filter 100 read count
# input datasets 
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/")
su <- read.csv(file = "spliced_unspliced/gtag_only/2022_run3_run4_5ss_gtag_only_spliced_unspliced_sort_fvals.csv",row.names = 1)
cn <- read.csv(file = "can_noncan/gtag_only/2022_run3_run4_5ss_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu <- read.csv(file = "can_noncan_unspliced/gtag_only/2022_run3_run4_5ss_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)

## filter 100 read count
#use cnu_f
wt_sum_3 <- (cnu$wt_sum>=100) + (cnu$wt_sum.1>=100)+(cnu$wt_sum.2>=100)
mt_sum_3 <- (cnu$mt_sum>=100) + (cnu$mt_sum.1>=100)+(cnu$mt_sum.2>=100)
cnu_f <- cbind(wt_sum_3,mt_sum_3,cnu)
cnu_f <- su_f[,c(1,2,4,5,12,15,22,25,32,35)]


set1 <- su_f_reads %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% select(wt_sum_3,mt_sum_3,mt,wt)

colnames(set1)[1] <-"wt_reads_over_100"
colnames(set1)[2] <-"mt_reads_over_100"
setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/summary/")
write.csv(set1,"5ss_candidates_reads_100.csv",col.names = NULL)

### use pvalue only
set1p <- su %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set2p <- cn %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set3p <- cnu %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()


### use pvalue plus reads 100
su$wt_sum_3 <- (su$wt_sum>=100) + (su$wt_sum.1>=100)+(su$wt_sum.2>=100)
su$mt_sum_3 <- (su$mt_sum>=100) + (su$mt_sum.1>=100)+(su$mt_sum.2>=100)
set1p_r <- su %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% filter(pass_3_fval_strict == 3) %>% unlist()

cn$wt_sum_3 <- (cn$wt_sum>=100) + (cn$wt_sum.1>=100)+(cn$wt_sum.2>=100)
cn$mt_sum_3 <- (cn$mt_sum>=100) + (cn$mt_sum.1>=100)+(cn$mt_sum.2>=100)
set2p_r <- cn %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()

cnu$wt_sum_3 <- (cnu$wt_sum>=100) + (cnu$wt_sum.1>=100)+(cnu$wt_sum.2>=100)
cnu$mt_sum_3 <- (cnu$mt_sum>=100) + (cnu$mt_sum.1>=100)+(cnu$mt_sum.2>=100)
set3p_r <- cnu %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()



### use hunglun's condition: fdr 0.05, reads 100, odds ratio 2 fold
qval_3 <- (su$FDR_qval<=0.05) + (su$FDR_qval.1<=0.05)+(su$FDR_qval.2<=0.05)
#wt_sum_3 <- (su$wt_sum>=100) + (su$wt_sum.1>=100)+(su$wt_sum.2>=100)
#mt_sum_3 <- (su$mt_sum>=100) + (su$mt_sum.1>=100)+(su$mt_sum.2>=100)
or_3 <- (su$odds_ratio>=2 | su$odds_ratio<=0.5) + (su$odds_ratio.1>=2 | su$odds_ratio.1<=0.5)+(su$odds_ratio.2>=2 | su$odds_ratio.2<=0.5)
#su_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,su)
su_f <- cbind(qval_3,or_3,su)

qval_3 <- (cn$FDR_qval<=0.05) + (cn$FDR_qval.1<=0.05)+(cn$FDR_qval.2<=0.05)
#wt_sum_3 <- (cn$wt_sum>=100) + (cn$wt_sum.1>=100)+(cn$wt_sum.2>=100)
#mt_sum_3 <- (cn$mt_sum>=100) + (cn$mt_sum.1>=100)+(cn$mt_sum.2>=100)
or_3 <- (cn$odds_ratio>=2 | cn$odds_ratio<=0.5) + (cn$odds_ratio.1>=2 | cn$odds_ratio.1<=0.5)+(cn$odds_ratio.2>=2 | cn$odds_ratio.2<=0.5)
#cn_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,cn)
cn_f <- cbind(qval_3,or_3,cn)

qval_3 <- (cnu$FDR_qval<=0.05) + (cnu$FDR_qval.1<=0.05)+(cnu$FDR_qval.2<=0.05)
#wt_sum_3 <- (cnu$wt_sum>=100) + (cnu$wt_sum.1>=100)+(cnu$wt_sum.2>=100)
#mt_sum_3 <- (cnu$mt_sum>=100) + (cnu$mt_sum.1>=100)+(cnu$mt_sum.2>=100)
or_3 <- (cnu$odds_ratio>=2 | cnu$odds_ratio<=0.5) + (cnu$odds_ratio.1>=2 | cnu$odds_ratio.1<=0.5)+(cnu$odds_ratio.2>=2 | cnu$odds_ratio.2<=0.5)
#cnu_f <- cbind(qval_3,wt_sum_3,mt_sum_3,or_3,cnu)
cnu_f <- cbind(qval_3,or_3,cnu)

# Generate 3 sets of 200 words
set <- su_f %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% select(mt) %>% unlist()
set1_ <- su_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set1 <- intersect(set1_,set)
set2_ <- cn_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set2 <- intersect(set2_,set)
set3_ <- cnu_f %>% filter(qval_3 == 3) %>% filter(or_3 == 3) %>% select(mt) %>% unlist()
set3 <- intersect(set3_,set)

# input the wt- mt pair 
gtable5 <-read.table("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/5ss-1_new_reference(fisher)", head=F,row.names = 1)
View(gtable5)
colnames(gtable5) <- c("wt","mt")

# mark significant candidates: 
gtable5$spliced_unspliced_sig <- NA
gtable5$spliced_unspliced_sig[which(gtable5$mt %in% set1p)] <- "su"
gtable5$canonical_noncan_sig <- NA
gtable5$canonical_noncan_sig[which(gtable5$mt %in% set2p)] <- "cn"
gtable5$canonical_noncan_un_sig <- NA
gtable5$canonical_noncan_un_sig[which(gtable5$mt %in% set3p)] <- "cnu"

# mark significant candidates pvalue with reads: 
gtable5$spliced_unspliced_fval_strict_reads <- NA
gtable5$spliced_unspliced_fval_strict_reads[which(gtable5$mt %in% set1p_r)] <- "su"
gtable5$canonical_noncan_sig_fval_strict_reads <- NA
gtable5$canonical_noncan_sig_fval_strict_reads[which(gtable5$mt %in% set2p_r)] <- "cn"
gtable5$canonical_noncan_un_sig_fval_strict_reads <- NA
gtable5$canonical_noncan_un_sig_fval_strict_reads[which(gtable5$mt %in% set3p_r)] <- "cnu"

# mark significant candidates by hunglun's criteria: 
gtable5$spliced_unspliced_fdr0.05_reads100_odds2fold <- NA
gtable5$spliced_unspliced_fdr0.05_reads100_odds2fold[which(gtable5$mt %in% set1)] <- "su"
gtable5$canonical_noncan_fdr0.05_reads100_odds2fold <- NA
gtable5$canonical_noncan_fdr0.05_reads100_odds2fold[which(gtable5$mt %in% set2)] <- "cn"
gtable5$canonical_noncan_un_fdr0.05_reads100_odds2fold <- NA
gtable5$canonical_noncan_un_fdr0.05_reads100_odds2fold[which(gtable5$mt %in% set3)] <- "cnu"


# check if most are same spliced form
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output")
s1 <- read.csv(file = "5_1_most_can_non_gtag_only.csv",row.names = 1)
s2 <- read.csv(file = "5_2_most_can_non_gtag_only.csv",row.names = 1)
s3 <- read.csv(file = "5_3_most_can_non_gtag_only.csv",row.names = 1)

most <- na.omit(intersect(intersect(s1$most,s2$most),s3$most))  # 2551
s1c <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/5_1_correct_can_non_gtag_only.csv",row.names = 1)
correct <- s1c[,c(1,2)]
correct$mc <- NA   
correct$mc[which(correct$canonical %in% most)] <- correct$canonical[which(correct$canonical %in% most)]
length(which(correct$canonical %in% most))  # 1680

gtable5 <- merge(gtable5,correct[,c(1,3)],by.x = "wt",by.y = "chrom")
gtable5 <- merge(gtable5,correct[,c(1,3)],by.x = "mt",by.y = "chrom")
colnames(gtable5)[12] <- "wt_can_most"
colnames(gtable5)[13] <- "mt_can_most"

gtable5$can_is_most_wt_mt <- TRUE
gtable5$can_is_most_wt_mt[which(is.na(gtable5$wt_can_most))] <- NA
gtable5$can_is_most_wt_mt[which(is.na(gtable5$mt_can_most))] <- NA

intron <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_intron_fixed2_final_add_GT_AG.csv",row.names = 1)

exon <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_exon_fixed2_final_add_GT_AG.csv",row.names = 1)

intron_f <- intron[,c(103,3,42,101,102,66,65,85,83,38,39)]
colnames(intron_f) <- c("name","origin","5/3ss","can_intron_start","can_intron_end","nchar","fasta_seq","mt_form_5_GT_1_start",
                        "wt_form_5_GT_start","len.left(5ss)","len.right(3ss)")

exon_f <- exon[,c(4,6,68,66,67,37,3,58,46,28,32)]
colnames(exon_f) <- c("name","origin","5/3ss","can_intron_start","can_intron_end","nchar","fasta_seq","mt_form_5_GT_1_start",
                      "wt_form_5_GT_start","len.left(5ss)","len.right(3ss)")
exon_f$origin <- paste0(exon_f$origin,"exon")
exon_f$`5/3ss` <- as.character(exon_f$`5/3ss`)
final <- bind_rows("exon"=exon_f,"intron"=intron_f, .id="exon/intron")
final$distance_5ss <- final$`len.left(5ss)`
final$distance_3ss <- final$`len.right(3ss)`

final3 <- final[which(str_detect(final$`5/3ss`,pattern = "3")),]
final5 <- final[which(str_detect(final$`5/3ss`,pattern = "5")),]
final3$distance_3ss[(which(str_detect(final3$origin,pattern = "inf")))] <- 77
final5$distance_5ss[(which(str_detect(final5$origin,pattern = "inf")))] <- 77


gtable5 <- merge(gtable5,final5[,c(2,5,6,7,8)],by.x="wt",by.y="name",all.x=TRUE)
colnames(gtable5)[15:18] <- c("wt_can_intron_start","wt_can_intron_end","wt_nchar","wt_fasta")

gtable5 <- merge(gtable5,final5[,c(2,5,6,7,8)],by.x="mt",by.y="name",all.x=TRUE)
colnames(gtable5)[19:22] <- c("mt_can_intron_start","mt_can_intron_end","mt_nchar","mt_fasta")

## add GT new 5ss
gtable5 <- merge(gtable5,final5[,c(2,10,9)],by.x="mt",by.y="name",all.x=TRUE)
gtable5$GT <- NA
gtable5$GT[intersect(which(!is.na(gtable5$wt_form_5_GT_start)),which(is.na(gtable5$mt_form_5_GT_1_start)))] <- "reduce_GT"  #328
gtable5$GT[intersect(which(is.na(gtable5$wt_form_5_GT_start)),which(!is.na(gtable5$mt_form_5_GT_1_start)))] <- "add_GT"  #318
intersect(which(!is.na(gtable5$wt_form_5_GT_start)),which(!is.na(gtable5$mt_form_5_GT_1_start)))   # 82 in total
gtable5$GT[intersect(which(!is.na(gtable5$wt_form_5_GT_start)),which(!is.na(gtable5$mt_form_5_GT_1_start)))] <- "reduce/add_GT"


#write.csv(gtable5,file = "summarize/gtag_only/5ss_gtag_only_candidates_100_read_count.csv")

# find add ag effect
cn_gt <- read.csv(file = "can_noncan/gtag_only/wo_add_gt/2022_run3_run4_5ss_gtag_only_mt1_gt_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu_gt <- read.csv(file = "can_noncan_unspliced/gtag_only/wo_add_gt/2022_run3_run4_5ss_gtag_only_mt1_gt_most_can_noncan_un_sort_fvals.csv",row.names = 1)

set4 <- cn_gt %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set5 <- cnu_gt %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()

gt_caused_sig <- unique(c(setdiff(set2p,set4),setdiff(set3p,set5)))  # 24

gtable5$add_gt_cause_significant <- NA
gtable5$add_gt_cause_significant[which(gtable5$mt %in% gt_caused_sig)] <- TRUE   

######### find if there is minus gt cause significant
gtable5$mt[intersect(intersect(unique(c(which(!is.na(gtable5$spliced_unspliced_sig)),which(!is.na(gtable5$canonical_noncan_sig)),which(!is.na(gtable5$canonical_noncan_un_sig)))),which(gtable5$can_is_most_wt_mt==TRUE)),which(gtable5$GT=="reduce_GT"))]  # 119

# find minus gt effect
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/")
cn_wo_gt <- read.csv(file = "can_noncan/gtag_only/w_minus_gt/2022_run3_run4_5ss_wt_minus_gt_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu_wo_gt <- read.csv(file = "can_noncan_unspliced/gtag_only/w_minus_gt/2022_run3_run4_5ss_wt_minus_gt_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)

set6 <- cn_wo_gt %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set7 <- cnu_wo_gt %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()

m_gt_caused_sig <- unique(c(setdiff(set2p,set6),setdiff(set3p,set7)))  # 168

gtable5$minus_gt_cause_significant <- NA
gtable5$minus_gt_cause_significant[which(gtable5$mt %in% m_gt_caused_sig)] <- TRUE   

# input information from the previous final table
# wt.mt to 5'ss distance, position (start,end), use V3 to to calculate 3ss position(extract seq for prediction tool)

gtable5$"wt(ex-in)_delta_gc(percent)" <- 100*((str_count(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_end+1,end =  gtable5$wt_can_intron_end+21), "G") + str_count(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_end+1,end =  gtable5$wt_can_intron_end+21), "C"))/str_length(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_end+1,end =  gtable5$wt_can_intron_end+21)) -
                                               (str_count(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_start+60,end =  gtable5$wt_can_intron_end), "G") + str_count(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_start+60,end =  gtable5$wt_can_intron_end), "C"))/str_length(str_sub(gtable5$wt_fasta,start = gtable5$wt_can_intron_start+60,end =  gtable5$wt_can_intron_end)))
gtable5$"mt(ex-in)_delta_gc(percent)" <- 100*((str_count(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_end+1,end =  gtable5$mt_can_intron_end+21), "G") + str_count(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_end+1,end =  gtable5$mt_can_intron_end+21), "C"))/str_length(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_end+1,end =  gtable5$mt_can_intron_end+21)) -
                                               (str_count(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_start+60,end =  gtable5$mt_can_intron_end), "G") + str_count(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_start+60,end =  gtable5$mt_can_intron_end), "C"))/str_length(str_sub(gtable5$mt_fasta,start = gtable5$mt_can_intron_start+60,end =  gtable5$mt_can_intron_end)))

wt_s_e <- na.omit(cnu[,c(2,8,10,18,20,28,30)])
wt_s_e$"wt_splicing_efficiency" <-100*((wt_s_e[,2]/wt_s_e[,3])+(wt_s_e[,4]/wt_s_e[,5])+(wt_s_e[,6]/wt_s_e[,7]))/3
mt_s_e <- na.omit(cnu[,c(3,11,13,21,23,31,33)])
mt_s_e$"mt_splicing_efficiency" <- 100*((mt_s_e[,2]/mt_s_e[,3])+(mt_s_e[,4]/mt_s_e[,5])+(mt_s_e[,6]/mt_s_e[,7]))/3

gtable5 <- merge(gtable5,distinct(wt_s_e[,c(1,8)]),by.x="wt",by.y="wt",all.x=TRUE)
gtable5 <- merge(gtable5,mt_s_e[,c(1,8)],by.x="mt",by.y="mt",all.x=TRUE)


#write.csv(gtable5,file = "summarize/gtag_only/2022_5ss_candidates_info.csv")
#save(gtable5, file = "summarize/gtag_only/2022_5ss_candidates_1020.RData")
gtable5_1 <- gtable5 %>%
  mutate(delta_efficiency = mt_splicing_efficiency - wt_splicing_efficiency,
         sig_fdr_reads_odds = ifelse(!is.na(spliced_unspliced_fdr0.05_reads100_odds2fold) |
                                       !is.na(canonical_noncan_fdr0.05_reads100_odds2fold) |
                                       !is.na(canonical_noncan_un_fdr0.05_reads100_odds2fold),"sig","non"))
write.csv(gtable5_1,file = "summarize/gtag_only/2025_5ss_candidates_info.csv")
save(gtable5_1, file = "summarize/gtag_only/2025_5ss_candidates_1020.RData")

##### pval and position of mt for splicing vulnerable region
intron_f2 <- intron[,c(103,3,42,4,5,6,7,8)]
colnames(intron_f2) <- c("name","origin","5/3ss","chr","ref_start","ref_end","ref","alt")
exon_f2 <- exon[,c(4,6,68,7,8,9,11,12)]
colnames(exon_f2) <- c("name","origin","5/3ss","chr","ref_start","ref_end","ref","alt")
exon_f2$origin <- paste0(exon_f2$origin,"exon")
exon_f2$`5/3ss` <- as.character(exon_f2$`5/3ss`)
final_2 <- bind_rows("exon"=exon_f2,"intron"=intron_f2, .id="exon/intron")

final3_2 <- final_2[which(str_detect(final_2$`5/3ss`,pattern = "3")),]
final5_2 <- final_2[which(str_detect(final_2$`5/3ss`,pattern = "5")),]

pval <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/can_noncan_unspliced/gtag_only/2022_run3_run4_5ss_gtag_only_most_can_noncan_un_sort_fvals.csv")
library(dplyr)
pvalmax <- pval %>% group_by(mt) %>% dplyr::summarise(max_pval=max(fval,fval.1,fval.2))
final_pval <- merge(final5_2,pvalmax,by.x = "name",by.y = "mt")
final_pval <- merge(final_pval,gtable5[,c(1,5)],by.x = "name",by.y = "mt")
#write.csv(final_pval,file = "summarize/gtag_only/2022_5ss_for_hotspot.csv")
write.csv(final_pval,file = "summarize/gtag_only/2025_5ss_for_hotspot.csv")
final_pval <- merge(final_pval,pval[,c(4,6,7,16,17,26,27)],by.x="name",by.y="mt")
final_pval$id <- NA
final_pval$id[(final_pval$max_pval==final_pval$fval)] <- "s1"
final_pval$id[(final_pval$max_pval==final_pval$fval.1)] <- "s2"
final_pval$id[(final_pval$max_pval==final_pval$fval.2)] <- "s3"
final_pval$max_pval_odds_ratio <- NA
final_pval$max_pval_odds_ratio[which(final_pval$id=="s1")] <- final_pval$odds_ratio[which(final_pval$id=="s1")]
final_pval$max_pval_odds_ratio[which(final_pval$id=="s2")] <- final_pval$odds_ratio.1[which(final_pval$id=="s2")]
final_pval$max_pval_odds_ratio[which(final_pval$id=="s3")] <- final_pval$odds_ratio.2[which(final_pval$id=="s3")]
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/summarize/gtag_only")
#write.csv(final_pval,file = "2022_5ss_for_hotspot_0115.csv")
write.csv(final_pval,file = "2025_5ss_for_hotspot_0115.csv")

####### cosmic intergration #####

setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/summary/")
tsc_5ss <- read.csv("2022_5ss_candidates_info.csv")
tsc_3ss <- read.csv("2022_3ss_candidates_info.csv")

old_cosmic_5ss <- read.csv("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/TSC_cosmic/hi_5ss_candidates_0702_freq_cosmic.csv")
old_cosmic_3ss <- read.csv("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/TSC_cosmic/hi_3ss_candidates_0702_freq_cosmic.csv")

old_cosmic_5ss <- old_cosmic_5ss[,c(2,31,37,38)]
old_cosmic_3ss <- old_cosmic_3ss[,c(2,31,33,34)]

merge_5ss_cosmic <- merge(tsc_5ss,old_cosmic_5ss,by="mt")
merge_3ss_cosmic <- merge(tsc_3ss,old_cosmic_3ss,by="mt")
setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/summary/")
# write.csv(merge_5ss_cosmic[,c(1:31,33:37)],"merge_5ss_cosmic.csv")
# write.csv(merge_3ss_cosmic[,c(1:31,33:37)],"merge_3ss_cosmic.csv")
# merge_5ss_cosmic <- read.csv("merge_5ss_cosmic.csv")
# merge_3ss_cosmic <- read.csv("merge_3ss_cosmic.csv")
write.csv(merge_5ss_cosmic[,c(1:31,33:37)],"2025_merge_5ss_cosmic.csv")
write.csv(merge_3ss_cosmic[,c(1:31,33:37)],"2025_merge_3ss_cosmic.csv")
merge_5ss_cosmic <- read.csv("2025_merge_5ss_cosmic.csv")
merge_3ss_cosmic <- read.csv("2025_merge_3ss_cosmic.csv")


setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/freq_taiwan_biobank/")

tsc2 <- read_xlsx("tsc2.xlsx",sheet="biobank")
tsc1 <- read_xlsx("tsc1.xlsx",sheet="biobank")


biobank_freq <- rbind(tsc1[1:7],tsc2)
colnames(biobank_freq) <- c("biobank_chr","biobank_pos","biobank_ID",
                            "biobank_ref","biobank_individuals",
                            "biobank_freq","biobank_gene")

merge_5ss_cosmic_biobank <-  merge(merge_5ss_cosmic,biobank_freq,by.x="mut_start.1",by.y="biobank_pos",all.x=T)
merge_3ss_cosmic_biobank <-  merge(merge_3ss_cosmic,biobank_freq,by.x="mut_start.1",by.y="biobank_pos",all.x=T)

# write.csv(merge_5ss_cosmic_biobank[,c(2,3,4:32,35:43,33,34)],"merge_5ss_cosmic_biobank.csv")
# write.csv(merge_3ss_cosmic_biobank[,c(2,3,4:32,35:43,33,34)],"merge_3ss_cosmic_biobank.csv")
write.csv(merge_5ss_cosmic_biobank[,c(2,3,4:32,35:43,33,34)],"2025_merge_5ss_cosmic_biobank.csv")
write.csv(merge_3ss_cosmic_biobank[,c(2,3,4:32,35:43,33,34)],"2025_merge_3ss_cosmic_biobank.csv")

###########################check intronic barcode affected ###########################
merge_5ss_cosmic_biobank <- read.csv("merge_5ss_cosmic_biobank.csv")
merge_3ss_cosmic_biobank <- read.csv("merge_3ss_cosmic_biobank.csv")

merge_all_cosmic_biobank <- rbind(merge_5ss_cosmic_biobank[,c(1:3,41)],merge_3ss_cosmic_biobank[,c(1:3,41)])
merge_sig_cosmic_biobank <- merge_all_cosmic_biobank[which(str_detect(merge_all_cosmic_biobank$`sig`,pattern = "sig")),]
uniq_merge_sig_cosmic_biobank <- merge_sig_cosmic_biobank[!duplicated(merge_sig_cosmic_biobank$mt),]

intron_info <- read.csv("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_intron_fixed_final_add_GT_AG 2.csv")
intron_info <- intron_info[,c(2,97,50,66)]
intron_sig <- merge(uniq_merge_sig_cosmic_biobank,intron_info,by.x="mt",by.y="new_name",all.x=TRUE)


# 
# write.csv(intron_sig,"intron_sig.csv",col.names = NULL)
# 








