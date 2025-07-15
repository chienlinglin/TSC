# Load library
library(dplyr)
library(stringr)
library(Biostrings)
library(readxl)
library(tidyr)

########## filter 100 read count
# input datasets 
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
su <- read.csv(file = "spliced_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_spliced_unspliced_sort_fvals.csv",row.names = 1)
cn <- read.csv(file = "can_nocan/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu <- read.csv(file = "can_noncan_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)

## filter 100 read count
# use su_f
wt_sum_3 <- (su$wt_sum>=100) + (su$wt_sum.1>=100)+(su$wt_sum.2>=100)
mt_sum_3 <- (su$mt_sum>=100) + (su$mt_sum.1>=100)+(su$mt_sum.2>=100)
su_f <- cbind(wt_sum_3,mt_sum_3,su)
su_f_reads <- su_f[,c(1,2,4,5,12,15,22,25,32,35)]
set1 <- su_f_reads %>% filter(wt_sum_3 == 3) %>% filter(mt_sum_3 == 3) %>% select(wt_sum_3,mt_sum_3,mt,wt)
colnames(set1)[1] <-"wt_reads_over_100"
colnames(set1)[2] <-"mt_reads_over_100"
setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/summary/")
write.csv(set1,"3ss_candidates_reads_100.csv",col.names = NULL)

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
gtable3 <-read.table("~/Desktop/phD/Rotation/2/TSC alignment/R script/5ss/new_library/3ss_new_reference(fisher)", head=F,row.names = 1)
View(gtable3)
colnames(gtable3) <- c("wt","mt")

# mark significant candidates
gtable3$spliced_unspliced_sig <- NA
gtable3$spliced_unspliced_sig[which(gtable3$mt %in% set1p)] <- "su"
gtable3$canonical_noncan_sig <- NA
gtable3$canonical_noncan_sig[which(gtable3$mt %in% set2p)] <- "cn"
gtable3$canonical_noncan_un_sig <- NA
gtable3$canonical_noncan_un_sig[which(gtable3$mt %in% set3p)] <- "cnu"

# mark significant candidates pvalue with reads: 
gtable3$spliced_unspliced_fval_strict_reads <- NA
gtable3$spliced_unspliced_fval_strict_reads[which(gtable3$mt %in% set1p_r)] <- "su"
gtable3$canonical_noncan_sig_fval_strict_reads <- NA
gtable3$canonical_noncan_sig_fval_strict_reads[which(gtable3$mt %in% set2p_r)] <- "cn"
gtable3$canonical_noncan_un_sig_fval_strict_reads <- NA
gtable3$canonical_noncan_un_sig_fval_strict_reads[which(gtable3$mt %in% set3p_r)] <- "cnu"

# mark significant candidates by hunglun's criteria
gtable3$spliced_unspliced_fdr0.05_reads100_odds2fold <- NA
gtable3$spliced_unspliced_fdr0.05_reads100_odds2fold[which(gtable3$mt %in% set1)] <- "su"
gtable3$canonical_noncan_fdr0.05_reads100_odds2fold <- NA
gtable3$canonical_noncan_fdr0.05_reads100_odds2fold[which(gtable3$mt %in% set2)] <- "cn"
gtable3$canonical_noncan_un_fdr0.05_reads100_odds2fold <- NA
gtable3$canonical_noncan_un_fdr0.05_reads100_odds2fold[which(gtable3$mt %in% set3)] <- "cnu"

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
# check if most are same spliced form
s1 <- read.csv(file = "3_1_most_can_non_gtag_only.csv",row.names = 1)
s2 <- read.csv(file = "3_2_most_can_non_gtag_only.csv",row.names = 1)
s3 <- read.csv(file = "3_3_most_can_non_gtag_only.csv",row.names = 1)

most <- na.omit(intersect(intersect(s1$most,s2$most),s3$most))  # 2184
s1c <- read.csv(file = "3_1_correct_can_non_gtag_only.csv",row.names = 1)
correct <- s1c[,c(2,3)]
correct$mc <- NA   
correct$mc[which(correct$canonical %in% most)] <- correct$canonical[which(correct$canonical %in% most)]
length(which(correct$canonical %in% most))  # 2024

gtable3 <- merge(gtable3,correct[,c(1,3)],by.x = "wt",by.y = "name")
gtable3 <- merge(gtable3,correct[,c(1,3)],by.x = "mt",by.y = "name")
colnames(gtable3)[12] <- "wt_can_most"
colnames(gtable3)[13] <- "mt_can_most"

gtable3$can_is_most_wt_mt <- TRUE
gtable3$can_is_most_wt_mt[which(is.na(gtable3$wt_can_most))] <- NA
gtable3$can_is_most_wt_mt[which(is.na(gtable3$mt_can_most))] <- NA

intron <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_intron_fixed2_final_add_GT_AG.csv",row.names = 1)
exon <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_exon_fixed2_final_add_GT_AG.csv",row.names = 1)
intron_f <- intron[,c(103,3,42,101,102,66,65,96,94,38,39)]
colnames(intron_f) <- c("name","origin","5/3ss","can_intron_start","can_intron_end","nchar","fasta_seq","mt_form_3_AG_1_end",
                        "wt_form_3_AG_1_end","len.left(5ss)","len.right(3ss)")
exon_f <- exon[,c(4,6,68,66,67,37,3,63,51,28,32)]
colnames(exon_f) <- c("name","origin","5/3ss","can_intron_start","can_intron_end","nchar","fasta_seq","mt_form_3_AG_1_end",
                        "wt_form_3_AG_1_end","len.left(5ss)","len.right(3ss)")
exon_f$origin <- paste0(exon_f$origin,"exon")
exon_f$`5/3ss` <- as.character(exon_f$`5/3ss`)
final <- bind_rows("exon"=exon_f,"intron"=intron_f, .id="exon/intron")
final$distance_5ss <- final$`len.left(5ss)`
final$distance_3ss <- final$`len.right(3ss)`

final3 <- final[which(str_detect(final$`5/3ss`,pattern = "3")),]
final5 <- final[which(str_detect(final$`5/3ss`,pattern = "5")),]
final3$distance_3ss[(which(str_detect(final3$origin,pattern = "inf")))] <- 77
final5$distance_5ss[(which(str_detect(final5$origin,pattern = "inf")))] <- 77


gtable3 <- merge(gtable3,final3[,c(2,5,6,7,8)],by.x="wt",by.y="name",all.x=TRUE)
colnames(gtable3)[15:18] <- c("wt_can_intron_start","wt_can_intron_end","wt_nchar","wt_fasta")

gtable3 <- merge(gtable3,final3[,c(2,5,6,7,8)],by.x="mt",by.y="name",all.x=TRUE)
colnames(gtable3)[19:22] <- c("mt_can_intron_start","mt_can_intron_end","mt_nchar","mt_fasta")

## add AG new 3ss
gtable3 <- merge(gtable3,final3[,c(2,10,9)],by.x="mt",by.y="name",all.x=TRUE)
gtable3$AG <- NA
gtable3$AG[intersect(which(!is.na(gtable3$wt_form_3_AG_1_end)),which(is.na(gtable3$mt_form_3_AG_1_end)))] <- "reduce_AG"  #305
gtable3$AG[intersect(which(!is.na(gtable3$mt_form_3_AG_1_end)),which(is.na(gtable3$wt_form_3_AG_1_end)))] <- "add_AG"  #189
intersect(which(!is.na(gtable3$wt_form_3_AG_1_end)),which(!is.na(gtable3$mt_form_3_AG_1_end)))    # 79 in total, 28 =,38(AG REALTED),1:A>AA
#[10]A>G  [14]A>-  [29] G>A  [94]G>A  [118]A>-  [119]A>G  [192]G>A  [212]G>A  [303]G>A  [316]G>A  [351]G>A  [363]G>A [402]G>A
#[428]G>A  [472]G>A  [495]A>-  [510]G>A  [534]G>A  [536]AA>-  [551]A>G  [584]AG>-  [585]->A  [617]A>G  [629]->A  [667]G>A
#[711]G>A  [716]->A  [722]G>A  [766]G>A  [776]G>A  [785]AA>GAG  [825]G>A  [828]AA>-  [890]->A  [900]->A  [953]G>A  [958]->GT
#[984]->A  [1009]G>A  [1073]AG>-  [1074]G>A  [1113]G>A  [1122]G>A  [1173]CA>-  [1177]G>A  [1197]G>A  [1232]->GGCA  [1368]->A
#[1375]G>A  [1384]->CTTA  [1612]->A  [1615]->A  [1625]G>A  [1690]CCA>AG  [1739]CA>-  [1740]->CACA  [1862]A>-  [1890]->A
#[1892]A>-  [1893]A>G  [1894]AA>-  [1911]G>A  [1930]G>A  [1948]AACA>-  [1998]G>A  [2014]->A  [2059]G>A  [2105]A>-  [2120]A>G
#[2127]A>-  [2136]G>A  [2155]A>-  [2156]A>G  [2159]G>A  [2181]G>A  [2209]->A  [2234]G>A  [2239]AG>CAGTCCT  [2292]G>A
gtable3$AG[intersect(which(!is.na(gtable3$wt_form_3_AG_1_end)),which(!is.na(gtable3$mt_form_3_AG_1_end)))] <- "reduce/add_AG"


#write.csv(gtable3,file = "summarize/3ss_hunglun_gtag_only_candidates_100_read_count.csv")

# find add ag effect
cn_ag <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_nocan/gtag_only/wo_add_ag/2022_run3_run4_3ss_gtag_only_mt1_ag_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu_ag <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/2022_run3_run4_3ss_gtag_only_mt1_ag_most_can_noncan_un_sort_fvals.csv",row.names = 1)

set4 <- cn_ag %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set5 <- cnu_ag %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()

ag_caused_sig <- unique(c(setdiff(set2p,set4),setdiff(set3p,set5)))  # 39

gtable3$add_ag_cause_significant <- NA
gtable3$add_ag_cause_significant[which(gtable3$mt %in% ag_caused_sig)] <- TRUE   

######### find if there is minus ag cause significant
gtable3$mt[intersect(intersect(unique(c(which(!is.na(gtable3$spliced_unspliced_sig)),which(!is.na(gtable3$canonical_noncan_sig)),which(!is.na(gtable3$canonical_noncan_un_sig)))),which(gtable3$can_is_most_wt_mt==TRUE)),which(gtable3$AG=="reduce_AG"))]  #110

# find minus ag effect
cn_wo_ag <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_nocan/gtag_only/wo_add_ag/2022_run3_run4_3ss_wt_minus_ag_gtag_only_most_can_noncan_sort_fvals.csv",row.names = 1)
cnu_wo_ag <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/w_minus_ag/2022_run3_run4_3ss_wt_minus_ag_gtag_only_most_can_noncan_un_sort_fvals.csv",row.names = 1)

set6 <- cn_wo_ag %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()
set7 <- cnu_wo_ag %>% filter(pass_3_fval_strict == 3) %>% select(mt) %>% unlist()

m_ag_caused_sig <- unique(c(setdiff(set2p,set6),setdiff(set3p,set7)))  # 134

gtable3$minus_ag_cause_significant <- NA
gtable3$minus_ag_cause_significant[which(gtable3$mt %in% m_ag_caused_sig)] <- TRUE   

# input information from the previous final table
# wt.mt to 3'ss distance, position (start,end), use V3 to to calculate 3ss position(extract seq for prediction tool)

gtable3$"wt(ex-in)_delta_gc(percent)" <- 100*((str_count(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_end+1,end =  gtable3$wt_can_intron_end+21), "G") + str_count(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_end+1,end =  gtable3$wt_can_intron_end+21), "C"))/str_length(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_end+1,end =  gtable3$wt_can_intron_end+21)) -
                                               (str_count(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_start+60,end =  gtable3$wt_can_intron_end), "G") + str_count(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_start+60,end =  gtable3$wt_can_intron_end), "C"))/str_length(str_sub(gtable3$wt_fasta,start = gtable3$wt_can_intron_start+60,end =  gtable3$wt_can_intron_end)))
gtable3$"mt(ex-in)_delta_gc(percent)" <- 100*((str_count(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_end+1,end =  gtable3$mt_can_intron_end+21), "G") + str_count(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_end+1,end =  gtable3$mt_can_intron_end+21), "C"))/str_length(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_end+1,end =  gtable3$mt_can_intron_end+21)) -
                                               (str_count(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_start+60,end =  gtable3$mt_can_intron_end), "G") + str_count(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_start+60,end =  gtable3$mt_can_intron_end), "C"))/str_length(str_sub(gtable3$mt_fasta,start = gtable3$mt_can_intron_start+60,end =  gtable3$mt_can_intron_end)))

wt_s_e <- na.omit(cnu[,c(2,8,10,18,20,28,30)])
wt_s_e$"wt_splicing_efficiency" <-100*((wt_s_e[,2]/wt_s_e[,3])+(wt_s_e[,4]/wt_s_e[,5])+(wt_s_e[,6]/wt_s_e[,7]))/3
mt_s_e <- na.omit(cnu[,c(3,11,13,21,23,31,33)])
mt_s_e$"mt_splicing_efficiency" <- 100*((mt_s_e[,2]/mt_s_e[,3])+(mt_s_e[,4]/mt_s_e[,5])+(mt_s_e[,6]/mt_s_e[,7]))/3

gtable3 <- merge(gtable3,distinct(wt_s_e[,c(1,8)]),by.x="wt",by.y="wt",all.x=TRUE)
gtable3 <- merge(gtable3,mt_s_e[,c(1,8)],by.x="mt",by.y="mt",all.x=TRUE)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
# write.csv(gtable3,file = "summarize/gtag_only/2022_3ss_candidates_info.csv")
# save(gtable3, file = "summarize/gtag_only/2022_3ss_candidates_1020.RData")

gtable3_1 <- gtable3 %>%
  mutate(delta_efficiency = mt_splicing_efficiency - wt_splicing_efficiency,
         sig_fdr_reads_odds = ifelse(!is.na(spliced_unspliced_fdr0.05_reads100_odds2fold) |
                                       !is.na(canonical_noncan_fdr0.05_reads100_odds2fold) |
                                       !is.na(canonical_noncan_un_fdr0.05_reads100_odds2fold),"sig","non"))
write.csv(gtable3_1,file = "summarize/gtag_only/2025_3ss_candidates_info.csv")
save(gtable3_1, file = "summarize/gtag_only/2025_3ss_candidates_1020.RData")

###### pval and position of mt for splicing vulnerable region ######
intron_f2 <- intron[,c(103,3,42,4,5,6,7,8)]
colnames(intron_f2) <- c("name","origin","5/3ss","chr","ref_start","ref_end","ref","alt")
exon_f2 <- exon[,c(4,6,68,7,8,9,11,12)]
colnames(exon_f2) <- c("name","origin","5/3ss","chr","ref_start","ref_end","ref","alt")
exon_f2$origin <- paste0(exon_f2$origin,"exon")
exon_f2$`5/3ss` <- as.character(exon_f2$`5/3ss`)
final_2 <- bind_rows("exon"=exon_f2,"intron"=intron_f2, .id="exon/intron")

final3_2 <- final_2[which(str_detect(final_2$`5/3ss`,pattern = "3")),]
final5_2 <- final_2[which(str_detect(final_2$`5/3ss`,pattern = "5")),]

pval <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/2022_run3_run4_3ss_gtag_only_most_can_noncan_un_sort_fvals.csv")
library(dplyr)
pvalmax <- pval %>% group_by(mt) %>% dplyr::summarise(max_pval=max(fval,fval.1,fval.2))
final_pval <- merge(final3_2,pvalmax,by.x = "name",by.y = "mt")
final_pval <- merge(final_pval,gtable3[,c(1,5)],by.x = "name",by.y = "mt")
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
write.csv(final_pval,file = "summarize/gtag_only/2022_3ss_for_hotspot.csv")

final_pval <- merge(final_pval,pval[,c(4,6,7,16,17,26,27)],by.x="name",by.y="mt")
final_pval$id <- NA
final_pval$id[(final_pval$max_pval==final_pval$fval)] <- "s1"
final_pval$id[(final_pval$max_pval==final_pval$fval.1)] <- "s2"
final_pval$id[(final_pval$max_pval==final_pval$fval.2)] <- "s3"
final_pval$max_pval_odds_ratio <- NA
final_pval$max_pval_odds_ratio[which(final_pval$id=="s1")] <- final_pval$odds_ratio[which(final_pval$id=="s1")]
final_pval$max_pval_odds_ratio[which(final_pval$id=="s2")] <- final_pval$odds_ratio.1[which(final_pval$id=="s2")]
final_pval$max_pval_odds_ratio[which(final_pval$id=="s3")] <- final_pval$odds_ratio.2[which(final_pval$id=="s3")]
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/summarize/gtag_only/")
write.csv(final_pval,file = "2022_3ss_for_hotspot_0115.csv")
