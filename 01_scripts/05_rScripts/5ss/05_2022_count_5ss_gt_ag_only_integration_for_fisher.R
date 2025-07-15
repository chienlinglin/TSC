library("dplyr")
library("tidyr")
library("stringr")

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/junction/")
exon <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_exon_fixed_final_add_GT_AG.csv",row.names = 1)
intron <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/hi_3ss_5ss_intron_fixed_final_add_GT_AG.csv",row.names = 1)

ex_f <- exon[,c(4,53,54,55)]
in_f <- intron[,c(96,94,95,42)]
colnames(ex_f) <- c("name","can_intron_start","can_intron_end","3/5ss")
colnames(in_f) <- c("name","can_intron_start","can_intron_end","3/5ss")

can1 <- rbind(ex_f,in_f)
can <- can1[which(str_detect(can1$`3/5ss`,pattern = "5")),c(1:3)]

can$canonical <- paste0(can$name,can$can_intron_start,can$can_intron_end)



ex_mt <- exon[,c(4,45,46,47,48,53,54,55)]
in_mt <- intron[,c(96,80,81,82,83,94,95,42)]
colnames(ex_mt) <- c("name","gt1_start","gt1_end","gt2_start","gt2_end","can_intron_start","can_intron_end","3/5ss")
colnames(in_mt) <- c("name","gt1_start","gt1_end","gt2_start","gt2_end","can_intron_start","can_intron_end","3/5ss")
mt <- rbind(ex_mt,in_mt)
mt_gt <- mt[which(str_detect(mt$`3/5ss`,pattern = "5")),c(1:7)]

mt_gt <- mt_gt[!is.na(mt_gt$gt1_start),]
mt_gt$gt1 <- paste0(mt_gt$name,mt_gt$gt1_start,mt_gt$can_intron_end)
mt_gt$gt2[!is.na(mt_gt$gt2_start)] <- paste0(mt_gt$name[!is.na(mt_gt$gt2_start)],
                                             mt_gt$gt2_start[!is.na(mt_gt$gt2_start)],
                                             mt_gt$can_intron_end[!is.na(mt_gt$gt2_start)])
mt_gt <- mt_gt[,c(1,8,9)]


barcode_gt <- intron[,c(96,84,85,86,87,94,95,42)]
colnames(barcode_gt) <- c("name","gt1_start","gt1_end","gt2_start","gt2_end","can_intron_start","can_intron_end","3/5ss")
barcode_gt <- barcode_gt[which(str_detect(barcode_gt$`3/5ss`,pattern = "5")),c(1:7)]
barcode_gt$gt <- paste0(barcode_gt$name,barcode_gt$can_intron_start,barcode_gt$gt1_end)
barcode_gt$gt2[!is.na(barcode_gt$gt2_start)] <- paste0(barcode_gt$name[!is.na(barcode_gt$gt2_start)],
                                                       barcode_gt$gt2_start[!is.na(barcode_gt$gt2_start)],
                                                       barcode_gt$can_intron_end[!is.na(barcode_gt$gt2_start)])
barcode_gt <- barcode_gt[,c(1,8,9)]


setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/junction/")
# combine 3 files: mt_info_S, mt_info_N.filtered, no_mt_info_S
# mt_info_spliced+unspliced_total for reference
S <- read.table(file = "2022_5ss_q60_pairend_mt_info.count")
colnames(S) <- c("chr","start","end","s1_t","s1_t_no_mt","s2_t","s2_t_no_mt","s3_t","s3_t_no_mt")
# mt_info_S
s1 <- read.table(file = "norm_spliced_5_1_gtag_only.bed",header = TRUE)
gt1 <- paste0(s1$chrom,s1$start,s1$end)
s1gt <- cbind(s1,gt1)
s1_t <-  
  group_by(s1,chrom,start,end) %>%      
  summarise(num2 = sum(num)) %>%
  mutate(temp = paste0(chrom,start,end)) %>%
  ungroup()
s1gt <- merge(s1gt,s1_t[,c(4,5)],by.x="gt1",by.y="temp")

s1_mt1 <- merge(s1gt,mt_gt,by.x = "gt1",by.y = "gt1")     #149
s1_mt2 <- merge(s1gt,mt_gt,by.x = "gt1",by.y = "gt2")     #0
s1_barcode <- merge(s1gt,barcode_gt,by.x = "gt1",by.y = "gt") #0
s1_barcode2 <- merge(s1gt,barcode_gt,by.x = "gt1",by.y = "gt2") #0
s1_mt1_f <- s1_mt1[,c(1,2,9)]
s1_mt2_f <- s1_mt2[,c(1,2,9)]
colnames(s1_mt1_f) <- c("mt1_+gt_spliced","chrom","mt1_+gt_num")
colnames(s1_mt2_f) <- c("mt2_+gt_spliced","chrom","mt2_+gt_num")

s2 <- read.table(file = "norm_spliced_5_2_gtag_only.bed",header = TRUE)
gt2 <- paste0(s2$chrom,s2$start,s2$end)
s2gt <- cbind(s2,gt2)
s2_t <-  
  group_by(s2,chrom,start,end) %>%      
  summarise(num2 = sum(num)) %>%
  mutate(temp = paste0(chrom,start,end)) %>%
  ungroup()
s2gt <- merge(s2gt,s2_t[,c(4,5)],by.x="gt2",by.y="temp")
s2_mt1 <- merge(s2gt,mt_gt,by.x = "gt2",by.y = "gt1")     # 143
s2_mt2 <- merge(s2gt,mt_gt,by.x = "gt2",by.y = "gt2")     # 2
s2_barcode <- merge(s2gt,barcode_gt,by.x = "gt2",by.y = "gt") #0
s2_barcode2 <- merge(s2gt,barcode_gt,by.x = "gt2",by.y = "gt2") #0
s2_mt1_f <- s2_mt1[,c(1,2,9)]
s2_mt2_f <- s2_mt2[,c(1,2,9)]
colnames(s2_mt1_f) <- c("mt1_+gt_spliced","chrom","mt1_+gt_num")
colnames(s2_mt2_f) <- c("mt2_+gt_spliced","chrom","mt2_+gt_num")

s3 <- read.table(file = "norm_spliced_5_3_gtag_only.bed",header = TRUE)
gt3 <- paste0(s3$chrom,s3$start,s3$end)
s3gt <- cbind(s3,gt3)
s3_t <-  
  group_by(s3,chrom,start,end) %>%      
  summarise(num2 = sum(num)) %>%
  mutate(temp = paste0(chrom,start,end)) %>%
  ungroup()
s3gt <- merge(s3gt,s3_t[,c(4,5)],by.x="gt3",by.y="temp")
s3_mt1 <- merge(s3gt,mt_gt,by.x = "gt3",by.y = "gt1")     # 137
s3_mt2 <- merge(s3gt,mt_gt,by.x = "gt3",by.y = "gt2")     # 1
s3_barcode <- merge(s3gt,barcode_gt,by.x = "gt3",by.y = "gt") #0
s3_barcode2 <- merge(s3gt,barcode_gt,by.x = "gt3",by.y = "gt2") #0
s3_mt1_f <- s3_mt1[,c(1,2,9)]
s3_mt2_f <- s3_mt2[,c(1,2,9)]
colnames(s3_mt1_f) <- c("mt1_+gt_spliced","chrom","mt1_+gt_num")
colnames(s3_mt2_f) <- c("mt2_+gt_spliced","chrom","mt2_+gt_num")



# mt_info_N.filtered
N_filtered <- read.table(file = "2022_5ss_q60_mt_info_N_filtered.count")
colnames(N_filtered) <- c("chr","start","end","s1_unspliced","s2_unspliced","s3_unspliced")


# no_mt_info_S
no_s1 <- read.table(file = "norm_spliced_no_mt_5_1_gtag_only.bed",header = TRUE)
no_s1$gt <- paste0(no_s1$chrom,no_s1$start,no_s1$end)
no_s1_mt1 <- merge(no_s1,mt_gt,by.x = "gt",by.y = "gt1")     # 2
no_s1_mt2 <- merge(no_s1,mt_gt,by.x = "gt",by.y = "gt2")     #0
no_s1_barcode <- merge(no_s1,barcode_gt,by.x = "gt",by.y = "gt") #0
no_s1_barcode2 <- merge(no_s1,barcode_gt,by.x = "gt",by.y = "gt2") #0

no_s2 <- read.table(file = "norm_spliced_no_mt_5_2_gtag_only.bed",header = TRUE)
no_s2$gt <- paste0(no_s2$chrom,no_s2$start,no_s2$end)
no_s2_mt1 <- merge(no_s2,mt_gt,by.x = "gt",by.y = "gt1")     # 1
no_s2_mt2 <- merge(no_s2,mt_gt,by.x = "gt",by.y = "gt2")     #0
no_s2_barcode <- merge(no_s2,barcode_gt,by.x = "gt",by.y = "gt") #0
no_s2_barcode2 <- merge(no_s2,barcode_gt,by.x = "gt",by.y = "gt2") #0

no_s3 <- read.table(file = "norm_spliced_no_mt_5_3_gtag_only.bed",header = TRUE)
no_s3$gt <- paste0(no_s3$chrom,no_s3$start,no_s3$end)
no_s3_mt1 <- merge(no_s3,mt_gt,by.x = "gt",by.y = "gt1")     # 2
no_s3_mt2 <- merge(no_s3,mt_gt,by.x = "gt",by.y = "gt2")     #0
no_s3_barcode <- merge(no_s3,barcode_gt,by.x = "gt",by.y = "gt") #0
no_s3_barcode2 <- merge(no_s3,barcode_gt,by.x = "gt",by.y = "gt2") #0


##################################################################################################################
############## correct as canonical ############################################################################## 
# mt_info_S: first file
can1 <- merge(can,distinct(s1gt[,c(2,7)]),by.x = "name",by.y = "chrom",all.x = TRUE)
can2 <- merge(can,distinct(s2gt[,c(2,7)]),by.x = "name",by.y = "chrom",all.x = TRUE)
can3 <- merge(can,distinct(s3gt[,c(2,7)]),by.x = "name",by.y = "chrom",all.x = TRUE)

s1_correct <- merge(can1,s1gt[,-7],by.y = "gt1",by.x = "canonical",all.x = TRUE)
s2_correct <- merge(can2,s2gt[,-7],by.y = "gt2",by.x = "canonical",all.x = TRUE)
s3_correct <- merge(can3,s3gt[,-7],by.y = "gt3",by.x = "canonical",all.x = TRUE)

s1_correct$chrom[is.na(s1_correct$chrom)] <- s1_correct$name[is.na(s1_correct$chrom)]
s2_correct$chrom[is.na(s2_correct$chrom)] <- s2_correct$name[is.na(s2_correct$chrom)]
s3_correct$chrom[is.na(s3_correct$chrom)] <- s3_correct$name[is.na(s3_correct$chrom)]

# mt_info_N.filtered: second file
s1_correct_2 <- merge(s1_correct,N_filtered[,c(1,4)],by.x = "name",by.y = "chr",all=TRUE)
s2_correct_2 <- merge(s2_correct,N_filtered[,c(1,5)],by.x = "name",by.y = "chr",all=TRUE)
s3_correct_2 <- merge(s3_correct,N_filtered[,c(1,6)],by.x = "name",by.y = "chr",all=TRUE)


# no_mt_info_S: third file 
s1_correct_2 <- merge(s1_correct_2,no_s1_mt1[,c(1,6)],by.x = "canonical",by.y = "gt",all=TRUE)
s2_correct_2 <- merge(s2_correct_2,no_s2_mt1[,c(1,6)],by.x = "canonical",by.y = "gt",all=TRUE)
s3_correct_2 <- merge(s3_correct_2,no_s3_mt1[,c(1,6)],by.x = "canonical",by.y = "gt",all=TRUE)
s1_correct_2$num2[is.na(s1_correct_2$num2)] <- 0
s2_correct_2$num2[is.na(s2_correct_2$num2)] <- 0
s3_correct_2$num2[is.na(s3_correct_2$num2)] <- 0

s1_correct_2$num.y[is.na(s1_correct_2$num.y)] <- 0
s2_correct_2$num.y[is.na(s2_correct_2$num.y)] <- 0
s3_correct_2$num.y[is.na(s3_correct_2$num.y)] <- 0
s1_correct_2$num2 <- s1_correct_2$num2 + s1_correct_2$num.y
s2_correct_2$num2 <- s2_correct_2$num2 + s2_correct_2$num.y
s3_correct_2$num2 <- s3_correct_2$num2 + s3_correct_2$num.y

s1_correct_2 <- s1_correct_2[,-c(2,14)]
s2_correct_2 <- s2_correct_2[,-c(2,14)]
s3_correct_2 <- s3_correct_2[,-c(2,14)]

#calculate unspliced unspliced+noncanonical
#s1_correct_2 <- s1_correct_2[-c(which(duplicated(s1_correct_2[,-c(8,9,10)]))),] ## no duplicate, original code -0
which(duplicated(s1_correct_2[,-c(9,10,11)]))
s1_correct_2$num_canonical <- s1_correct_2$num2
s1_correct_2$sum_spliced <- s1_correct_2$sum
s1_correct_2$sum_spliced[which(is.na(s1_correct_2$sum_spliced))] <- 0
s1_correct_2$noncanonical <- s1_correct_2$sum_spliced - s1_correct_2$num_canonical
s1_correct_2$non_un <- s1_correct_2$noncanonical + s1_correct_2$s1_unspliced
s1_correct_2 <- s1_correct_2[which(!is.na(s1_correct_2$chrom)),]


#s2_correct_2 <- s2_correct_2[-c(which(duplicated(s2_correct_2[,-c(8,9,10)]))),] ## no duplicate, original code -0
which(duplicated(s2_correct_2[,-c(9,10,11)]))
s2_correct_2$num_canonical <- s2_correct_2$num2
s2_correct_2$sum_spliced <- s2_correct_2$sum
s2_correct_2$sum_spliced[which(is.na(s2_correct_2$sum_spliced))] <- 0
s2_correct_2$noncanonical <- s2_correct_2$sum_spliced - s2_correct_2$num_canonical
s2_correct_2$non_un <- s2_correct_2$noncanonical + s2_correct_2$s2_unspliced
s2_correct_2 <- s2_correct_2[which(!is.na(s2_correct_2$chrom)),]


#s3_correct_2 <- s3_correct_2[-c(which(duplicated(s3_correct_2[,-c(8,9,10)]))),] ## no duplicate, original code -0
which(duplicated(s3_correct_2[,-c(9,10,11)]))
s3_correct_2$num_canonical <- s3_correct_2$num2
s3_correct_2$sum_spliced <- s3_correct_2$sum
s3_correct_2$sum_spliced[which(is.na(s3_correct_2$sum_spliced))] <- 0
s3_correct_2$noncanonical <- s3_correct_2$sum_spliced - s3_correct_2$num_canonical
s3_correct_2$non_un <- s3_correct_2$noncanonical + s3_correct_2$s3_unspliced
s3_correct_2 <- s3_correct_2[which(!is.na(s3_correct_2$chrom)),]


# combine add gt file
s1_correct_2 <- merge(s1_correct_2,s1_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s2_correct_2 <- merge(s2_correct_2,s2_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s3_correct_2 <- merge(s3_correct_2,s3_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s1_correct_2$"mt1_+gt_form" <- NA
s2_correct_2$"mt1_+gt_form" <- NA
s3_correct_2$"mt1_+gt_form" <- NA


s1_correct_2 <- merge(s1_correct_2,s1_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s2_correct_2 <- merge(s2_correct_2,s2_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s3_correct_2 <- merge(s3_correct_2,s3_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
s1_correct_2$"mt2_+gt_form" <- NA
s2_correct_2$"mt2_+gt_form" <- NA
s3_correct_2$"mt2_+gt_form" <- NA

# add gt s1
s1_correct_2$"mt1_+gt_form"[which(s1_correct_2$`mt1_+gt_spliced` == s1_correct_2$canonical)] <- "can"
s1_correct_2$"mt1_+gt_form"[which(s1_correct_2$`mt1_+gt_spliced` != s1_correct_2$canonical)] <- "noncan"
s1_correct_2$"mt2_+gt_form"[which(s1_correct_2$`mt2_+gt_spliced` == s1_correct_2$canonical)] <- "can"
s1_correct_2$"mt2_+gt_form"[which(s1_correct_2$`mt2_+gt_spliced` != s1_correct_2$canonical)] <- "noncan"

for (i in 1:nrow(s1_correct_2)){
  if (is.na(s1_correct_2$`mt1_+gt_form`[i])) {a <- s1_correct_2$sum_spliced[i]; b <- s1_correct_2$num_canonical[i];
  c <- s1_correct_2$noncanonical[i]; d <- s1_correct_2$non_un[i]}
  else if (s1_correct_2$`mt1_+gt_form`[i] == "can"){b <- s1_correct_2$num_canonical[i]-s1_correct_2$`mt1_+gt_num`[i];c <- s1_correct_2$noncanonical[i];
  a <- b+c;  d <- c+s1_correct_2$s1_unspliced[i]}
  else if (s1_correct_2$`mt1_+gt_form`[i] == "noncan"){b <- s1_correct_2$num_canonical[i]; c <- s1_correct_2$noncanonical[i]-s1_correct_2$`mt1_+gt_num`[i]; 
  a <- b+c; d <- c+s1_correct_2$s1_unspliced[i]}
  s1_correct_2$"mt1_+gt_spliced_num"[i] <- a
  s1_correct_2$"mt1_+gt_can"[i] <- b
  s1_correct_2$"mt1_+gt_noncan"[i] <- c
  s1_correct_2$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(s1_correct_2)){
  if (is.na(s1_correct_2$`mt2_+gt_form`[i])) {a <- s1_correct_2$sum_spliced[i]; b <- s1_correct_2$num_canonical[i];
  c <- s1_correct_2$noncanonical[i]; d <- s1_correct_2$non_un[i]}
  else if (s1_correct_2$`mt2_+gt_form`[i] == "can"){b <- s1_correct_2$num_canonical[i]-s1_correct_2$`mt2_+gt_num`[i];
  c <- s1_correct_2$noncanonical[i]; a <- b+c; d <- c+s1_correct_2$s1_unspliced[i]}
  else if (s1_correct_2$`mt2_+gt_form`[i] == "noncan"){ b <- s1_correct_2$num_canonical[i];
  c <- s1_correct_2$noncanonical[i]-s1_correct_2$`mt2_+gt_num`[i];a <- b+c; d <- c+s1_correct_2$s1_unspliced[i]}
  s1_correct_2$"mt2_+gt_spliced_num"[i] <- a
  s1_correct_2$"mt2_+gt_can"[i] <- b
  s1_correct_2$"mt2_+gt_noncan"[i] <- c
  s1_correct_2$"mt2_+gt_non_un"[i] <- d
}


# add gt s2
s2_correct_2$"mt1_+gt_form"[which(s2_correct_2$`mt1_+gt_spliced` == s2_correct_2$canonical)] <- "can"
s2_correct_2$"mt1_+gt_form"[which(s2_correct_2$`mt1_+gt_spliced` != s2_correct_2$canonical)] <- "noncan"
s2_correct_2$"mt2_+gt_form"[which(s2_correct_2$`mt2_+gt_spliced` == s2_correct_2$canonical)] <- "can"
s2_correct_2$"mt2_+gt_form"[which(s2_correct_2$`mt2_+gt_spliced` != s2_correct_2$canonical)] <- "noncan"

for (i in 1:nrow(s2_correct_2)){
  if (is.na(s2_correct_2$`mt1_+gt_form`[i])) {a <- s2_correct_2$sum_spliced[i]; b <- s2_correct_2$num_canonical[i];
  c <- s2_correct_2$noncanonical[i]; d <- s2_correct_2$non_un[i]}
  else if (s2_correct_2$`mt1_+gt_form`[i] == "can"){ b <- s2_correct_2$num_canonical[i]-s2_correct_2$`mt1_+gt_num`[i];
  c <- s2_correct_2$noncanonical[i];a <- b+c; d <- c+s2_correct_2$s2_unspliced[i]}
  else if (s2_correct_2$`mt1_+gt_form`[i] == "noncan"){ b <- s2_correct_2$num_canonical[i];
  c <- s2_correct_2$noncanonical[i]-s2_correct_2$`mt1_+gt_num`[i];a <- b+c; d <- c+s2_correct_2$s2_unspliced[i]}
  s2_correct_2$"mt1_+gt_spliced_num"[i] <- a
  s2_correct_2$"mt1_+gt_can"[i] <- b
  s2_correct_2$"mt1_+gt_noncan"[i] <- c
  s2_correct_2$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(s2_correct_2)){
  if (is.na(s2_correct_2$`mt2_+gt_form`[i])) {a <- s2_correct_2$sum_spliced[i]; b <- s2_correct_2$num_canonical[i];
  c <- s2_correct_2$noncanonical[i]; d <- s2_correct_2$non_un[i]}
  else if (s2_correct_2$`mt2_+gt_form`[i] == "can"){ b <- s2_correct_2$num_canonical[i]-s2_correct_2$`mt2_+gt_num`[i];
  c <- s2_correct_2$noncanonical[i];a <- b+c; d <- c+s2_correct_2$s2_unspliced[i]}
  else if (s2_correct_2$`mt2_+gt_form`[i] == "noncan"){ b <- s2_correct_2$num_canonical[i];
  c <- s2_correct_2$noncanonical[i]-s2_correct_2$`mt2_+gt_num`[i];a <- b+c; d <- c+s2_correct_2$s2_unspliced[i]}
  s2_correct_2$"mt2_+gt_spliced_num"[i] <- a
  s2_correct_2$"mt2_+gt_can"[i] <- b
  s2_correct_2$"mt2_+gt_noncan"[i] <- c
  s2_correct_2$"mt2_+gt_non_un"[i] <- d
}

# add gt s3
s3_correct_2$"mt1_+gt_form"[which(s3_correct_2$`mt1_+gt_spliced` == s3_correct_2$canonical)] <- "can"
s3_correct_2$"mt1_+gt_form"[which(s3_correct_2$`mt1_+gt_spliced` != s3_correct_2$canonical)] <- "noncan"
s3_correct_2$"mt2_+gt_form"[which(s3_correct_2$`mt2_+gt_spliced` == s3_correct_2$canonical)] <- "can"
s3_correct_2$"mt2_+gt_form"[which(s3_correct_2$`mt2_+gt_spliced` != s3_correct_2$canonical)] <- "noncan"

for (i in 1:nrow(s3_correct_2)){
  if (is.na(s3_correct_2$`mt1_+gt_form`[i])) {a <- s3_correct_2$sum_spliced[i]; b <- s3_correct_2$num_canonical[i];
  c <- s3_correct_2$noncanonical[i]; d <- s3_correct_2$non_un[i]}
  else if (s3_correct_2$`mt1_+gt_form`[i] == "can"){b <- s3_correct_2$num_canonical[i]-s3_correct_2$`mt1_+gt_num`[i];
  c <- s3_correct_2$noncanonical[i]; a <- b+c; d <- c+s3_correct_2$s3_unspliced[i]}
  else if (s3_correct_2$`mt1_+gt_form`[i] == "noncan"){b <- s3_correct_2$num_canonical[i];
  c <- s3_correct_2$noncanonical[i]-s3_correct_2$`mt1_+gt_num`[i];a <- b+c;  d <- c+s3_correct_2$s3_unspliced[i]}
  s3_correct_2$"mt1_+gt_spliced_num"[i] <- a
  s3_correct_2$"mt1_+gt_can"[i] <- b
  s3_correct_2$"mt1_+gt_noncan"[i] <- c
  s3_correct_2$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(s3_correct_2)){
  if (is.na(s3_correct_2$`mt2_+gt_form`[i])) {a <- s3_correct_2$sum_spliced[i]; b <- s3_correct_2$num_canonical[i];
  c <- s3_correct_2$noncanonical[i]; d <- s3_correct_2$non_un[i]}
  else if (s3_correct_2$`mt2_+gt_form`[i] == "can"){b <- s3_correct_2$num_canonical[i]-s3_correct_2$`mt2_+gt_num`[i];
  c <- s3_correct_2$noncanonical[i];a <- b+c;  d <- c+s3_correct_2$s3_unspliced[i]}
  else if (s3_correct_2$`mt2_+gt_form`[i] == "noncan"){b <- s3_correct_2$num_canonical[i];
  c <- s3_correct_2$noncanonical[i]-s3_correct_2$`mt2_+gt_num`[i]; a <- b+c; d <- c+s3_correct_2$s3_unspliced[i]}
  s3_correct_2$"mt2_+gt_spliced_num"[i] <- a
  s3_correct_2$"mt2_+gt_can"[i] <- b
  s3_correct_2$"mt2_+gt_noncan"[i] <- c
  s3_correct_2$"mt2_+gt_non_un"[i] <- d
}




#write csv correct can_non
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output")
write.csv(s1_correct_2,file = "5_1_correct_can_non_gtag_only.csv",row.names = s1_correct_2[,1])
write.csv(s2_correct_2,file = "5_2_correct_can_non_gtag_only.csv",row.names = s2_correct_2[,1])
write.csv(s3_correct_2,file = "5_3_correct_can_non_gtag_only.csv",row.names = s3_correct_2[,1])




##################################################################################################################
############## most as canonical ############################################################################## 
library(tidyr)
#read normalized table
# mt_info_S: previous already read in
#s1 <- read.table(file = "norm_spliced_5_1_gtag_only.bed",header = TRUE)
#s2 <- read.table(file = "norm_spliced_5_2_gtag_only.bed",header = TRUE)
#s3 <- read.table(file = "norm_spliced_5_3_gtag_only.bed",header = TRUE)
reference <- read.table(file = "/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/5ss-1_new_reference(fisher)")

######################################################################################################
######################### s1 #########################################################################
# find most
most1 <- group_by(s1gt, chrom) %>%  
  top_n(1,num2) %>%                     
  ungroup()

## separate wt/mt pairs
most1_mt <- filter(most1, grepl('to',most1$chrom))
most1_wt <- dplyr::setdiff(most1,most1_mt) 
## separate wt/mt exon/intron 4 pairs

ex_mt5 <- read.table(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/new_5_ss_exon_mut")
n <- group_by(most1_wt,chrom) %>% count(chrom)
n1 <- n[n$n==1,]
n2 <- n[n$n>1,]
repeated_wt <- merge(most1_wt,n2,by.x = "chrom",by.y = "chrom")
repeated_wt2 <- repeated_wt[c(1,3,6,8,9,12,14),]
no_repeat_wt <- merge(most1_wt,n1,by.x = "chrom",by.y = "chrom")
most1_wt2 <- rbind(no_repeat_wt,repeated_wt2)
most1_mt_exon <- merge(most1_mt,ex_mt5,by.x="chrom",by.y="x")
most1_mt_intron <- dplyr::setdiff(most1_mt,most1_mt_exon)
most1_wt_exon <- merge(most1_wt2,ex_mt5,by.x="chrom",by.y="x")
most1_wt_intron <- dplyr::setdiff(most1_wt2,most1_wt_exon)



##count intron mutation mutant character of ref>alt
most1_mt_intron$name <- most1_mt_intron$chrom
intron_size <- separate(most1_mt_intron,"name", into=c("l","r"),sep="_") ## some of the tsc1 var should be modified
intron_size$l <- str_replace(intron_size$l,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr16:",replacement ="") ## new
intron_size <- separate(intron_size,"r", into=c("ref","alt"),sep="to") ## new
View(intron_size)

ref_nchar <- nchar(str_extract(intron_size$ref,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(intron_size$alt,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
intron_size$ref_nchar <- as.numeric(ref_nchar)
intron_size$alt_nchar <- as.numeric(alt_nchar)
intron_size$ref_nchar[is.na(intron_size$ref_nchar)] <- 0
intron_size$alt_nchar[is.na(intron_size$alt_nchar)] <- 0
View(intron_size)


##calculate canonical intron size
intron_size$mt_can <- intron_size$alt_nchar-intron_size$ref_nchar
t <- t(data.frame(strsplit(most1_wt_intron$chrom,"_")))
t <- data.frame(t,most1_wt_intron$start,most1_wt_intron$end)
t$X1 <- str_replace(t$X1,pattern = "chr9:",replacement ="") ## new
intron_size1 <- merge(intron_size,t,by.x="l",by.y="X1")
intron_size1$junc_end <- intron_size1$mt_can+intron_size1$most1_wt_intron.end 

intron_size1$canonical <- paste0(intron_size1$chrom,intron_size1$most1_wt_intron.start,intron_size1$mt_can+intron_size1$most1_wt_intron.end)
dat_intron_mt <- intron_size1[,c(3,8,20)]
dat_intron_mt <- distinct(dat_intron_mt)

##count exon mutation mutant character of ref>alt
most1_mt_exon$name <- most1_mt_exon$chrom
exon_size <- separate(most1_mt_exon,"name", into=c("l","r"),sep="_")
exon_size <- separate(exon_size,"r", into=c("ref","alt"),sep="to")
View(exon_size)

ref_nchar <- nchar(str_extract(exon_size$ref,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(exon_size$alt,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
exon_size$ref_nchar <- as.numeric(ref_nchar)
exon_size$alt_nchar <- as.numeric(alt_nchar)
exon_size$ref_nchar[is.na(exon_size$ref_nchar)] <- 0
exon_size$alt_nchar[is.na(exon_size$alt_nchar)] <- 0
View(exon_size)

##calculate canonical exon size
exon_size$mt_can <- exon_size$alt_nchar-exon_size$ref_nchar
t2 <- t(data.frame(strsplit(most1_wt_exon$chrom,"_")))
t2 <- data.frame(t2,most1_wt_exon$start,most1_wt_exon$end)
exon_size1 <- merge(exon_size,t2,by.x="l",by.y="X1")
exon_size1$junc_start <- exon_size1$mt_can+exon_size1$most1_wt_exon.start
exon_size1$junc_end <- exon_size1$mt_can+exon_size1$most1_wt_exon.end 

exon_size1$canonical <- paste0(exon_size1$chrom,exon_size1$junc_start,exon_size1$junc_end)
dat_exon_mt <- exon_size1[,c(2,8,21)]
dat_exon_mt <- distinct(dat_exon_mt)


dat_ns5_1_intron_mt <- merge(dat_intron_mt,s1gt[,-c(2,7)],by.x="canonical",by.y="gt1",all.x = TRUE)
dat_ns5_1_exon_mt <- merge(dat_exon_mt,s1gt[,-c(2,7)],by.x="canonical",by.y="gt1",all.x = TRUE)
dat_ns5_1_intron_mt$num2[is.na(dat_ns5_1_intron_mt$num2)] <- 0
dat_ns5_1_exon_mt$num2[is.na(dat_ns5_1_exon_mt$num2)] <- 0

colnames(dat_ns5_1_exon_mt)[1] <- "most"
colnames(dat_ns5_1_intron_mt)[1] <- "most"
colnames(most1_wt2)[2] <- "most"
most1_wt2 <- most1_wt2[,-10]
can <- rbind(most1_wt2,dat_ns5_1_exon_mt,dat_ns5_1_intron_mt)


#can_intron_mt <- data.frame(dat_s1_intron_mt[,2:8])
#can <- rbind(most1_wt2[,1:7],can_intron_mt)

## # mt_info_S: first file (most as canonical)

# mt_info_N.filtered: second file
can_1_N_filtered <- merge(can,N_filtered[,c(1,4)],by.x = "chrom",by.y = "chr",all=TRUE)


# no_mt_info_S: third file
can_1_N_filtered <- merge(can_1_N_filtered,no_s1_mt1[,c(1,6)],by.x = "most",by.y = "gt",all=TRUE)
can_1_N_filtered$num.y[is.na(can_1_N_filtered$num.y)] <- 0
colnames(can_1_N_filtered) <- str_replace(string = colnames(can_1_N_filtered),pattern = "num.y",replacement = "no_mt_info_num")
can_1_N_filtered$num2 <- can_1_N_filtered$num2 + can_1_N_filtered$no_mt_info_num
can_1_N_filtered <- can_1_N_filtered[,-11]
## all three noncanonical and only 1,4,1, can be ignored
can_1_N_filtered <- can_1_N_filtered[!is.na(can_1_N_filtered$chrom),]

#calculate unspliced unspliced+noncanonical
#can_1_N_filtered <- can_1_N_filtered[-c(which(duplicated(can_1_N_filtered[,-c(5,6,8)]))),]
which(duplicated(can_1_N_filtered[,-c(5,6,8)]))

can_1_N_filtered$num_canonical <- can_1_N_filtered$num2
can_1_N_filtered$num_canonical[is.na(can_1_N_filtered$num_canonical)] <- 0
can_1_N_filtered$sum_spliced <- can_1_N_filtered$sum
can_1_N_filtered$sum_spliced[is.na(can_1_N_filtered$sum_spliced)] <- 0
can_1_N_filtered$noncanonical <- can_1_N_filtered$sum_spliced - can_1_N_filtered$num_canonical
can_1_N_filtered$non_un <- can_1_N_filtered$noncanonical + can_1_N_filtered$s1_unspliced

# combine add gt file
can_1_N_filtered <- merge(can_1_N_filtered,s1_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_1_N_filtered$"mt1_+gt_form" <- NA
can_1_N_filtered <- merge(can_1_N_filtered,s1_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_1_N_filtered$"mt2_+gt_form" <- NA


# add gt s1
can_1_N_filtered$"mt1_+gt_form"[which(can_1_N_filtered$`mt1_+gt_spliced` == can_1_N_filtered$most)] <- "can"
can_1_N_filtered$"mt1_+gt_form"[which(can_1_N_filtered$`mt1_+gt_spliced` != can_1_N_filtered$most)] <- "noncan"
can_1_N_filtered$"mt2_+gt_form"[which(can_1_N_filtered$`mt2_+gt_spliced` == can_1_N_filtered$most)] <- "can"
can_1_N_filtered$"mt2_+gt_form"[which(can_1_N_filtered$`mt2_+gt_spliced` != can_1_N_filtered$most)] <- "noncan"

for (i in 1:nrow(can_1_N_filtered)){
  if (is.na(can_1_N_filtered$`mt1_+gt_form`[i])) {a <- can_1_N_filtered$sum_spliced[i]; b <- can_1_N_filtered$num_canonical[i];
  c <- can_1_N_filtered$noncanonical[i]; d <- can_1_N_filtered$non_un[i]}
  else if (can_1_N_filtered$`mt1_+gt_form`[i] == "can"){b <- can_1_N_filtered$num_canonical[i]-can_1_N_filtered$`mt1_+gt_num`[i];
  c <- can_1_N_filtered$noncanonical[i]; a <- b+c; d <- c+can_1_N_filtered$s1_unspliced[i]}
  else if (can_1_N_filtered$`mt1_+gt_form`[i] == "noncan"){ b <- can_1_N_filtered$num_canonical[i];
  c <- can_1_N_filtered$noncanonical[i]-can_1_N_filtered$`mt1_+gt_num`[i]; a <- b+c;d <- c+can_1_N_filtered$s1_unspliced[i]}
  can_1_N_filtered$"mt1_+gt_spliced_num"[i] <- a
  can_1_N_filtered$"mt1_+gt_can"[i] <- b
  can_1_N_filtered$"mt1_+gt_noncan"[i] <- c
  can_1_N_filtered$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(can_1_N_filtered)){
  if (is.na(can_1_N_filtered$`mt2_+gt_form`[i])) {a <- can_1_N_filtered$sum_spliced[i]; b <- can_1_N_filtered$num_canonical[i];
  c <- can_1_N_filtered$noncanonical[i]; d <- can_1_N_filtered$non_un[i]}
  else if (can_1_N_filtered$`mt2_+gt_form`[i] == "can"){b <- can_1_N_filtered$num_canonical[i]-can_1_N_filtered$`mt2_+gt_num`[i];
  c <- can_1_N_filtered$noncanonical[i];a <- b+c;  d <- c+can_1_N_filtered$s1_unspliced[i]}
  else if (can_1_N_filtered$`mt2_+gt_form`[i] == "noncan"){b <- can_1_N_filtered$num_canonical[i];
  c <- can_1_N_filtered$noncanonical[i]-can_1_N_filtered$`mt2_+gt_num`[i]; a <- b+c; d <- c+can_1_N_filtered$s1_unspliced[i]}
  can_1_N_filtered$"mt2_+gt_spliced_num"[i] <- a
  can_1_N_filtered$"mt2_+gt_can"[i] <- b
  can_1_N_filtered$"mt2_+gt_noncan"[i] <- c
  can_1_N_filtered$"mt2_+gt_non_un"[i] <- d
}
# repeats in chrom
dup_chroms <- can_1_N_filtered %>%
  count(chrom) %>%
  filter(n > 1) %>%
  pull(chrom)
can_1_N_filtered <- can_1_N_filtered %>%
  filter(!(chrom %in% dup_chroms) | (chrom %in% dup_chroms & start == 65))

#write csv correct can_non
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output")
write.csv(can_1_N_filtered,file = "5_1_most_can_non_gtag_only.csv",row.names = can_1_N_filtered[,1])


######################################################################################################
######################### s2 #########################################################################
# find most
most1 <- group_by(s2gt, chrom) %>%  
  top_n(1,num2) %>%                     
  ungroup()

## separate wt/mt pairs
most1_mt <- filter(most1, grepl('to',most1$chrom))
most1_wt <- dplyr::setdiff(most1,most1_mt) 
## separate wt/mt exon/intron 4 pairs

ex_mt5 <- read.table(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/new_5_ss_exon_mut")
n <- group_by(most1_wt,chrom) %>% count(chrom)
n1 <- n[n$n==1,]
n2 <- n[n$n>1,]
repeated_wt <- merge(most1_wt,n2,by.x = "chrom",by.y = "chrom")
repeated_wt2 <- repeated_wt[c(1,3,5,8,9,12,13),]
no_repeat_wt <- merge(most1_wt,n1,by.x = "chrom",by.y = "chrom")
most1_wt2 <- rbind(no_repeat_wt,repeated_wt2)
most1_mt_exon <- merge(most1_mt,ex_mt5,by.x="chrom",by.y="x")
most1_mt_intron <- dplyr::setdiff(most1_mt,most1_mt_exon)
most1_wt_exon <- merge(most1_wt2,ex_mt5,by.x="chrom",by.y="x")
most1_wt_intron <- dplyr::setdiff(most1_wt2,most1_wt_exon)



##count intron mutation mutant character of ref>alt
most1_mt_intron$name <- most1_mt_intron$chrom
intron_size <- separate(most1_mt_intron,"name", into=c("l","r"),sep="_") ## some of the tsc1 var should be modified
intron_size$l <- str_replace(intron_size$l,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr16:",replacement ="") ## new
intron_size <- separate(intron_size,"r", into=c("ref","alt"),sep="to") ## new
View(intron_size)

ref_nchar <- nchar(str_extract(intron_size$ref,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(intron_size$alt,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
intron_size$ref_nchar <- as.numeric(ref_nchar)
intron_size$alt_nchar <- as.numeric(alt_nchar)
intron_size$ref_nchar[is.na(intron_size$ref_nchar)] <- 0
intron_size$alt_nchar[is.na(intron_size$alt_nchar)] <- 0
View(intron_size)


##calculate canonical intron size
intron_size$mt_can <- intron_size$alt_nchar-intron_size$ref_nchar
t <- t(data.frame(strsplit(most1_wt_intron$chrom,"_")))
t <- data.frame(t,most1_wt_intron$start,most1_wt_intron$end)
t$X1 <- str_replace(t$X1,pattern = "chr9:",replacement ="") ## new
intron_size1 <- merge(intron_size,t,by.x="l",by.y="X1")
intron_size1$junc_end <- intron_size1$mt_can+intron_size1$most1_wt_intron.end 

intron_size1$canonical <- paste0(intron_size1$chrom,intron_size1$most1_wt_intron.start,intron_size1$mt_can+intron_size1$most1_wt_intron.end)
dat_intron_mt <- intron_size1[,c(3,8,20)]
dat_intron_mt <- distinct(dat_intron_mt)

##count exon mutation mutant character of ref>alt
most1_mt_exon$name <- most1_mt_exon$chrom
exon_size <- separate(most1_mt_exon,"name", into=c("l","r"),sep="_")
exon_size <- separate(exon_size,"r", into=c("ref","alt"),sep="to")
View(exon_size)

ref_nchar <- nchar(str_extract(exon_size$ref,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(exon_size$alt,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
exon_size$ref_nchar <- as.numeric(ref_nchar)
exon_size$alt_nchar <- as.numeric(alt_nchar)
exon_size$ref_nchar[is.na(exon_size$ref_nchar)] <- 0
exon_size$alt_nchar[is.na(exon_size$alt_nchar)] <- 0
View(exon_size)

##calculate canonical exon size
exon_size$mt_can <- exon_size$alt_nchar-exon_size$ref_nchar
t2 <- t(data.frame(strsplit(most1_wt_exon$chrom,"_")))
t2 <- data.frame(t2,most1_wt_exon$start,most1_wt_exon$end)
exon_size1 <- merge(exon_size,t2,by.x="l",by.y="X1")
exon_size1$junc_start <- exon_size1$mt_can+exon_size1$most1_wt_exon.start
exon_size1$junc_end <- exon_size1$mt_can+exon_size1$most1_wt_exon.end 

exon_size1$canonical <- paste0(exon_size1$chrom,exon_size1$junc_start,exon_size1$junc_end)
dat_exon_mt <- exon_size1[,c(2,8,21)]
dat_exon_mt <- distinct(dat_exon_mt)


dat_ns5_2_intron_mt <- merge(dat_intron_mt,s2gt[,-c(2,7)],by.x="canonical",by.y="gt2",all.x = TRUE)
dat_ns5_2_exon_mt <- merge(dat_exon_mt,s2gt[,-c(2,7)],by.x="canonical",by.y="gt2",all.x = TRUE)
dat_ns5_2_intron_mt$num2[is.na(dat_ns5_2_intron_mt$num2)] <- 0
dat_ns5_2_exon_mt$num2[is.na(dat_ns5_2_exon_mt$num2)] <- 0

colnames(dat_ns5_2_exon_mt)[1] <- "most"
colnames(dat_ns5_2_intron_mt)[1] <- "most"
colnames(most1_wt2)[2] <- "most"
most1_wt2 <- most1_wt2[,-10]
can <- rbind(most1_wt2,dat_ns5_2_exon_mt,dat_ns5_2_intron_mt)


## # mt_info_S: first file (most as canonical)

# mt_info_N.filtered: second file
can_2_N_filtered <- merge(can,N_filtered[,c(1,5)],by.x = "chrom",by.y = "chr",all=TRUE)


# no_mt_info_S: third file
can_2_N_filtered <- merge(can_2_N_filtered,no_s2_mt1[,c(1,6)],by.x = "most",by.y = "gt",all=TRUE)
can_2_N_filtered$num.y[is.na(can_2_N_filtered$num.y)] <- 0
colnames(can_2_N_filtered) <- str_replace(string = colnames(can_2_N_filtered),pattern = "num.y",replacement = "no_mt_info_num")
can_2_N_filtered$num2 <- can_2_N_filtered$num2 + can_2_N_filtered$no_mt_info_num
can_2_N_filtered <- can_2_N_filtered[,-11]
## only one noncanonical and only 10, can be ignored
can_2_N_filtered <- can_2_N_filtered[!is.na(can_2_N_filtered$chrom),]

#calculate unspliced unspliced+noncanonical
#can_2_N_filtered <- can_2_N_filtered[-c(which(duplicated(can_2_N_filtered[,-c(5,6,8)]))),]
which(duplicated(can_2_N_filtered[,-c(5,6,8)]))
can_2_N_filtered$num_canonical <- can_2_N_filtered$num2
can_2_N_filtered$num_canonical[is.na(can_2_N_filtered$num_canonical)] <- 0
can_2_N_filtered$sum_spliced <- can_2_N_filtered$sum
can_2_N_filtered$sum_spliced[is.na(can_2_N_filtered$sum_spliced)] <- 0
can_2_N_filtered$noncanonical <- can_2_N_filtered$sum_spliced - can_2_N_filtered$num_canonical
can_2_N_filtered$non_un <- can_2_N_filtered$noncanonical + can_2_N_filtered$s2_unspliced

# combine add gt file
can_2_N_filtered <- merge(can_2_N_filtered,s2_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_2_N_filtered$"mt1_+gt_form" <- NA
can_2_N_filtered <- merge(can_2_N_filtered,s2_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_2_N_filtered$"mt2_+gt_form" <- NA


# add gt s2
can_2_N_filtered$"mt1_+gt_form"[which(can_2_N_filtered$`mt1_+gt_spliced` == can_2_N_filtered$most)] <- "can"
can_2_N_filtered$"mt1_+gt_form"[which(can_2_N_filtered$`mt1_+gt_spliced` != can_2_N_filtered$most)] <- "noncan"
can_2_N_filtered$"mt2_+gt_form"[which(can_2_N_filtered$`mt2_+gt_spliced` == can_2_N_filtered$most)] <- "can"
can_2_N_filtered$"mt2_+gt_form"[which(can_2_N_filtered$`mt2_+gt_spliced` != can_2_N_filtered$most)] <- "noncan"

for (i in 1:nrow(can_2_N_filtered)){
  if (is.na(can_2_N_filtered$`mt1_+gt_form`[i])) {a <- can_2_N_filtered$sum_spliced[i]; b <- can_2_N_filtered$num_canonical[i];
  c <- can_2_N_filtered$noncanonical[i]; d <- can_2_N_filtered$non_un[i]}
  else if (can_2_N_filtered$`mt1_+gt_form`[i] == "can"){b <- can_2_N_filtered$num_canonical[i]-can_2_N_filtered$`mt1_+gt_num`[i];
  c <- can_2_N_filtered$noncanonical[i]; a <- b+c; d <- c+can_2_N_filtered$s2_unspliced[i]}
  else if (can_2_N_filtered$`mt1_+gt_form`[i] == "noncan"){ b <- can_2_N_filtered$num_canonical[i];
  c <- can_2_N_filtered$noncanonical[i]-can_2_N_filtered$`mt1_+gt_num`[i]; a <- b+c;d <- c+can_2_N_filtered$s2_unspliced[i]}
  can_2_N_filtered$"mt1_+gt_spliced_num"[i] <- a
  can_2_N_filtered$"mt1_+gt_can"[i] <- b
  can_2_N_filtered$"mt1_+gt_noncan"[i] <- c
  can_2_N_filtered$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(can_2_N_filtered)){
  if (is.na(can_2_N_filtered$`mt2_+gt_form`[i])) {a <- can_2_N_filtered$sum_spliced[i]; b <- can_2_N_filtered$num_canonical[i];
  c <- can_2_N_filtered$noncanonical[i]; d <- can_2_N_filtered$non_un[i]}
  else if (can_2_N_filtered$`mt2_+gt_form`[i] == "can"){b <- can_2_N_filtered$num_canonical[i]-can_2_N_filtered$`mt2_+gt_num`[i];
  c <- can_2_N_filtered$noncanonical[i];a <- b+c;  d <- c+can_2_N_filtered$s2_unspliced[i]}
  else if (can_2_N_filtered$`mt2_+gt_form`[i] == "noncan"){b <- can_2_N_filtered$num_canonical[i];
  c <- can_2_N_filtered$noncanonical[i]-can_2_N_filtered$`mt2_+gt_num`[i]; a <- b+c; d <- c+can_2_N_filtered$s2_unspliced[i]}
  can_2_N_filtered$"mt2_+gt_spliced_num"[i] <- a
  can_2_N_filtered$"mt2_+gt_can"[i] <- b
  can_2_N_filtered$"mt2_+gt_noncan"[i] <- c
  can_2_N_filtered$"mt2_+gt_non_un"[i] <- d
}


#write csv correct can_non
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output")
write.csv(can_2_N_filtered,file = "5_2_most_can_non_gtag_only.csv",row.names = can_2_N_filtered[,1])


######################################################################################################
######################### s3 #########################################################################
# find most
most1 <- group_by(s3gt, chrom) %>%  
  top_n(1,num2) %>%                     
  ungroup()

## separate wt/mt pairs
most1_mt <- filter(most1, grepl('to',most1$chrom))
most1_wt <- dplyr::setdiff(most1,most1_mt) 
## separate wt/mt exon/intron 4 pairs

ex_mt5 <- read.table(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/new_5_ss_exon_mut")
n <- group_by(most1_wt,chrom) %>% count(chrom)
n1 <- n[n$n==1,]
n2 <- n[n$n>1,]
repeated_wt <- merge(most1_wt,n2,by.x = "chrom",by.y = "chrom")
repeated_wt2 <- repeated_wt[c(1,4,5,8,10,12,14),]
no_repeat_wt <- merge(most1_wt,n1,by.x = "chrom",by.y = "chrom")
most1_wt2 <- rbind(no_repeat_wt,repeated_wt2)
most1_mt_exon <- merge(most1_mt,ex_mt5,by.x="chrom",by.y="x")
most1_mt_intron <- dplyr::setdiff(most1_mt,most1_mt_exon)
most1_wt_exon <- merge(most1_wt2,ex_mt5,by.x="chrom",by.y="x")
most1_wt_intron <- dplyr::setdiff(most1_wt2,most1_wt_exon)



##count intron mutation mutant character of ref>alt
most1_mt_intron$name <- most1_mt_intron$chrom
intron_size <- separate(most1_mt_intron,"name", into=c("l","r"),sep="_") ## some of the tsc1 var should be modified
intron_size$l <- str_replace(intron_size$l,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr9:",replacement ="") ## new
intron_size$r <- str_replace(intron_size$r,pattern = "chr16:",replacement ="") ## new
intron_size <- separate(intron_size,"r", into=c("ref","alt"),sep="to") ## new
View(intron_size)

ref_nchar <- nchar(str_extract(intron_size$ref,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(intron_size$alt,"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
intron_size$ref_nchar <- as.numeric(ref_nchar)
intron_size$alt_nchar <- as.numeric(alt_nchar)
intron_size$ref_nchar[is.na(intron_size$ref_nchar)] <- 0
intron_size$alt_nchar[is.na(intron_size$alt_nchar)] <- 0
View(intron_size)


##calculate canonical intron size
intron_size$mt_can <- intron_size$alt_nchar-intron_size$ref_nchar
t <- t(data.frame(strsplit(most1_wt_intron$chrom,"_")))
t <- data.frame(t,most1_wt_intron$start,most1_wt_intron$end)
t$X1 <- str_replace(t$X1,pattern = "chr9:",replacement ="") ## new
intron_size1 <- merge(intron_size,t,by.x="l",by.y="X1")
intron_size1$junc_end <- intron_size1$mt_can+intron_size1$most1_wt_intron.end 

intron_size1$canonical <- paste0(intron_size1$chrom,intron_size1$most1_wt_intron.start,intron_size1$mt_can+intron_size1$most1_wt_intron.end)
dat_intron_mt <- intron_size1[,c(3,8,20)]
dat_intron_mt <- distinct(dat_intron_mt)

##count exon mutation mutant character of ref>alt
most1_mt_exon$name <- most1_mt_exon$chrom
exon_size <- separate(most1_mt_exon,"name", into=c("l","r"),sep="_")
exon_size <- separate(exon_size,"r", into=c("ref","alt"),sep="to")
View(exon_size)

ref_nchar <- nchar(str_extract(exon_size$ref,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(exon_size$alt,"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
exon_size$ref_nchar <- as.numeric(ref_nchar)
exon_size$alt_nchar <- as.numeric(alt_nchar)
exon_size$ref_nchar[is.na(exon_size$ref_nchar)] <- 0
exon_size$alt_nchar[is.na(exon_size$alt_nchar)] <- 0
View(exon_size)

##calculate canonical exon size
exon_size$mt_can <- exon_size$alt_nchar-exon_size$ref_nchar
t2 <- t(data.frame(strsplit(most1_wt_exon$chrom,"_")))
t2 <- data.frame(t2,most1_wt_exon$start,most1_wt_exon$end)
exon_size1 <- merge(exon_size,t2,by.x="l",by.y="X1")
exon_size1$junc_start <- exon_size1$mt_can+exon_size1$most1_wt_exon.start
exon_size1$junc_end <- exon_size1$mt_can+exon_size1$most1_wt_exon.end 


exon_size1$canonical <- paste0(exon_size1$chrom,exon_size1$junc_start,exon_size1$junc_end)
dat_exon_mt <- exon_size1[,c(2,8,21)]
dat_exon_mt <- distinct(dat_exon_mt)


dat_ns5_3_intron_mt <- merge(dat_intron_mt,s3gt[,-c(2,7)],by.x="canonical",by.y="gt3",all.x = TRUE)
dat_ns5_3_exon_mt <- merge(dat_exon_mt,s3gt[,-c(2,7)],by.x="canonical",by.y="gt3",all.x = TRUE)
dat_ns5_3_intron_mt$num2[is.na(dat_ns5_3_intron_mt$num2)] <- 0
dat_ns5_3_exon_mt$num2[is.na(dat_ns5_3_exon_mt$num2)] <- 0

colnames(dat_ns5_3_exon_mt)[1] <- "most"
colnames(dat_ns5_3_intron_mt)[1] <- "most"
colnames(most1_wt2)[2] <- "most"
most1_wt2 <- most1_wt2[,-10]
can <- rbind(most1_wt2,dat_ns5_3_exon_mt,dat_ns5_3_intron_mt)


## # mt_info_S: first file (most as canonical)

# mt_info_N.filtered: second file
can_3_N_filtered <- merge(can,N_filtered[,c(1,6)],by.x = "chrom",by.y = "chr",all=TRUE)


# no_mt_info_S: third file
can_3_N_filtered <- merge(can_3_N_filtered,no_s3_mt1[,c(1,6)],by.x = "most",by.y = "gt",all=TRUE)
can_3_N_filtered$num.y[is.na(can_3_N_filtered$num.y)] <- 0
colnames(can_3_N_filtered) <- str_replace(string = colnames(can_3_N_filtered),pattern = "num.y",replacement = "no_mt_info_num")
can_3_N_filtered$num2 <- can_3_N_filtered$num2 + can_3_N_filtered$no_mt_info_num
can_3_N_filtered <- can_3_N_filtered[,-11]
## all three noncanonical and only 1 can be ignored
can_3_N_filtered <- can_3_N_filtered[!is.na(can_3_N_filtered$chrom),]

#calculate unspliced unspliced+noncanonical
#can_3_N_filtered <- can_3_N_filtered[-c(which(duplicated(can_3_N_filtered[,-c(5,6,8)]))),]
which(duplicated(can_3_N_filtered[,-c(5,6,8)]))
can_3_N_filtered$num_canonical <- can_3_N_filtered$num2
can_3_N_filtered$num_canonical[is.na(can_3_N_filtered$num_canonical)] <- 0
can_3_N_filtered$sum_spliced <- can_3_N_filtered$sum
can_3_N_filtered$sum_spliced[is.na(can_3_N_filtered$sum_spliced)] <- 0
can_3_N_filtered$noncanonical <- can_3_N_filtered$sum_spliced - can_3_N_filtered$num_canonical
can_3_N_filtered$non_un <- can_3_N_filtered$noncanonical + can_3_N_filtered$s3_unspliced

# combine add gt file
can_3_N_filtered <- merge(can_3_N_filtered,s3_mt1_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_3_N_filtered$"mt1_+gt_form" <- NA
can_3_N_filtered <- merge(can_3_N_filtered,s3_mt2_f,by.x = "chrom",by.y = "chrom",all.x = TRUE)
can_3_N_filtered$"mt2_+gt_form" <- NA


# add gt s3
can_3_N_filtered$"mt1_+gt_form"[which(can_3_N_filtered$`mt1_+gt_spliced` == can_3_N_filtered$most)] <- "can"
can_3_N_filtered$"mt1_+gt_form"[which(can_3_N_filtered$`mt1_+gt_spliced` != can_3_N_filtered$most)] <- "noncan"
can_3_N_filtered$"mt2_+gt_form"[which(can_3_N_filtered$`mt2_+gt_spliced` == can_3_N_filtered$most)] <- "can"
can_3_N_filtered$"mt2_+gt_form"[which(can_3_N_filtered$`mt2_+gt_spliced` != can_3_N_filtered$most)] <- "noncan"

for (i in 1:nrow(can_3_N_filtered)){
  if (is.na(can_3_N_filtered$`mt1_+gt_form`[i])) {a <- can_3_N_filtered$sum_spliced[i]; b <- can_3_N_filtered$num_canonical[i];
  c <- can_3_N_filtered$noncanonical[i]; d <- can_3_N_filtered$non_un[i]}
  else if (can_3_N_filtered$`mt1_+gt_form`[i] == "can"){b <- can_3_N_filtered$num_canonical[i]-can_3_N_filtered$`mt1_+gt_num`[i];
  c <- can_3_N_filtered$noncanonical[i]; a <- b+c; d <- c+can_3_N_filtered$s3_unspliced[i]}
  else if (can_3_N_filtered$`mt1_+gt_form`[i] == "noncan"){ b <- can_3_N_filtered$num_canonical[i];
  c <- can_3_N_filtered$noncanonical[i]-can_3_N_filtered$`mt1_+gt_num`[i]; a <- b+c;d <- c+can_3_N_filtered$s3_unspliced[i]}
  can_3_N_filtered$"mt1_+gt_spliced_num"[i] <- a
  can_3_N_filtered$"mt1_+gt_can"[i] <- b
  can_3_N_filtered$"mt1_+gt_noncan"[i] <- c
  can_3_N_filtered$"mt1_+gt_non_un"[i] <- d
}

for (i in 1:nrow(can_3_N_filtered)){
  if (is.na(can_3_N_filtered$`mt2_+gt_form`[i])) {a <- can_3_N_filtered$sum_spliced[i]; b <- can_3_N_filtered$num_canonical[i];
  c <- can_3_N_filtered$noncanonical[i]; d <- can_3_N_filtered$non_un[i]}
  else if (can_3_N_filtered$`mt2_+gt_form`[i] == "can"){b <- can_3_N_filtered$num_canonical[i]-can_3_N_filtered$`mt2_+gt_num`[i];
  c <- can_3_N_filtered$noncanonical[i];a <- b+c;  d <- c+can_3_N_filtered$s3_unspliced[i]}
  else if (can_3_N_filtered$`mt2_+gt_form`[i] == "noncan"){b <- can_3_N_filtered$num_canonical[i];
  c <- can_3_N_filtered$noncanonical[i]-can_3_N_filtered$`mt2_+gt_num`[i]; a <- b+c; d <- c+can_3_N_filtered$s3_unspliced[i]}
  can_3_N_filtered$"mt2_+gt_spliced_num"[i] <- a
  can_3_N_filtered$"mt2_+gt_can"[i] <- b
  can_3_N_filtered$"mt2_+gt_noncan"[i] <- c
  can_3_N_filtered$"mt2_+gt_non_un"[i] <- d
}


#write csv correct can_non
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output")
write.csv(can_3_N_filtered,file = "5_3_most_can_non_gtag_only.csv",row.names = can_3_N_filtered[,1])


