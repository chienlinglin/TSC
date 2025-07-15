library(dplyr)
library(stringr)
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/junction/")
exon <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/barcode_mt_3'5'ss/hi_3ss_5ss_exon_fixed_final_add_GT_AG.csv",row.names = 1)
intron <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/barcode_mt_3'5'ss/hi_3ss_5ss_intron_fixed_final_add_GT_AG.csv",row.names = 1)

ex_f <- exon[,c(4,53,54,55)]
in_f <- intron[,c(96,94,95,42)]
colnames(ex_f) <- c("name","can_intron_start","can_intron_end","3/5ss")
colnames(in_f) <- c("name","can_intron_start","can_intron_end","3/5ss")

can1 <- rbind(ex_f,in_f)
can <- can1[which(str_detect(can1$`3/5ss`,pattern = "3")),c(1:3)]

m1 <- read.table(file = "norm_spliced_3_1.bed",header = TRUE)
m2 <- read.table(file = "norm_spliced_3_2.bed",header = TRUE)
m3 <- read.table(file = "norm_spliced_3_3.bed",header = TRUE)

most1 <-  
  group_by(m1, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most1$chrom) # 36
setdiff(most1$chrom,can$name) # 0
most1$chrom[duplicated(most1$chrom)] # 28

can_m1 <- merge(can,most1,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m1)){
  can_m1$ag_start[i] <- min(can_m1$"can_intron_end"[i],can_m1$end[i])-1
  can_m1$ag_end[i] <- min(can_m1$"can_intron_end"[i],can_m1$end[i])
}


can1_m1_rm <- group_by(can_m1, name,ag_start,ag_end) %>% 
  summarise(min(end)) %>%     
  group_by(name) %>%
  arrange(ag_start) %>%
  dplyr::slice(1) %>%
  ungroup()

can1_m1_rm <- na.omit(can1_m1_rm)
write.table(can1_m1_rm[,c(1:3)], quote= FALSE,file = "3_1_spliced_ag.bed",sep="\t",row.names= FALSE,col.names = FALSE)


## s2
most2 <-  
  group_by(m2, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most2$chrom) # 43
setdiff(most2$chrom,can$name) #0
most2$chrom[duplicated(most2$chrom)] # 43

can_m2 <- merge(can,most2,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m2)){
  can_m2$ag_start[i] <- min(can_m2$"can_intron_end"[i],can_m2$end[i])-1
  can_m2$ag_end[i] <- min(can_m2$"can_intron_end"[i],can_m2$end[i])
}


can2_m2_rm <- group_by(can_m2, name,ag_start,ag_end) %>% 
  summarise(min(end)) %>%     
  group_by(name) %>%
  arrange(ag_start) %>%
  dplyr::slice(1) %>%
  ungroup()

can2_m2_rm <- na.omit(can2_m2_rm)
write.table(can2_m2_rm[,c(1:3)], quote= FALSE,file = "3_2_spliced_ag.bed",sep="\t",row.names= FALSE,col.names = FALSE)



## s3
most3 <-  
  group_by(m3, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most3$chrom) # 27
setdiff(most3$chrom,can$name) #0
most3$chrom[duplicated(most3$chrom)] # 41

can_m3 <- merge(can,most3,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m3)){
  can_m3$ag_start[i] <- min(can_m3$"can_intron_end"[i],can_m3$end[i])-1
  can_m3$ag_end[i] <- min(can_m3$"can_intron_end"[i],can_m3$end[i])
}


can3_m3_rm <- group_by(can_m3, name,ag_start,ag_end) %>% 
  summarise(min(end)) %>%     
  group_by(name) %>%
  arrange(ag_start) %>%
  dplyr::slice(1) %>%
  ungroup()

can3_m3_rm <- na.omit(can3_m3_rm)
write.table(can3_m3_rm[,c(1:3)], quote= FALSE,file = "3_3_spliced_ag.bed",sep="\t",row.names= FALSE,col.names = FALSE)


