library(dplyr)
library(stringr)
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/junction")
exon <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/barcode_mt_3'5'ss/hi_3ss_5ss_exon_fixed_final_add_GT_AG.csv",row.names = 1)
intron <- read.csv(file = "~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/barcode_mt_3'5'ss/hi_3ss_5ss_intron_fixed_final_add_GT_AG.csv",row.names = 1)

ex_f <- exon[,c(4,53,54,55)]
in_f <- intron[,c(64,94,95,42)]
colnames(ex_f) <- c("name","can_intron_start","can_intron_end","3/5ss")
colnames(in_f) <- c("name","can_intron_start","can_intron_end","3/5ss")
in_f$name <- str_replace(in_f$name,pattern = ">",replacement = "")

can1 <- rbind(ex_f,in_f)
can <- can1[which(str_detect(can1$`3/5ss`,pattern = "5")),c(1:3)]

m1 <- read.table(file = "norm_spliced_5_1.bed",header = TRUE)
m2 <- read.table(file = "norm_spliced_5_2.bed",header = TRUE)
m3 <- read.table(file = "norm_spliced_5_3.bed",header = TRUE)

most1 <-  
  group_by(m1, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most1$chrom) # 32
setdiff(most1$chrom,can$name) #0
most1$chrom[duplicated(most1$chrom)] # 23

can_m1 <- merge(can,most1,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m1)){
  can_m1$gt_start[i] <- max(can_m1$"can_intron_start"[i],can_m1$start[i])
  can_m1$gt_end[i] <- max(can_m1$"can_intron_start"[i],can_m1$start[i])+1
}


can1_m1_rm <- group_by(can_m1, name,gt_start,gt_end) %>% 
  summarise(max(start)) %>%     
  group_by(name) %>%
  arrange(desc(gt_start)) %>%
  slice(1) %>%
  ungroup()

can1_m1_rm <- na.omit(can1_m1_rm)
write.table(can1_m1_rm[,c(1:3)], quote= FALSE,file = "5_1_spliced_gt.bed",sep="\t",row.names= FALSE,col.names = FALSE)


## s2
most2 <-  
  group_by(m2, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most2$chrom) # 33
setdiff(most2$chrom,can$name) # 0
most2$chrom[duplicated(most2$chrom)] # 25

can_m2 <- merge(can,most2,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m2)){
  can_m2$gt_start[i] <- max(can_m2$"can_intron_start"[i],can_m2$start[i])
  can_m2$gt_end[i] <- max(can_m2$"can_intron_start"[i],can_m2$start[i])+1
}


can2_m2_rm <- group_by(can_m2, name,gt_start,gt_end) %>% 
  summarise(max(start)) %>%     
  group_by(name) %>%
  arrange(desc(gt_start)) %>%
  slice(1) %>%
  ungroup()

can2_m2_rm <- na.omit(can2_m2_rm)
write.table(can2_m2_rm[,c(1:3)], quote= FALSE,file = "5_2_spliced_gt.bed",sep="\t",row.names= FALSE,col.names = FALSE)



## s3
most3 <-  
  group_by(m3, chrom) %>%      
  top_n(1,num) %>%                     
  ungroup()

setdiff(can$name,most3$chrom) # 27
setdiff(most3$chrom,can$name) # 0
most3$chrom[duplicated(most3$chrom)] # 15

can_m3 <- merge(can,most3,by.x = "name",by.y = "chrom",all = TRUE)

for (i in 1:nrow(can_m3)){
  can_m3$gt_start[i] <- max(can_m3$"can_intron_start"[i],can_m3$start[i])
  can_m3$gt_end[i] <- max(can_m3$"can_intron_start"[i],can_m3$start[i])+1
}


can3_m3_rm <- group_by(can_m3, name,gt_start,gt_end) %>% 
  summarise(max(start)) %>%     
  group_by(name) %>%
  arrange(desc(gt_start)) %>%
  slice(1) %>%
  ungroup()

can3_m3_rm <- na.omit(can3_m3_rm)
write.table(can3_m3_rm[,c(1:3)], quote= FALSE,file = "5_3_spliced_gt.bed",sep="\t",row.names= FALSE,col.names = FALSE)

# > sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS 15.4.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] zh_TW.UTF-8/zh_TW.UTF-8/zh_TW.UTF-8/C/zh_TW.UTF-8/zh_TW.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] stringr_1.5.0 dplyr_1.1.2  
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.4      utf8_1.2.3       R6_2.5.1         lifecycle_1.0.3  magrittr_2.0.3   pillar_1.9.0     stringi_1.7.12  
# [8] rlang_1.1.1      cli_3.6.1        rstudioapi_0.14  vctrs_0.6.2      generics_0.1.3   tools_4.2.0      glue_1.6.2      
# [15] compiler_4.2.0   pkgconfig_2.0.3  tidyselect_1.2.0 tibble_3.2.1   