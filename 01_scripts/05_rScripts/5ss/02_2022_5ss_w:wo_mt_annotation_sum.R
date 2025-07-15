library(dplyr)
library(tidyr)
library(stringr)
library(stringi)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/junction")
t1_ori <- read.table(file = "5-1_S4_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
t2_ori <- read.table(file = "5-2_S5_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
t3_ori <- read.table(file = "5-3_S6_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
#exon_m <- bind_rows("ex_5_1_m"=ex_5_1_m,"ex_5_2_m"=ex_5_2_m,"ex_5_3_m"=ex_5_3_m, .id="id")

#in_5_1_m <- read.table(file = "5-1_S4_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_5_2_m <- read.table(file = "5-2_S5_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_5_3_m <- read.table(file = "5-3_S6_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#intron_m <- bind_rows("in_5_1_m"=in_5_1_m,"in_5_2_m"=in_5_2_m,"in_5_3_m"=in_5_3_m, .id="id")

# intron exon file same, because giving total fasta file and exon/intron gtf didn't work

#t1 <- exon_m[which(str_detect(exon_m$id,pattern = "5_1")),]
#t2 <- exon_m[which(str_detect(exon_m$id,pattern = "5_2")),]
#t3 <- exon_m[which(str_detect(exon_m$id,pattern = "5_3")),]

t1_ori$strand_add <- t1_ori$strand
t1_ori$splice_site_add <- t1_ori$splice_site
t2_ori$strand_add <- t2_ori$strand
t2_ori$splice_site_add <- t2_ori$splice_site
t3_ori$strand_add <- t3_ori$strand
t3_ori$splice_site_add <- t3_ori$splice_site


for (i in 1:nrow(t1_ori)){
  if (str_detect(t1_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t1_ori$strand_add[i] <- "+"; t1_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t1_ori$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t1_ori$strand_add[i] <- "-"; t1_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 

for (i in 1:nrow(t2_ori)){
  if (str_detect(t2_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t2_ori$strand_add[i] <- "+"; t2_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t2$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t2_ori$strand_add[i] <- "-"; t2_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 

for (i in 1:nrow(t3_ori)){
  if (str_detect(t3_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t3_ori$strand_add[i] <- "+"; t3_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t3_ori$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t3_ori$strand_add[i] <- "-"; t3_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 
# remove the barcode affected junction 
hi_mt_info <- read.table(file = "~/Desktop/phD/TSC_summary/06_others/build oligo library/barcode_mt_3'5'ss/hi_mt_info.bed",
                 header = FALSE) # read the mt info (barcode position)
colnames(hi_mt_info) <- c("chrom","start","end","library","totla_length","wt")
ss5_mt_info_barcodes <- hi_mt_info %>%
  filter(library=="5ss") %>%
  transmute(chrom = chrom,
            barcode_start = start,
            barcode_end = end,
            library = library,
            totla_length = totla_length,
            wt = wt,
            avoid_junction_from = 1,
            avoid_junction_to = end + 2)
# function for barcode regions filtering (in 5'ss library)
filter_junction <- function(data_ori, ss5_info) {
  data_ori %>%
    inner_join(ss5_info, by = "chrom") %>%
    mutate(
      start = as.numeric(start),
      avoid_junction_from = as.numeric(avoid_junction_from),
      avoid_junction_to = as.numeric(avoid_junction_to),
      in_range = start >= avoid_junction_from & start <= avoid_junction_to
    ) %>%
    filter(in_range) %>%
    pull(name) %>%
    unique() -> names_to_remove
  
  data_ori %>% filter(!(name %in% names_to_remove))
}
t1 <- filter_junction(t1_ori, ss5_mt_info_barcodes)
t2 <- filter_junction(t2_ori, ss5_mt_info_barcodes)
t3 <- filter_junction(t3_ori, ss5_mt_info_barcodes)

#
result <-  
  group_by(t1,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%  
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_5_1.bed",sep="\t",row.names= FALSE,col.names = TRUE)


result <-  
  group_by(t2,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_5_2.bed",sep="\t",row.names= FALSE,col.names = TRUE)


result <-  
  group_by(t3,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_5_3.bed",sep="\t",row.names= FALSE,col.names = TRUE)



###############################################################################################################
####### no mt info ############################################################################################

t1_n_ori <- read.table(file = "5-1_S4_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
t2_n_ori <- read.table(file = "5-2_S5_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
t3_n_ori <- read.table(file = "5-3_S6_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
#exon_n <- bind_rows("ex_5_1_n"=ex_5_1_n,"ex_5_2_n"=ex_5_2_n,"ex_5_3_n"=ex_5_3_n, .id="id")

#in_5_1_n <- read.table(file = "5-1_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_5_2_n <- read.table(file = "5-2_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_5_3_n <- read.table(file = "5-3_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#intron_n <- bind_rows("in_5_1_n"=in_5_1_n,"in_5_2_n"=in_5_2_n,"in_5_3_n"=in_5_3_n, .id="id")

# intron exon file same, because giving total fasta file and exon/intron gtf didn't work

#t1_n <- exon_n[which(str_detect(exon_n$id,pattern = "5_1")),]
#t2_n <- exon_n[which(str_detect(exon_n$id,pattern = "5_2")),]
#t3_n <- exon_n[which(str_detect(exon_n$id,pattern = "5_3")),]

t1_n_ori$strand_add <- t1_n_ori$strand
t1_n_ori$splice_site_add <- t1_n_ori$splice_site
t2_n_ori$strand_add <- t2_n_ori$strand
t2_n_ori$splice_site_add <- t2_n_ori$splice_site
t3_n_ori$strand_add <- t3_n_ori$strand
t3_n_ori$splice_site_add <- t3_n_ori$splice_site



for (i in 1:nrow(t1_n_ori)){
  if (str_detect(t1_n_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t1_n_ori$strand_add[i] <- "+"; t1_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t1_n_ori$splice_site[i], fixed("ct-ac", ignore_case=TRUE))){t1_n_ori$strand_add[i] <- "-"; t1_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 

for (i in 1:nrow(t2_n_ori)){
  if (str_detect(t2_n_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t2_n_ori$strand_add[i] <- "+"; t2_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t2_n_ori$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t2_n_ori$strand_add[i] <- "-"; t2_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 

for (i in 1:nrow(t3_n_ori)){
  if (str_detect(t3_n_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t3_n_ori$strand_add[i] <- "+"; t3_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t3_n_ori$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t3_n_ori$strand_add[i] <- "-"; t3_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  print(i)
} 

# function for barcode regions filtering (in 5'ss library)
filter_junction <- function(data_ori, ss5_info) {
  data_ori %>%
    inner_join(ss5_info, by = "chrom") %>%
    mutate(
      start = as.numeric(start),
      avoid_junction_from = as.numeric(avoid_junction_from),
      avoid_junction_to = as.numeric(avoid_junction_to),
      in_range = start >= avoid_junction_from & start <= avoid_junction_to
    ) %>%
    filter(in_range) %>%
    pull(name) %>%
    unique() -> names_to_remove
  
  data_ori %>% filter(!(name %in% names_to_remove))
}
t1_n <- filter_junction(t1_n_ori, ss5_mt_info_barcodes)
t2_n <- filter_junction(t2_n_ori, ss5_mt_info_barcodes)
t3_n <- filter_junction(t3_n_ori, ss5_mt_info_barcodes)

#
result <-  
  group_by(t1_n,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%   
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_5_1.bed",sep="\t",row.names= FALSE,col.names = TRUE)


result <-  
  group_by(t2_n,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_5_2.bed",sep="\t",row.names= FALSE,col.names = TRUE)


result <-  
  group_by(t3_n,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(score)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%  
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_5_3.bed",sep="\t",row.names= FALSE,col.names = TRUE)

######## sessionInfo()
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
#   [1] stringi_1.7.12 stringr_1.5.0  tidyr_1.3.0    dplyr_1.1.2   
# 
# loaded via a namespace (and not attached):
#   [1] fansi_1.0.4      utf8_1.2.3       R6_2.5.1         lifecycle_1.0.3  magrittr_2.0.3   pillar_1.9.0     rlang_1.1.1     
# [8] cli_3.6.1        rstudioapi_0.14  vctrs_0.6.2      generics_0.1.3   tools_4.2.0      glue_1.6.2       purrr_1.0.1     
# [15] compiler_4.2.0   pkgconfig_2.0.3  tidyselect_1.2.0 tibble_3.2.1    