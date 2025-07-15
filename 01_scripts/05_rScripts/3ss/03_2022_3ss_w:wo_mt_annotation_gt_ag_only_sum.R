library(dplyr)
library(tidyr)
library(stringr)
library(stringi)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/junction/")
t1_ori <- read.table(file = "3-1_S1_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
t2_ori <- read.table(file = "3-2_S2_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
t3_ori <- read.table(file = "3-3_S3_q60_sorted_pairend_mt_info_S.juncbed_exon.annotation",header = TRUE)
#exon_m <- bind_rows("ex_3_1_m"=ex_3_1_m,"ex_3_2_m"=ex_3_2_m,"ex_3_3_m"=ex_3_3_m, .id="id")

#in_3_1_m <- read.table(file = "3-1_t_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_3_2_m <- read.table(file = "3-2_t_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_3_3_m <- read.table(file = "3-3_t_q60_sorted_pairend_mt_info_S.juncbed_intron.annotation",header = TRUE)
#intron_m <- bind_rows("in_3_1_m"=in_3_1_m,"in_3_2_m"=in_3_2_m,"in_3_3_m"=in_3_3_m, .id="id")

# > which(ex_3_1_m$chrom != in_3_1_m$chrom)
#integer(0)
# > which(ex_3_1_m$start != in_3_1_m$start)
#integer(0)
# > which(ex_3_1_m$end != in_3_1_m$end)


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
  else if (str_detect(t2_ori$splice_site[i],
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
ss3_mt_info_barcodes <- hi_mt_info %>%
  filter(library=="3ss") %>%
  transmute(chrom = chrom,
            barcode_start = start,
            barcode_end = end,
            library = library,
            totla_length = totla_length,
            wt = wt,
            avoid_junction_from = start -2,
            avoid_junction_to = totla_length)
# function for barcode regions filtering (in 3'ss library)
filter_junction <- function(data_ori, ss3_info) {
  data_ori %>%
    inner_join(ss3_info, by = "chrom") %>%
    mutate(
      start = as.numeric(start),
      avoid_junction_from = as.numeric(avoid_junction_from),
      avoid_junction_to = as.numeric(avoid_junction_to),
      in_range = end >= avoid_junction_from & end <= avoid_junction_to
    ) %>%
    filter(in_range) %>%
    pull(name) %>%
    unique() -> names_to_remove
  
  data_ori %>% filter(!(name %in% names_to_remove))
}
t1 <- filter_junction(t1_ori, ss3_mt_info_barcodes)
t2 <- filter_junction(t2_ori, ss3_mt_info_barcodes)
t3 <- filter_junction(t3_ori, ss3_mt_info_barcodes)

#
t1$num_gtag <- t1$score
t1$num_gtag[which(t1$splice_site_add!="GT-AG/CT-AC")] <- 0 
t1_f <- t1[which(t1$num_gtag!=0),]

result <-  
  group_by(t1_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%  
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_3_1_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)


t2$num_gtag <- t2$score
t2$num_gtag[which(t2$splice_site_add!="GT-AG/CT-AC")] <- 0 
t2_f <- t2[which(t2$num_gtag!=0),]

result <-  
  group_by(t2_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_3_2_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)



t3$num_gtag <- t3$score
t3$num_gtag[which(t3$splice_site_add!="GT-AG/CT-AC")] <- 0 
t3_f <- t3[which(t3$num_gtag!=0),]

result <-  
  group_by(t3_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_3_3_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)



###############################################################################################################
####### no mt info ############################################################################################

t1_n_ori <- read.table(file = "3-1_S1_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
t2_n_ori <- read.table(file = "3-2_S2_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
t3_n_ori <- read.table(file = "3-3_S3_q60_sorted_pairend_no_mt_info_S.juncbed_exon.annotation",header = TRUE)
#exon_n <- bind_rows("ex_3_1_n"=ex_3_1_n,"ex_3_2_n"=ex_3_2_n,"ex_3_3_n"=ex_3_3_n, .id="id")

#in_3_1_n <- read.table(file = "3-1_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_3_2_n <- read.table(file = "3-2_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#in_3_3_n <- read.table(file = "3-3_t_q60_sorted_pairend_no_mt_info_S.juncbed_intron.annotation",header = TRUE)
#intron_n <- bind_rows("in_3_1_n"=in_3_1_n,"in_3_2_n"=in_3_2_n,"in_3_3_n"=in_3_3_n, .id="id")

# > which(ex_3_1_n$chrom != in_3_1_n$chrom)
#integer(0)
# > which(ex_3_1_n$start != in_3_1_n$start)
#integer(0)
# > which(ex_3_1_n$end != in_3_1_n$end)
#integer(0)

#t1_n <- exon_n[which(str_detect(exon_n$id,pattern = "3_1")),]
#t2_n <- exon_n[which(str_detect(exon_n$id,pattern = "3_2")),]
#t3_n <- exon_n[which(str_detect(exon_n$id,pattern = "3_3")),]

t1_n_ori$strand_add <- t1_n_ori$strand
t1_n_ori$splice_site_add <- t1_n_ori$splice_site
t2_n_ori$strand_add <- t2_n_ori$strand
t2_n_ori$splice_site_add <- t2_n_ori$splice_site
t3_n_ori$strand_add <- t3_n_ori$strand
t3_n_ori$splice_site_add <- t3_n_ori$splice_site


for (i in 1:nrow(t1_n_ori)){
  if (str_detect(t1_n_ori$splice_site[i], 
                 fixed("gt-ag", ignore_case=TRUE))){t1_n_ori$strand_add[i] <- "+"; t1_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
  else if (str_detect(t1_n_ori$splice_site[i], 
                      fixed("ct-ac", ignore_case=TRUE))){t1_n_ori$strand_add[i] <- "-"; t1_n_ori$splice_site_add[i] <- "GT-AG/CT-AC"}
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
# remove the barcode affected junction 
hi_mt_info <- read.table(file = "~/Desktop/phD/TSC_summary/06_others/build oligo library/barcode_mt_3'5'ss/hi_mt_info.bed",
                         header = FALSE) # read the mt info (barcode position)
colnames(hi_mt_info) <- c("chrom","start","end","library","totla_length","wt")
ss3_mt_info_barcodes <- hi_mt_info %>%
  filter(library=="3ss") %>%
  transmute(chrom = chrom,
            barcode_start = start,
            barcode_end = end,
            library = library,
            totla_length = totla_length,
            wt = wt,
            avoid_junction_from = start -2,
            avoid_junction_to = totla_length)
# function for barcode regions filtering (in 3'ss library)
filter_junction <- function(data_ori, ss3_info) {
  data_ori %>%
    inner_join(ss3_info, by = "chrom") %>%
    mutate(
      start = as.numeric(start),
      avoid_junction_from = as.numeric(avoid_junction_from),
      avoid_junction_to = as.numeric(avoid_junction_to),
      in_range = end >= avoid_junction_from & end <= avoid_junction_to
    ) %>%
    filter(in_range) %>%
    pull(name) %>%
    unique() -> names_to_remove
  
  data_ori %>% filter(!(name %in% names_to_remove))
}
t1_n <- filter_junction(t1_n_ori, ss3_mt_info_barcodes)
t2_n <- filter_junction(t2_n_ori, ss3_mt_info_barcodes)
t3_n <- filter_junction(t3_n_ori, ss3_mt_info_barcodes)

#
t1_n$num_gtag <- t1_n$score
t1_n$num_gtag[which(t1_n$splice_site_add!="GT-AG/CT-AC")] <- 0 
t1_n_f <- t1_n[which(t1_n$num_gtag!=0),]

result <-  
  group_by(t1_n_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%   
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_3_1_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)


t2_n$num_gtag <- t2_n$score
t2_n$num_gtag[which(t2_n$splice_site_add!="GT-AG/CT-AC")] <- 0 
t2_n_f <- t2_n[which(t2_n$num_gtag!=0),]

result <-  
  group_by(t2_n_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>% 
  mutate(proportion = 100*num/sum) %>%
  ungroup()  

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_3_2_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)


t3_n$num_gtag <- t3_n$score
t3_n$num_gtag[which(t3_n$splice_site_add!="GT-AG/CT-AC")] <- 0 
t3_n_f <- t3_n[which(t3_n$num_gtag!=0),]

result <-  
  group_by(t3_n_f,chrom,start,end,splice_site_add) %>%      
  summarise(num = sum(num_gtag)) %>%                     
  ungroup() %>%                                                           
  group_by(chrom) %>%                                                       
  mutate(sum = sum(num)) %>%  
  mutate(proportion = 100*num/sum) %>%
  ungroup()

result$start <- result$start+1
result$end <- result$end-1
View(result) 
write.table(result, quote= FALSE,file = "norm_spliced_no_mt_3_3_gtag_only.bed",sep="\t",row.names= FALSE,col.names = TRUE)

