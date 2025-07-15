library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
options(stringsAsFactors = FALSE)

setwd("~/Documents/yiting/little Hi/build oligo library/")
tsc1_5 <- read.csv("TSC1_exon3_95-115.csv")
tsc1_3 <- read.csv("TSC1_exon5_65-115.csv")
tsc2_5 <- read.csv("TSC2_exon3_95-uc059-115.csv")
tsc2_3 <- read.csv("TSC2_exon5_65-uc059-115.csv")

intron_3 <- read.table("intron_barcode/intron_3ss_total")
intron_5 <- read.table("intron_barcode/intron_5ss_total")

setwd("~/Documents/yiting/little Hi/build oligo library/barcode_mt_3'5'ss/")

intron_3 <- intron_3[,-15]
intron_5 <- intron_5[,-c(5,6)]
intron <- bind_rows("3ss"=intron_3,"5ss"=intron_5, .id="3'/5'ss")

colnames(tsc1_5) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_exon","exon_end",
                        "exon_start","exon_name","ref_left_exon_length","ref_right_exon_length","exon_length",
                        "exon_in_seq","chr_extract","extract115_start","extract115_end","chr_extract115_start",
                      "zero","strand","wt_seq","extract115_left","extract115_right","mt_seq","len.left",
                      "mt_group_count","Name","F_primer","R_primer")

colnames(tsc1_3) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_exon","exon_end",
                      "exon_start","exon_name","ref_left_exon_length","ref_right_exon_length","exon_length",
                      "exon_in_seq","chr_extract","extract115_start","extract115_end","chr_extract115_start",
                      "zero","strand","wt_seq","extract115_left","extract115_right","mt_seq","len.right","one")

colnames(tsc2_5) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_exon","exon_start",
                      "exon_end","exon_name","ref_left_exon_length","ref_right_exon_length","exon_length",
                      "exon_in_seq","chr_extract","extract115_start","extract115_end","chr_extract115_start",
                      "zero","strand","wt_seq","extract115_left","extract115_right","mt_seq","len.left","one")

colnames(tsc2_3) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_exon","exon_start",
                      "exon_end","exon_name","ref_left_exon_length","ref_right_exon_length","exon_length",
                      "exon_in_seq","chr_extract","extract115_start","extract115_end","chr_extract115_start",
                      "extract_length","strand","wt_seq","extract115_left","extract115_right","mt_seq","len.right","one")


######intron_ref_colnames(tsc1) <- c("wt","variant","mt","wt/mt","ref_num","Name",
######                    "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
######                    "right","L_barcode_R","seq_len","exon")

exon <- bind_rows("tsc1_5ss"=tsc1_5,"tsc1_3ss"=tsc1_3,"tsc2_5ss"=tsc2_5,"tsc2_3ss"=tsc2_3, .id="id")
exon <- exon[,-c(2,9,17,21,28,33,34)]
exon$Name <- str_replace(exon$mt.name,pattern = ">",replacement = "to")

for (i in 1:nrow(exon)){
  if (str_detect(exon$id[i], pattern = fixed("3ss"))==TRUE){
    fprimer <- "TGGGGCCTGAAGACCTGACT";rprimer <- "GCTGCTGTGCTCCTTTTCCG"}
  else if (str_detect(exon$id[i], pattern = fixed("5ss"))==TRUE){
    fprimer <- "AGAGACAGGAGCCTTCGGTG";rprimer <- "AGGTAGGAGGATAGGCCAGG"}
  
  print(i)
  exon$F_primer[i] <- fprimer
  exon$R_primer[i] <- rprimer
}

exon$len_ref <- nchar(as.character(exon$ref))
exon$len_alt <- nchar(as.character(exon$alt))
exon$final_seq <- paste0(exon$F_primer,exon$mt_seq,exon$R_primer)
exon$final_length <- nchar(exon$final_seq)

##### add fasta info #########################################################
setwd("~/Documents/yiting/little Hi/new_library/reordered")
exon_5 <- read.table(file = "reordered_5_ss_exon_splicing.fasta")
fasta_exon5 <- data.frame(exon_5[seq(1,4313,by=2),1],exon_5[seq(2,4314,by=2),1])
colnames(fasta_exon5) <- c("fasta_name","fasta_seq")

exon_3 <- read.table(file = "reordered_3_ss_exon_splicing.fasta")
fasta_exon3 <- data.frame(exon_3[seq(1,3745,by=2),1],exon_3[seq(2,3746,by=2),1])
colnames(fasta_exon3) <- c("fasta_name","fasta_seq")

fasta_exon <- rbind(fasta_exon3,fasta_exon5)
fasta_exon$final_name <- str_replace(fasta_exon$fasta_name,pattern = ">",replacement = "")
fasta_exon$original_name <- str_replace(fasta_exon$final_name,pattern = "to",replacement = ">")
fasta_exon$original_name <- as.character(fasta_exon$original_name)
fasta_exon$original_name <- str_replace(fasta_exon$original_name,pattern = "-1$",replacement = "")
fasta_exon <- separate(fasta_exon,"original_name",into=c("wt","original_name"),sep="_")

fasta_exon[(fasta_exon$original_name=="")==TRUE,5] <- fasta_exon[(fasta_exon$original_name=="")==TRUE,3]


###### length(which(str_detect(fasta_exon$original_name,pattern = ">")))  3901

exon$name <- exon$mt.name
exon_final <- merge(fasta_exon,exon,by.y = "name",by.x = "original_name",all.x = TRUE)

#### CHECKING NUMBERS OF CANDIDATES COMPATIBILITY ###########
exon_both <- merge(fasta_exon,exon,by.y = "name",by.x = "original_name")
nrow(setdiff(exon_final,exon_both))
##### add 129 wt between exon_both and exon_final ##################################

length(setdiff(exon$name,fasta_exon$original_name))
setdiff(exon$name,fasta_exon$original_name)
###"chr16:2071873-2071873:->CTGCAGGCCCCGCCGT" "chr16:2075823-2075823:->GCCAGGCTGCCGCACCTCTA" "chr16:2081611-2081611:->GGAACACCAGCTGGCTG"
### may be sequences that are too long: all three of them have large insertion

exon_final_add_wt <- group_by(exon_final, wt) %>% 
  arrange(original_name, .by_group=TRUE) 


for (i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$id[i])){rowi <- i}
  exon_final_add_wt[rowi,c(7:34)] <- exon_final_add_wt[rowi+1,c(7:34)]
  print(rowi)
}

exon_final_add_wt$final_seq[is.na(exon_final_add_wt$id)] <- paste0(exon_final_add_wt$F_primer[is.na(exon_final_add_wt$id)],
                                                                   exon_final_add_wt$extract115_left[is.na(exon_final_add_wt$id)],
                                                                   exon_final_add_wt$ref[is.na(exon_final_add_wt$id)],
                                                                   exon_final_add_wt$extract115_right[is.na(exon_final_add_wt$id)],
                                                                   exon_final_add_wt$R_primer[is.na(exon_final_add_wt$id)])
exon_final_add_wt$final_length[is.na(exon_final_add_wt$id)] <- nchar(exon_final_add_wt$final_seq[is.na(exon_final_add_wt$id)])  

######### len.left comes from 5ss; len.right comes from 3ss
exon_final_add_wt$nchar <- nchar(as.character(exon_final_add_wt$fasta_seq))
exon_final_add_wt$"wt/mt_start(left)" <- NA
exon_final_add_wt$"wt/mt_end(right)" <- NA

for (i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$id[i])){
    mt_s <- NA ;mt_e <- NA}
  else if (str_detect(exon_final_add_wt$id[i], pattern = fixed("3ss"))==TRUE){
    mt_e <-exon_final_add_wt$nchar[i]-143-exon_final_add_wt$len.right[i] ;mt_s <- mt_e-exon_final_add_wt$len_alt[i]-1}
  else if (str_detect(exon_final_add_wt$id[i], pattern = fixed("5ss"))==TRUE){
    mt_s <- 40+exon_final_add_wt$len.left[i];mt_e <- mt_s+exon_final_add_wt$len_alt[i]+1}
  print(i)
  exon_final_add_wt$"wt/mt_start(left)"[i] <- mt_s
  exon_final_add_wt$"wt/mt_end(right)"[i] <- mt_e
}

for (i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$id[i])) {s <- "";e <- ""}
  else {s <- str_sub(exon_final_add_wt$fasta_seq[i],start = exon_final_add_wt$"wt/mt_start(left)"[i],end=exon_final_add_wt$"wt/mt_start(left)"[i]); 
    e <- str_sub(exon_final_add_wt$fasta_seq[i],start = exon_final_add_wt$"wt/mt_end(right)"[i],end=exon_final_add_wt$"wt/mt_end(right)"[i])}
  exon_final_add_wt$"wt_form_GT/AG"[i] <- str_replace(paste0(s,exon_final_add_wt$ref[i],e),pattern = "NA",replacement = "")
  exon_final_add_wt$"mt_form_GT/AG"[i] <- str_replace(paste0(s,exon_final_add_wt$alt[i],e),pattern = "NA",replacement = "")
  print(i); flush.console()
}

exon_final_add_wt$"wt_form_GT/AG"[which(is.na(exon_final_add_wt$`wt/mt_start(left)`))] <- ""
exon_final_add_wt$"mt_form_GT/AG"[which(is.na(exon_final_add_wt$`wt/mt_start(left)`))] <- ""

exon_final_add_wt$"wt_form_5_GT" <- NA
exon_final_add_wt$"wt_form_3_AG" <- NA

exon_final_add_wt$"mt_form_5_GT" <- NA
exon_final_add_wt$"mt_form_3_AG" <- NA

for (i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$id[i])){gt_wt <- NA; gt_mt <- NA}
  else if (str_detect(exon_final_add_wt$id[i],pattern = fixed("5ss"))==TRUE){
    gt_wt <- str_locate_all(exon_final_add_wt$`wt_form_GT/AG`[i],pattern = "GT")
    gt_mt <- str_locate_all(exon_final_add_wt$`mt_form_GT/AG`[i],pattern = "GT")
  }
  print(i)
  exon_final_add_wt$wt_form_5_GT[i] <- gt_wt
  exon_final_add_wt$mt_form_5_GT[i] <- gt_mt
}


for (i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$id[i])){ag_wt <- NA; ag_mt <- NA}
  else if (str_detect(exon_final_add_wt$id[i],pattern = "3ss")==TRUE){
    ag_wt <- str_locate_all(exon_final_add_wt$`wt_form_GT/AG`[i],pattern = fixed("AG"))
    ag_mt <- str_locate_all(exon_final_add_wt$`mt_form_GT/AG`[i],pattern = fixed("AG"))
  }
  print(i)
  exon_final_add_wt$"wt_form_3_AG"[i] <- ag_wt
  exon_final_add_wt$"mt_form_3_AG"[i] <- ag_mt
}




exon_final_add_wt$wt_form_5_GT_pos <- NA
for(i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$wt_form_5_GT[i])){pos <- NA}
  else pos <- sapply(exon_final_add_wt$wt_form_5_GT[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  exon_final_add_wt$wt_form_5_GT_pos[i] <- pos
}

exon_final_add_wt$wt_form_3_AG_pos <- NA
for(i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$wt_form_3_AG[i])){pos <- NA}
  else pos <- sapply(exon_final_add_wt$wt_form_3_AG[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  exon_final_add_wt$wt_form_3_AG_pos[i] <- pos
}

exon_final_add_wt <- data.frame(exon_final_add_wt)
exon_final_add_wt$wt_form_5_GT_pos <- as.character(exon_final_add_wt$wt_form_5_GT_pos)
exon_final_add_wt$wt_form_3_AG_pos <- as.character(exon_final_add_wt$wt_form_3_AG_pos)
which(str_detect(exon_final_add_wt$wt_form_5_GT_pos,pattern = ","))
which(str_detect(exon_final_add_wt$wt_form_3_AG_pos,pattern = ","))

exon_final_add_wt <- separate(exon_final_add_wt, "wt_form_5_GT_pos",into = c("wt_form_5_GT_pos1","wt_form_5_GT_pos2"),sep=",")
exon_final_add_wt <- separate(exon_final_add_wt, "wt_form_3_AG_pos",into = c("wt_form_3_AG_pos1","wt_form_3_AG_pos2"),sep=",")

exon_final_add_wt$"wt_form_5_GT_1_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$wt_form_5_GT_pos1)
exon_final_add_wt$"wt_form_5_GT_1_end" <- exon_final_add_wt$"wt_form_5_GT_1_start"+1
exon_final_add_wt$"wt_form_5_GT_2_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$wt_form_5_GT_pos2)
exon_final_add_wt$"wt_form_5_GT_2_end" <- exon_final_add_wt$"wt_form_5_GT_2_start"+1

exon_final_add_wt$"wt_form_3_AG_1_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$wt_form_3_AG_pos1)
exon_final_add_wt$"wt_form_3_AG_1_end" <- exon_final_add_wt$"wt_form_3_AG_1_start"+1
exon_final_add_wt$"wt_form_3_AG_2_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$wt_form_3_AG_pos2)
exon_final_add_wt$"wt_form_3_AG_2_end" <- exon_final_add_wt$"wt_form_3_AG_2_start"+1



exon_final_add_wt$mt_form_5_GT_pos <- NA
for(i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$mt_form_5_GT[i])){pos <- NA}
  else pos <- sapply(exon_final_add_wt$mt_form_5_GT[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  exon_final_add_wt$mt_form_5_GT_pos[i] <- pos
}

exon_final_add_wt$mt_form_3_AG_pos <- NA
for(i in 1:nrow(exon_final_add_wt)){
  if (is.na(exon_final_add_wt$mt_form_3_AG[i])){pos <- NA}
  else pos <- sapply(exon_final_add_wt$mt_form_3_AG[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  exon_final_add_wt$mt_form_3_AG_pos[i] <- pos
}



exon_final_add_wt <- data.frame(exon_final_add_wt)
exon_final_add_wt$mt_form_5_GT_pos <- as.character(exon_final_add_wt$mt_form_5_GT_pos)
exon_final_add_wt$mt_form_3_AG_pos <- as.character(exon_final_add_wt$mt_form_3_AG_pos)
which(str_detect(exon_final_add_wt$mt_form_5_GT_pos,pattern = ","))
which(str_detect(exon_final_add_wt$mt_form_3_AG_pos,pattern = ","))

exon_final_add_wt <- separate(exon_final_add_wt, "mt_form_5_GT_pos",into = c("mt_form_5_GT_pos1","mt_form_5_GT_pos2"),sep=",")
exon_final_add_wt <- separate(exon_final_add_wt, "mt_form_3_AG_pos",into = c("mt_form_3_AG_pos1","mt_form_3_AG_pos2"),sep=",")

exon_final_add_wt$"mt_form_5_GT_1_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$mt_form_5_GT_pos1)
exon_final_add_wt$"mt_form_5_GT_1_end" <- exon_final_add_wt$"mt_form_5_GT_1_start"+1
exon_final_add_wt$"mt_form_5_GT_2_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$mt_form_5_GT_pos2)
exon_final_add_wt$"mt_form_5_GT_2_end" <- exon_final_add_wt$"mt_form_5_GT_2_start"+1

exon_final_add_wt$"mt_form_3_AG_1_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$mt_form_3_AG_pos1)
exon_final_add_wt$"mt_form_3_AG_1_end" <- exon_final_add_wt$"mt_form_3_AG_1_start"+1
exon_final_add_wt$"mt_form_3_AG_2_start" <- exon_final_add_wt$wt.mt_start.left.-1 + as.numeric(exon_final_add_wt$mt_form_3_AG_pos2)
exon_final_add_wt$"mt_form_3_AG_2_end" <- exon_final_add_wt$"mt_form_3_AG_2_start"+1

exon_final_add_wt$"can_intron_start" <- NA
exon_final_add_wt$"can_intron_end" <- NA
exon_final_add_wt$"3/5ss" <- NA

############ check all exon kept same intron length #######
exon_final_add_wt$extract115_end - exon_final_add_wt$exon_end

exon_final_add_wt$id_wt <- NA

for (i in 1:nrow(exon_final_add_wt)){
  if ((is.na(exon_final_add_wt$id[i])==TRUE & is.na(exon_final_add_wt$len.left[i])==TRUE)==TRUE) {wt <- "3ss_wt"}
  else if ((is.na(exon_final_add_wt$id[i])==TRUE & is.na(exon_final_add_wt$len.right[i])==TRUE)==TRUE) {wt <- "5ss_wt"}
  else wt <- "mt"
  print(i)
  exon_final_add_wt$id_wt[i] <- wt
}


for (i in 1:nrow(exon_final_add_wt)){
  if (str_detect(exon_final_add_wt$id_wt[i], pattern = fixed("3ss_wt"))== TRUE) {in_s <- 28;in_e <- 132;ss <- 3}
  else if (str_detect(exon_final_add_wt$id_wt[i], pattern = fixed("5ss_wt"))== TRUE) {in_e <- exon_final_add_wt$nchar[i]-24;in_s <-in_e-173;ss <- 5}
  else if (str_detect(exon_final_add_wt$id[i], pattern = fixed("3ss"))== TRUE){
    in_s <- 28;in_e <- 132;ss <- 3}
  else if (str_detect(exon_final_add_wt$id[i], pattern = fixed("5ss"))== TRUE){
    in_e <- exon_final_add_wt$nchar[i]-24;in_s <- in_e-173;ss <- 5}
  print(i)
  exon_final_add_wt$"can_intron_start"[i] <- in_s
  exon_final_add_wt$"can_intron_end"[i] <- in_e
  exon_final_add_wt$"3/5ss"[i] <- ss
}

######### exclude list in dataframe for writing csv #################

exon_final_add_wt <- exon_final_add_wt[,-c(42,43,44,45)]
setwd("~/Documents/yiting/little Hi/build oligo library/barcode_mt_3'5'ss")
write.csv(exon_final_add_wt,file = "hi_3ss_5ss_exon_fixed2_final_add_GT_AG.csv")

##############################################################################
##################### intron_part ############################################
####### intron file name 3'/5'ss same as the distribution ##############################################################
#### only intronfar has assigned 3'/5'ss wrong: wrong when synthesizing oligos(fasta same as synthesized oligo) ######

setwd("~/Documents/yiting/little Hi/build oligo library")
in3_tsc1 <- read.csv(file = "TSC1_intron3_87-112.csv",row.names = 1)
in3_tsc2 <- read.csv(file = "TSC2_intron3_87-112-intrononly.csv",row.names = 1)
in5_tsc1 <- read.csv(file = "TSC1_intron5_87-112.csv",row.names = 1)
in5_tsc2 <- read.csv(file = "TSC2_intron5_87-112-intrononly.csv",row.names = 1)
inf3_tsc1 <- read_excel("TSC1_3__intron_far29.xlsx")
inf5_tsc1 <- read.csv(file = "TSC1_5__intron_far-29.csv",row.names = 1)
inf3_tsc2 <- read.csv(file = "TSC2_intron3_far-30-intrononly.csv",row.names = 1)
inf5_tsc2 <- read.csv(file = "TSC2_intron5_far-30-intrononly.csv",row.names = 1)


colnames(in3_tsc1) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_end",
                      "intron_start","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                      "intron_in_seq","chr_extract","extract112_start","extract112_end","chr_extract112_start",
                      "extract_length","strand","wt_seq_112","extract112_left","extract112_right","mt_seq_112")

colnames(in3_tsc2) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_start",
                        "intron_end","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                        "intron_in_seq","chr_extract","extract112_start","extract112_end","chr_extract112_start",
                        "extract_length","strand","wt_seq_115","wt_seq_112","extract112_left","extract112_right","mt_seq_112")

colnames(in5_tsc1) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_end",
                        "intron_start","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                        "intron_in_seq","chr_extract","extract112_start","extract112_end","chr_extract112_start",
                        "extract_length","strand","wt_seq_112","extract112_left","extract112_right","mt_seq_112","mt_group_count")

colnames(in5_tsc2) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_start",
                        "intron_end","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                        "intron_in_seq","chr_extract","extract112_start","extract112_end","chr_extract112_start",
                        "extract_length","strand","wt_seq_112","extract112_left","extract112_right","mt_seq_112")

colnames(inf3_tsc1) <- c("x","X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_end",
                        "intron_start","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                        "chr_extract","extract_29_start","extract_29_end","chr_extract29_start","intron_far_29","strand","wt_seq_29",
                        "extract29_left14","extract29_right14","mt_seq_29","exon20","intron63","wt_seq_112","mt_seq_112")
inf3_tsc1 <- inf3_tsc1[1:2,]
inf3_tsc1$intron_far_29 <- as.character(inf3_tsc1$intron_far_29)

colnames(inf5_tsc1) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_end",
                         "intron_start","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                         "chr_extract","extract_29_start","extract_29_end","chr_extract29_start","intron_far_29","strand",
                         "wt_seq_29","extract29_left14","extract29_right14","mt_seq_29","intron63","exon20","wt_seq_112","mt_seq_112")
inf5_tsc1$intron_far_29 <- as.character(inf5_tsc1$intron_far_29)

colnames(inf3_tsc2) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_start",
                         "intron_end","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                         "chr_extract","extract_29_start","extract_29_end","chr_extract29_start","zero","strand","wt_seq_29",
                         "extract29_left14","extract29_right14","mt_seq_29","exon20","intron63","wt_seq_112","mt_seq_112")

colnames(inf5_tsc2) <- c("X","chr","ref_start","ref_end","mt.name","ref","alt","chr_intron","intron_start",
                         "intron_end","intron_name","ref_left_intron_length","ref_right_intron_length","intron_length",
                         "chr_extract","extract_29_start","extract_29_end","chr_extract29_start","zero","strand",
                         "wt_seq_29","extract29_left14","extract29_right14","mt_seq_29","intron63","exon20","wt_seq_112","mt_seq_112")

intron_8file <- bind_rows("in3_tsc1"=in3_tsc1,"in3_tsc2"=in3_tsc2,"in5_tsc1"=in5_tsc1,"in5_tsc2"=in5_tsc2,
                     "inf3_tsc1"=inf3_tsc1,"inf3_tsc2"=inf3_tsc2,"inf5_tsc1"=inf5_tsc1,"inf5_tsc2"=inf5_tsc2, .id="id_intron_3'/5'_tsc1/2")

intron_8file <- intron_8file[,-c(2,30,40)]
intron_8file$len_extract112_left <- nchar(intron_8file$extract112_left)
intron_8file$len_extract112_right <- nchar(intron_8file$extract112_right)
intron_8file$len_ref <- nchar(intron_8file$ref)
intron_8file$len_alt <- nchar(intron_8file$alt)


for (i in 1:nrow(intron)){
  if (str_detect(intron$`3'/5'ss`[i], pattern = fixed("3"))==TRUE){
    extract112 <- str_sub(intron$final_seq[i],start = 21,end = intron$seq_len[i]-20-intron$Barcode_length[i])}
  else if (str_detect(intron$`3'/5'ss`[i], pattern = fixed("5"))==TRUE){
    extract112 <- str_sub(intron$final_seq[i],start = 21+intron$Barcode_length[i],end = intron$seq_len[i]-20)}
  print(i)
  intron$mt_seq_112[i] <- extract112
  intron$mt_seq_112[i] <- extract112
}

  
intron_final <- merge(intron_8file,intron,by.x = "mt_seq_112",by.y = "mt_seq_112",all.y = TRUE)
### which((intron_final$mt.name == intron_final$mt)==FALSE) integer(0)


###############################################################################
####### check numbers ########################################################
intron_both <- merge(intron_8file,intron,by.x = "mt_seq_112",by.y = "mt_seq_112")
added <- setdiff(intron_final,intron_both)
######## added 156 are all wt: correct #######################################
##############################################################################
intron_final <- intron_final[,-6]

###############################################################################
####################### add fasta file ########################################
setwd("~/Documents/yiting/little Hi/new_library/reordered")
fasta_in3 <- read.table(file = "reordered_3_ss_intron_splicing.fasta")
fasta_in3 <- data.frame(fasta_in3[seq(1,1303,by=2),1],fasta_in3[seq(2,1304,by=2),1])
colnames(fasta_in3) <- c("fasta_name","fasta_seq")

fasta_in5 <- read.table(file = "reordered_5_ss_intron_splicing.fasta")
fasta_in5 <- data.frame(fasta_in5[seq(1,1265,by=2),1],fasta_in5[seq(2,1266,by=2),1])
colnames(fasta_in5) <- c("fasta_name","fasta_seq")

fasta_intron <- rbind(fasta_in3,fasta_in5)
fasta_intron$nchar <- nchar(fasta_intron$fasta_seq)
fasta_intron$name <- str_replace(fasta_intron$fasta_name,pattern = ">",replacement = "")
fasta_intron$name <- str_replace(fasta_intron$name,pattern = "to",replacement = ">")
fasta_intron$name <- as.character(fasta_intron$name)
fasta_intron$name <- str_replace(fasta_intron$name,pattern = "-1$",replacement = "")

intron_total <- merge(intron_final,fasta_intron,by.x = "Name",by.y = "name")

################################## check ############################################
### nrow(intron_total[str_detect(intron_total$Name,pattern = ">"),])  mt:1129/ wt:156 correct
### length(which(is.na(intron_total$`id_intron_3'/5'_tsc1/2`)))  wt:156

intron_total$wt.mt_position_left <- NA
intron_total$wt.mt_position_right <- NA
intron_total$barcode_position_start <- NA
intron_total$barcode_position_end <- NA

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$`id_intron_3'/5'_tsc1/2`[i])==TRUE & (intron_total$`3'/5'ss`[i]=="3ss")==TRUE){
    mt_s <- NA ;mt_e <- NA; bc_e <- intron_total$nchar[i]-144;bc_s <- bc_e-intron_total$Barcode_length[i]+1}
  else if (is.na(intron_total$`id_intron_3'/5'_tsc1/2`[i])==TRUE & (intron_total$`3'/5'ss`[i]=="5ss")==TRUE){
    mt_s <- NA ;mt_e <- NA; bc_s <- 41;bc_e <- bc_s+intron_total$Barcode_length[i]-1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("in3_tsc1"))== TRUE){
    mt_s <- 87+intron_total$len_extract112_left[i] ;mt_e <- mt_s+intron_total$len_alt[i]+1;
    bc_e <- intron_total$nchar[i]-144;bc_s <- bc_e-intron_total$Barcode_length[i]+1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("in3_tsc2"))==TRUE){
    mt_s <- 87+intron_total$len_extract112_left[i] ;mt_e <- mt_s+intron_total$len_alt[i]+1;
    bc_e <- intron_total$nchar[i]-144;bc_s <- bc_e-intron_total$Barcode_length[i]+1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("in5_tsc1"))==TRUE){
    mt_e <- intron_total$nchar[i]-182-intron_total$len_extract112_right[i];mt_s <- mt_e-intron_total$len_alt[i]-1;
    bc_s <- 41;bc_e <- bc_s+intron_total$Barcode_length[i]-1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("in5_tsc2"))==TRUE){
    mt_e <- intron_total$nchar[i]-182-intron_total$len_extract112_right[i];mt_s <- mt_e-intron_total$len_alt[i]-1;
    bc_s <- 41;bc_e <- bc_s+intron_total$Barcode_length[i]-1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("inf3_tsc1"))==TRUE){
    mt_s <- 87+97;mt_e <- mt_s+intron_total$len_alt[i]+1;
    bc_e <- intron_total$nchar[i]-144;bc_s <- bc_e-intron_total$Barcode_length[i]+1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("inf3_tsc2"))==TRUE){
    mt_s <- 87+97;mt_e <- mt_s+intron_total$len_alt[i]+1;
    bc_e <- intron_total$nchar[i]-144;bc_s <- bc_e-intron_total$Barcode_length[i]+1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("inf5_tsc1"))==TRUE){
    mt_e <- intron_total$nchar[i]-183-96;mt_s <- mt_e-intron_total$len_alt[i]-1;
    bc_s <- 41;bc_e <- bc_s +intron_total$Barcode_length[i]-1}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("inf5_tsc2"))==TRUE){
    mt_e <- intron_total$nchar[i]-183-96;mt_s <- mt_e-intron_total$len_alt[i]-1;
    bc_s <- 41;bc_e <- bc_s +intron_total$Barcode_length[i]-1}
  print(i)
  intron_total$wt.mt_position_left[i] <- mt_s
  intron_total$wt.mt_position_right[i] <- mt_e
  intron_total$barcode_position_start[i] <- bc_s
  intron_total$barcode_position_end[i] <- bc_e
}

##################### check ############################
which((intron_total$alt == str_sub(intron_total$fasta_seq,start = intron_total$wt.mt_position_left+1,end = intron_total$wt.mt_position_right-1))==FALSE)
which((intron_total$Barcode == str_sub(intron_total$fasta_seq,start = intron_total$barcode_position_start,end = intron_total$barcode_position_end))==FALSE)
## different str_sub barcode 737 738 907 908 ( 738/908 TWO BARCODE SWITCHED??? ALL FROM INF5_TSC2, in check_12_file)
## mts from inf5_tsc2: which(str_detect(intron_total$`id_intron_3'/5'_tsc1/2`,pattern = "inf5_tsc2"))
#  21  89  91  93  95  97  99 172 224 300 302 356 440 568 669 738 908

intron_total$"wt_form_GT/AG" <- str_replace_all(paste0(str_sub(intron_total$fasta_seq,start = intron_total$wt.mt_position_left,end = intron_total$wt.mt_position_left),intron_total$ref,str_sub(intron_total$fasta_seq,start = intron_total$wt.mt_position_right,end = intron_total$wt.mt_position_right)),pattern = "NA",replacement = "")
intron_total$"mt_form_GT/AG" <- str_replace_all(paste0(str_sub(intron_total$fasta_seq,start = intron_total$wt.mt_position_left,end = intron_total$wt.mt_position_left),intron_total$alt,str_sub(intron_total$fasta_seq,start = intron_total$wt.mt_position_right,end = intron_total$wt.mt_position_right)),pattern = "NA",replacement = "")
intron_total$"barcode_form_GT/AG" <- str_replace_all(str_sub(intron_total$fasta_seq,start = intron_total$barcode_position_start-1,end = intron_total$barcode_position_end+1),pattern = "NA",replacement = "")


intron_total$wt_form_5_GT <- NA
intron_total$wt_form_3_AG <- NA
intron_total$mt_form_5_GT <- NA
intron_total$mt_form_3_AG <- NA
intron_total$barcode_form_5_GT <- NA
intron_total$barcode_form_3_AG <- NA

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$`id_intron_3'/5'_tsc1/2`[i])){
    ag_wt <- NA ;ag_mt <- NA ;ag_b <- NA; gt_wt <- NA ; gt_mt <- NA ;gt_b<- NA}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("3"))==TRUE){
    ag_wt <- str_locate_all(intron_total$"wt_form_GT/AG"[i],pattern = "AG");
    ag_mt <- str_locate_all(intron_total$"mt_form_GT/AG"[i],pattern = "AG");
    ag_b <- str_locate_all(intron_total$"barcode_form_GT/AG"[i],pattern = "AG");gt_wt <- NA;gt_mt <- NA;gt_b <- NA}
  else if (str_detect(intron_total$`id_intron_3'/5'_tsc1/2`[i], pattern = fixed("5"))==TRUE){
    gt_wt <- str_locate_all(intron_total$"wt_form_GT/AG"[i],pattern = "GT");
    gt_mt <- str_locate_all(intron_total$"mt_form_GT/AG"[i],pattern = "GT");
    gt_b <- str_locate_all(intron_total$"barcode_form_GT/AG"[i],pattern = "GT");ag_wt <- NA;ag_mt <- NA;ag_b <- NA}
  print(i)
  intron_total$wt_form_5_GT[i] <- gt_wt
  intron_total$mt_form_5_GT[i] <- gt_mt
  intron_total$barcode_form_5_GT[i] <- gt_b
  intron_total$wt_form_3_AG[i] <- ag_wt
  intron_total$mt_form_3_AG[i] <- ag_mt
    intron_total$barcode_form_3_AG[i] <- ag_b
}

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$wt_form_5_GT[i])){
    pos <- NA}
  else    pos <- sapply(intron_total$wt_form_5_GT[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$wt_form_5_GT_pos[i] <- pos 
}


for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$mt_form_5_GT[i])){
    pos <- NA}
  else    pos <- sapply(intron_total$mt_form_5_GT[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$mt_form_5_GT_pos[i] <- pos 
}

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$barcode_form_5_GT[i])){
    pos <- NA}
  else   pos <- sapply(intron_total$barcode_form_5_GT[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$barcode_form_5_GT_pos[i] <- pos 
}

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$wt_form_3_AG[i])){
    pos <- NA}
  else    pos <- sapply(intron_total$wt_form_3_AG[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$wt_form_3_AG_pos[i] <- pos 
}

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$mt_form_3_AG[i])){
    pos <- NA}
  else    pos <- sapply(intron_total$mt_form_3_AG[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$mt_form_3_AG_pos[i] <- pos 
}

for (i in 1:nrow(intron_total)){
  if (is.na(intron_total$barcode_form_3_AG[i])){
    pos <- NA}
  else    pos <- sapply(intron_total$barcode_form_3_AG[i],simplify = "vector",FUN = function(x)paste0(x[,1], collapse = ", "))
  print(i)
  intron_total$barcode_form_3_AG_pos[i] <- pos 
}


intron_total[which(nchar(intron_total$mt_form_5_GT_pos)>1),77:80]
intron_total[which(nchar(intron_total$barcode_form_5_GT_pos)>1),77:80]
intron_total[which(nchar(intron_total$barcode_form_3_AG_pos)>1),77:80]
############ mt/barcode_form_3_AG only one, no need to separate column #############
intron_total[which(nchar(intron_total$wt_form_5_GT_pos)>1),77:80]
intron_total[which(nchar(intron_total$wt_form_3_AG_pos)>1),77:80]
intron_total[which(nchar(intron_total$mt_form_3_AG_pos)>1),77:80]


intron_total <- separate(intron_total, "mt_form_5_GT_pos",into = c("mt_form_5_GT_pos1","mt_form_5_GT_pos2"),sep=",")
intron_total <- separate(intron_total,"barcode_form_5_GT_pos",into = c("barcode_form_5_GT_pos1","barcode_form_5_GT_pos2"),sep=",")
intron_total <- separate(intron_total,"barcode_form_3_AG_pos",into = c("barcode_form_3_AG_pos1","barcode_form_3_AG_pos2"),sep=",")

intron_total$wt_form_5_GT_start <- intron_total$wt.mt_position_left-1 + as.numeric(intron_total$`wt_form_5_GT_pos`)
intron_total$wt_form_5_GT_end <- intron_total$wt_form_5_GT_start+1

intron_total$mt_form_5_GT_1_start <- intron_total$wt.mt_position_left-1 + as.numeric(intron_total$`mt_form_5_GT_pos1`)
intron_total$mt_form_5_GT_1_end <- intron_total$mt_form_5_GT_1_start+1
intron_total$mt_form_5_GT_2_start <- intron_total$wt.mt_position_left-1 + as.numeric(intron_total$`mt_form_5_GT_pos2`)
intron_total$mt_form_5_GT_2_end <- intron_total$mt_form_5_GT_2_start+1

intron_total$barcode_form_5_GT_1_start <- intron_total$barcode_position_start-2 + as.numeric(intron_total$`barcode_form_5_GT_pos1`)
intron_total$barcode_form_5_GT_1_end <- intron_total$barcode_form_5_GT_1_start+1
intron_total$barcode_form_5_GT_2_start <- intron_total$barcode_position_start-2 + as.numeric(intron_total$`barcode_form_5_GT_pos2`)
intron_total$barcode_form_5_GT_2_end <- intron_total$barcode_form_5_GT_2_start+1

intron_total$wt_form_3_AG_start <- intron_total$wt.mt_position_left-1 + as.numeric(intron_total$`wt_form_3_AG_pos`)
intron_total$wt_form_3_AG_end <- intron_total$wt_form_3_AG_start+1

intron_total$mt_form_3_AG_start <- intron_total$wt.mt_position_left-1 + as.numeric(intron_total$`mt_form_3_AG_pos`)
intron_total$mt_form_3_AG_end <- intron_total$mt_form_3_AG_start+1

intron_total$barcode_form_3_AG_1_start <- intron_total$barcode_position_start-2 + as.numeric(intron_total$`barcode_form_3_AG_pos1`)
intron_total$barcode_form_3_AG_1_end <- intron_total$barcode_form_3_AG_1_start+1
intron_total$barcode_form_3_AG_2_start <- intron_total$barcode_position_start-2 + as.numeric(intron_total$`barcode_form_3_AG_pos2`)
intron_total$barcode_form_3_AG_2_end <- intron_total$barcode_form_3_AG_2_start+1



intron_total$"can_intron_start" <- NA
intron_total$"can_intron_end" <- NA

for (i in 1:nrow(intron_total)){
  if (str_detect(intron_total$`3'/5'ss`[i], pattern = fixed("3ss"))==TRUE){
    in_s <- 28;in_e <- intron_total$nchar[i]-164-intron_total$Barcode_length[i]}
  else if (str_detect(intron_total$`3'/5'ss`[i], pattern = fixed("5ss"))==TRUE){
    in_s <- 60+1+intron_total$Barcode_length[i];in_e <- intron_total$nchar[i]-24}
  print(i)
  intron_total$"can_intron_start"[i] <- in_s
  intron_total$"can_intron_end"[i] <- in_e
}

## CHECK  which((intron_total$can_intron_end-intron_total$can_intron_start+1 == intron_total$intron)==FALSE)
### 280 281 282 283 670 671 672
### all 7 are in3_tsc2(1wt-3snp,1wt-2snp), seq_len:154, 154-20-20-2(barcode_length)=112(mt_seq_length), 
### exonlength should be (112-92)+20+2(barcode_length)=42 but showed 43(nchar(t$exon) in build_3ss_intron_barcode.R) => hi previous file wrong
## CALCULATION CORRECT

############## remove list ##########################################################
intron_total <- intron_total[,-c(74:79)]
intron_total$new_name <- str_replace(intron_total$fasta_name,pattern = ">",replacement = "")
setwd("~/Documents/yiting/little Hi/build oligo library/barcode_mt_3'5'ss")
write.csv(intron_total,file = "hi_3ss_5ss_intron_fixed2_final_add_GT_AG.csv")


