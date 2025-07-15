library(readxl)
library(dplyr)

setwd("~/Documents/little Hi/build oligo library/intron_barcode")
tsc1 <- read_excel("TSC1_intron5_87-112_code.xlsx")
tsc2 <- read_excel("TSC2_5_intron_code_ESE_checked.xlsx")
intron_far <- read_excel("TSC1_TSC2_5__intron_far_coding.xlsx")

## check no replication
which(duplicated(na.omit(c(tsc1$Name,tsc2$Name,intron_far$Name))))

colnames(tsc1) <- c("wt","variant","mt","chr","position","wt/mt","ref_num","Name",
                    "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                    "right","seq_len","exon")
colnames(tsc2)<- c("wt","variant","mt","wt/mt","ref_num","Name",
                   "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                   "right","seq_len","Barcode_length","exon")
colnames(intron_far)<- c("wt","variant","wt/mt","d1","ref_num","Name",
                   "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                   "right","seq_len","exon")
intron_far <- intron_far[,-c(4,5)]

tsc1 <- data.frame(apply(tsc1, 2, as.character))
tsc2 <- data.frame(apply(tsc2, 2, as.character))
intron_far <- data.frame(apply(intron_far, 2, as.character))

t <- bind_rows("tsc1"=tsc1,"tsc2"=tsc2,"tsc_intron_far"=intron_far, .id="tsc")

t$Barcode_length <- nchar(t$Barcode)
t$exon_length <- nchar(t$exon)
t$exon1 <- 20 + t$exon_length
t$intron <- as.numeric(t$seq_len) - t$exon_length + 139
t$exon2 <- 24
t$total <- t$exon1 + t$intron + t$exon2

write.table(t, file = "intron_5ss_total")
