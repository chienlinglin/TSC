library(readxl)
library(dplyr)

setwd("~/Documents/little Hi/build oligo library/intron_barcode")
tsc1 <- read_excel("TSC1_intron3_87-112_code.xlsx")
tsc2 <- read_excel("TSC2_3_intron_code.xlsx")
intron_far <- read_excel("TSC1_TSC2_3__intron_coding.xlsx")

## check no replication
which(duplicated(na.omit(c(tsc1$V4,tsc2$V4,intron_far$V4))))

colnames(tsc1) <- c("wt","variant","mt","wt/mt","ref_num","Name",
                    "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                    "right","L_barcode_R","seq_len","exon")
colnames(tsc2)<- c("wt","variant","mt","d1","ref_num","Name",
                   "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                   "right","seq_len","exon")
colnames(intron_far)<- c("wt","variant","mt","Name",
                         "Barcode","wt_alt_112","F_priming","R_priming","final_seq","ESE check left",
                         "right","seq_len","exon")
tsc2 <- tsc2[,-4]

tsc1 <- data.frame(apply(tsc1, 2, as.character))
tsc2 <- data.frame(apply(tsc2, 2, as.character))
intron_far <- data.frame(apply(intron_far, 2, as.character))

t <- bind_rows("tsc1"=tsc1,"tsc2"=tsc2,"tsc_intron_far"=intron_far, .id="tsc")

t$Barcode_length <- nchar(t$Barcode)
t$exon_length <- nchar(t$exon)
t$exon1 <- 27
t$intron <- 40 + as.numeric(t$seq_len) - t$exon_length
t$exon2 <- t$exon_length +124
t$total <- t$exon1 + t$intron + t$exon2

write.table(t, file = "intron_3ss_total")
