library("dplyr")
library("tidyr")
library("stringr")
setwd("~/Documents/little Hi/new_library")
ref <- read.table("5ss-1_new_reference(fisher)")
ex_mt5 <- read.table(file = "new_5_ss_exon_mut")

## separate library into exon mt/ intron mt
ex_ref <- merge(ref,ex_mt5,by.x="V2",by.y="x")
ex_ref <- data.frame(ex_ref[,2],ex_ref[,1],ex_ref[,3])
colnames(ex_ref) <- c("V1","V2","V3")
intron_ref <- setdiff(ref,ex_ref)

##replicate the mutant, duplicate the column
intron_ref$V4 = intron_ref$V3

##count character of ref>alt
intron_size <- separate(intron_ref,"V4", into=c("ref","alt"),sep="to")
View(intron_size)
in_ref_nchar <- nchar(str_extract(intron_size[,4],"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
in_alt_nchar <- nchar(str_extract(intron_size[,5],"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)

intron_size$in_ref_nchar <- as.numeric(in_ref_nchar)
intron_size$in_alt_nchar <- as.numeric(in_alt_nchar)
intron_size[is.na(data.frame(intron_size))] <- 0
View(intron_size)

##calculate canonical intron size
intron_size$wt_can <- 251
intron_size$mt_can <- (intron_size$wt_can+intron_size$in_alt_nchar-intron_size$in_ref_nchar)

dat_intron <- data.frame(matrix(nrow = 1102,ncol = 2))
dat_intron[,1] <- as.character(rbind(as.character(intron_size[,3]),as.character(intron_size[,2])))
dat_intron[,2] <- as.numeric(rbind(intron_size[,9],intron_size[,8]))

setwd("~/Documents/little Hi/new_library")
lib <- read.table(file = "5ss-1_new_lib.bed")
dat_intron <- merge(dat_intron, lib, by.x="X1" ,by.y= "V1")
dat_intron$X3 <- 60+ dat_intron$V3 -335 -(dat_intron$X2-(92+159))
dat_intron$X4 <- dat_intron$X2+dat_intron$X3
View(dat_intron)

##### intronic barcode
dat_intron$barcode_length <- dat_intron$V3 -335 -(dat_intron$X2-(92+159))
dat_intron$barcode_start <- 41
dat_intron$barcode_end <- 40 + barcode_length

## exon mutant canonical calculation
ex_ref$V4 = ex_ref$V3
exon_size <- separate(ex_ref,"V4", into=c("l","r"),sep="_chr")
exon_size <- separate(exon_size,"r", into=c("ref","alt"),sep="to")

ex_ref_nchar <- nchar(str_extract(exon_size[,5],"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)
ex_alt_nchar <- nchar(str_extract(exon_size[,6],"[A-Z]+"), type = "chars", allowNA = FALSE, keepNA = NA)

exon_size$ex_ref_nchar <- as.numeric(ex_ref_nchar)
exon_size$ex_alt_nchar <- as.numeric(ex_alt_nchar)
exon_size[is.na(data.frame(exon_size))] <- 0
View(exon_size)

exon_size$can_start <- 140 + exon_size$ex_alt_nchar - exon_size$ex_ref_nchar
exon_size$can_end <- 314 + exon_size$ex_alt_nchar - exon_size$ex_ref_nchar

dat_exon <- data.frame(matrix(nrow = 4186,ncol = 2))
dat_exon[,1] <- as.character(rbind(exon_size[,3],exon_size[,2]))
dat_exon[,2] <- 174
dat_exon$X3 <- exon_size$can_start
dat_exon$X4 <- exon_size$can_end
View(dat_exon)

dat_intron <- dat_intron[,-c(3,4)]
dat <- distinct(rbind(dat_exon,dat_intron))
setwd("~/Documents/little Hi/home2_cllin_run3&run4/5ss/wrong_spliced_fisher/")
write.table(dat,file = "new0602_5ss_canonical_ref")
