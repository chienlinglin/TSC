library("dplyr")
setwd("~/Documents/little Hi/new_library")
ref <- read.table("3ss_new_reference(fisher)")
ex_mt3 <- read.table(file = "new_3_ss_exon_mut")

## separate library into exon mt/ intron mt
ex_ref <- merge(ref,ex_mt3,by.x="V2",by.y="x")
ex_ref <- data.frame(ex_ref[,2],ex_ref[,1],ex_ref[,3])
colnames(ex_ref) <- c("V1","V2","V3")
intron_ref <- setdiff(ref,ex_ref)

##replicate the mutant, duplicate the column
intron_ref$V4 = intron_ref$V3

##count character of ref>alt
library("tidyr")
intron_size <- separate(intron_ref,"V4", into=c("l","r"),sep="_chr")
intron_size <- separate(intron_size,"r", into=c("ref","alt"),sep="to")
View(intron_size)
library("stringr")
ref_nchar <- nchar(str_extract(intron_size[,5],"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)
alt_nchar <- nchar(str_extract(intron_size[,6],"[aA-zZ]+"), type = "chars", allowNA = FALSE, keepNA = NA)

intron_size$ref_nchar <- as.numeric(ref_nchar)
intron_size$alt_nchar <- as.numeric(alt_nchar)
intron_size[is.na(data.frame(intron_size))] <- 0
View(intron_size)

##calculate canonical intron size
intron_size$wt_can <- 152
intron_size$mt_can <- (intron_size$wt_can+intron_size$alt_nchar-intron_size$ref_nchar)

dat_intron <- data.frame(matrix(nrow = 1156,ncol = 2))
dat_intron[,1] <- as.character(rbind(as.character(intron_size[,3]),as.character(intron_size[,2])))
dat_intron[,2] <- as.numeric(rbind(intron_size[,10],intron_size[,9]))
dat_intron$X3 <- 27
dat_intron$X4 <- dat_intron$X2+dat_intron$X3
View(dat_intron)

##### intronic barcode
setwd("~/Documents/little Hi/new_library/reordered")
lib_bed <- read.table(file = "reordered_3ss_new_lib.bed")
lib_bed_i <- merge(lib_bed,dat_intron,by.x = "V1",by.y = "X1",all = FALSE)

setwd("~/Documents/little Hi/build oligo library/intron_barcode")

lib_bed_i$barcode_length <- lib_bed_i$V3 - (27+60+92+20+144) -(lib_bed_i$X2-152)
lib_bed_i$barcode_start <- lib_bed_i$X4 +20 +1
lib_bed_i$barcode_end <- lib_bed_i$X4 + 20 + barcode_length


dat_exon <- data.frame(matrix(nrow = 3616,ncol = 2))
dat_exon[,1] <- as.character(rbind(ex_ref[,3],ex_ref[,2]))
dat_exon[,2] <- 105 
dat_exon$X3 <- 27
dat_exon$X4 <- dat_exon[,2]+dat_exon$X3

dat <- distinct(rbind(dat_exon,dat_intron))
setwd("~/Documents/little Hi/home2_cllin_run3&run4/3ss/fisher")
write.table(dat,file = "3ss_canonical_ref")
