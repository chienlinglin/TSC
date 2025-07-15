library("Biostrings")

setwd("~/Documents/little Hi/new_library")
test <- read.table("exon_far_non-splice_form_mod2.fasta",header=FALSE,stringsAsFactors= FALSE)
View(test)

## paste sequence by every 5 rows
chromosome <- data.frame(matrix(nrow=1860,ncol=10))
for (i in seq(1,16732,by=9)){
  chrom <- as.character(test[i,1])
  seq1 <- nchar(as.character(test[i+1,1]))
  seq2 <- nchar(as.character(test[i+2,1]))
  seq3 <- nchar(as.character(test[i+3,1]))
  seq4 <- nchar(as.character(test[i+4,1]))
  seq5 <- nchar(as.character(test[i+5,1]))
  seq6 <- nchar(as.character(test[i+6,1]))
  seq7 <- nchar(as.character(test[i+7,1]))
  seq8 <- nchar(as.character(test[i+8,1]))
  chromosome[(i+8)/9,1] <- chrom
  chromosome[(i+8)/9,2] <- seq1
  chromosome[(i+8)/9,3] <- seq2
  chromosome[(i+8)/9,4] <- seq3
  chromosome[(i+8)/9,5] <- seq4
  chromosome[(i+8)/9,6] <- seq5
  chromosome[(i+8)/9,7] <- seq6
  chromosome[(i+8)/9,8] <- seq7
  chromosome[(i+8)/9,9] <- seq8
  chromosome[(i+8)/9,10] <- sum(c(seq1,seq2,seq3,seq4,seq5,seq6,seq7,seq8))
}

View(chromosome)
lib_bed <- data.frame(chromosome$X1,1,chromosome$X10)
View(lib_bed)
lib_bed <- separate(lib_bed,"chromosome.X1",sep=">",into = c("sign","chromosome"))
lib_bed_f <- data.frame(lib_bed[,2],lib_bed[,3],lib_bed[,4])
colnames(lib_bed_f) <- c("chromosome","start","end")
View(lib_bed_f)  
write.table(lib_bed_f, quote= FALSE,file = "Exon_far_non-splice_new_lib.bed",sep="\t",row.names= FALSE,col.names = FALSE)



## yunlin's methos to reverse complement the sequences
library("Biostrings")
setwd("~/Documents/Hung Lun/library(intron,exon)/new pasted")
temp <- readDNAStringSet(filepath = "hunglun_BP_intron.fasta")
temp_rv <- reverseComplement(temp)
writeXStringSet(temp_rv,"hunglun_comrev_biostrings.fasta")
