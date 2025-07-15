## Build GTF File  
###### 5'ss intron ######

#### Global
## function
# GTF Attribute
lh_gff_attr <- function(range_tbl, nm, typ = c("gene", "mRNA", "exon", "LibSeq")){
  ID <- range_tbl[["ID"]]
  num <- strsplit(nm, "_")[[1]][2]
  attr <- switch(typ, 
                 gene = {return(sprintf("gene_id %s;", ID))},
                 mRNA = {return(sprintf("gene_id %s;transcript_id %s;", ID, paste0("m", ID)))},
                 exon = {return(sprintf("gene_id %s;transcript_id %s;", paste0("m", ID), paste0("m", ID, ":", num)))},
                 LibSeq = {return(sprintf("gene_id %s;transcript_id %s;", paste0("m", ID), paste0("m", ID, ":LibSeq")))}
  )
  return(attr)
}





## library
library(Biostrings)
library(rtracklayer)


## option
options(stringsAsFactors = FALSE)


#### Work
## Load FASTA
setwd("~/Documents/little Hi/new_library/reordered")
lh_fasta <- read.table("reordered_5_ss_intron_splicing.fasta")
fasta <- data.frame(matrix(nrow = 633,ncol = 2))
fasta[,1] <- lh_fasta[seq(1,1265,by=2),1]
fasta[,2] <- lh_fasta[seq(2,1266,by=2),1]
fasta$X1 <- str_replace(fasta$X1,">","")
fasta$nchar <- nchar(fasta$X2)
fasta$lib <- str_sub(fasta$X2,start=21,end = fasta$nchar-(139+24))

## Position Table
setwd("~/Documents/little Hi/build oligo library/intron_barcode")
intron_5ss <- read.table(file = "intron_5ss_total")
intron_5ss$Name <- str_replace(intron_5ss$Name,">","to")
intron_5ss$Name <- str_replace(intron_5ss$Name,">","to")

intron_5ss <- merge(fasta,intron_5ss,by.x="lib",by.y="final_seq")

##check
which((intron_5ss$nchar == intron_5ss$total)==FALSE)

lh_range <- data.frame(gene_s = 1,
                       gene_e = intron_5ss$total,
                       mRNA_s = 1,
                       mRNA_e = intron_5ss$total,
                       exon_1_s = 1,
                       exon_1_e = intron_5ss$exon1,
                       intr_1_s = intron_5ss$exon1+1,
                       intr_1_e = intron_5ss$total-24,
                       exon_2_s = intron_5ss$total-23,
                       exon_2_e = intron_5ss$total,
                       LibSeq_s = 41,
                       LibSeq_e = intron_5ss$total-(159+24),
                       ID = sprintf("LH_TSC%04d", c(1:length(intron_5ss$X1))),
                       row.names = intron_5ss$X1,
                       stringsAsFactors = FALSE)

## Build gff2 like data.frame
pre_gff2 <- list()
for(nm in c("gene", "mRNA", "exon_1", "exon_2", "LibSeq")){
  pre_gff2[[nm]] <- data.frame(row.names(lh_range), 
                               "LH_Lib",
                               strsplit(nm, "_")[[1]][1],
                               lh_range[[paste0(nm, "_s")]],
                               lh_range[[paste0(nm, "_e")]],
                               ".",
                               "+",
                               ".",
                               lh_gff_attr(lh_range, nm = nm, typ = strsplit(nm, "_")[[1]][1]),
                               stringsAsFactors = FALSE
  )
}

pre_gff2 <- do.call(rbind, pre_gff2)

#pre_gff2 <- pre_gff2[order(pre_gff2[[1]]),]

setwd("~/Documents/little Hi/new_library/reordered")
write.table(pre_gff2, file = "5ss_intron_tmp.gtf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
tmp <- import.gff2("5ss_intron_tmp.gtf")

export.gff2(tmp, "5ss_intron.gtf")
