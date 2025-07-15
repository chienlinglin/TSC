## Build GTF File  
###### 3'ss intron ######

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
lh_fasta <- read.table("reordered_3_ss_intron_splicing.fasta")
fasta <- data.frame(matrix(nrow = 652,ncol = 2))
fasta[,1] <- lh_fasta[seq(1,1303,by=2),1]
fasta[,2] <- lh_fasta[seq(2,1304,by=2),1]
fasta$X1 <- str_replace(fasta$X1,">","")

## Position Table
setwd("~/Documents/little Hi/build oligo library/intron_barcode")
intron_3ss <- read.table(file = "intron_3ss_total")
intron_3ss$Name <- str_replace(intron_3ss$Name,">","to")
intron_3ss <- merge(fasta,intron_3ss,by.x="X1",by.y="Name")

##check
which((nchar(intron_3ss$X2) == intron_3ss$total)==FALSE)

lh_range <- data.frame(gene_s = 1,
                       gene_e = intron_3ss$total,
                       mRNA_s = 1,
                       mRNA_e = intron_3ss$total,
                       exon_1_s = 1,
                       exon_1_e = 27,
                       intr_1_s = 28,
                       intr_1_e = 27+intron_3ss$intron,
                       exon_2_s = 28+intron_3ss$intron,
                       exon_2_e = intron_3ss$total,
                       LibSeq_s = 88,
                       LibSeq_e = 87+intron_3ss$seq_len-nchar(paste0(intron_3ss$F_priming,intron_3ss$R_priming)),
                       ID = sprintf("LH_TSC%04d", c(1:length(intron_3ss$X1))),
                       row.names = intron_3ss$X1,
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
write.table(pre_gff2, file = "3ss_intron_tmp.gtf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
tmp <- import.gff2("3ss_intron_tmp.gtf")

export.gff2(tmp, "3ss_intron.gtf")
