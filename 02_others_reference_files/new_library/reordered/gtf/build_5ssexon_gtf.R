## Build GTF File  
###### 5'ss exon ######

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
lh_fasta <- readDNAStringSet("reordered_5_ss_exon_splicing.fasta")

## Position Table
lh_range <- data.frame(gene_s = 1,
                       gene_e = lh_fasta@ranges@width,
                       mRNA_s = 1,
                       mRNA_e = lh_fasta@ranges@width,
                       exon_1_s = 1,
                       exon_1_e = 40,
                       intr_1_s = 41,
                       intr_1_e = lh_fasta@ranges@width-24,
                       exon_2_s = lh_fasta@ranges@width-24+1,
                       exon_2_e = lh_fasta@ranges@width,
                       LibSeq_s = 41,
                       LibSeq_e = lh_fasta@ranges@width-(159+24),
                       ID = sprintf("LH_TSC%04d", c(1:length(lh_fasta@ranges@NAMES))),
                       row.names = lh_fasta@ranges@NAMES,
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

write.table(pre_gff2, file = "5ss_exon_tmp.gtf", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
tmp <- import.gff2("5ss_exon_tmp.gtf")

export.gff2(tmp, "5ss_exon.gtf")
