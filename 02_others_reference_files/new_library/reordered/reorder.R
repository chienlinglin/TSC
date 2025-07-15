library(dplyr)


###############################################################
########### 3ss ###############################################

setwd("~/Documents/little Hi/new_library")
s3 <- read.table(file = "3_ss_splicing.fasta")
length(which(str_detect(s3$V1,">")))
fasta <- data.frame(matrix(nrow=2525,ncol=2))

i <- seq(1,15145,by=6)
fasta[,1] <- s3[i,1]
fasta[,2] <- paste(s3[i+1,1],s3[i+2,1],s3[i+3,1],s3[i+4,1],s3[i+5,1],sep = "")
fasta$nchar <- nchar(fasta[,2])
fasta$name <- fasta[,1]
fasta <- separate(fasta,"name", into=c("l","r"),sep="_chr")  

fasta$l[is.na(fasta$r)] <- str_replace(fasta$l[is.na(fasta$r)],pattern = "_",replacement = "")
fasta$r[is.na(fasta$r)] <- 0

test <- group_by(fasta, l) %>%  
  arrange(desc(.$r), .by_group=TRUE) %>% 
  mutate(count = n()) %>%
  mutate(str_detect(dplyr::first(X1),"16:")) %>%
  ungroup()

max(test$count)
## 52

###### separate tsc1/2 by chromosome
fasta_3s_tsc1 <- dplyr::filter(test,`str_detect(dplyr::first(X1), "16:")`==FALSE)
fasta_3s_tsc2 <- dplyr::filter(test,`str_detect(dplyr::first(X1), "16:")`==TRUE)

fasta_3s_tsc1$r[fasta_3s_tsc1$r==0] <- Inf
fasta_3s_tsc1 <- group_by(fasta_3s_tsc1, l) %>%  
  arrange(r, .by_group=TRUE) %>% 
  ungroup()

final_3ss <- rbind(fasta_3s_tsc1,fasta_3s_tsc2)
final_3ss$X1 <- as.character(final_3ss$X1) 

##### 2 column into one column fasta file
fasta_3ss <- data.frame(matrix(nrow=5050,ncol=1))
fasta_3ss[seq(1, 5049, by = 2),1] <- final_3ss[,1]
fasta_3ss[seq(2, 5050, by = 2),1] <- final_3ss[,2]

##### build lib bed
start <- 1
lib_bed_3ss <- data.frame(final_3ss$X1,start,final_3ss$nchar)
lib_bed_3ss$final_3ss.X1 <- str_replace(lib_bed_3ss$final_3ss.X1,pattern = ">",replacement = "")

##### check reference for fisher test
ref_3s <- final_3ss[,c(1,4)]
ref_3s_mt <- filter(ref_3s,str_detect(X1,"to") == TRUE)
ref_3s_wt <- setdiff(ref_3s,ref_3s_mt)

write.table(fasta_3ss,file = "reordered_3_ss_splicing.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(lib_bed_3ss, file = "reordered_3ss_new_lib.bed",sep="\t",row.names= FALSE,col.names = FALSE,quote = FALSE)

##### separate 3'ss exon/ intron fasta for building gtf file
setwd("~/Documents/little Hi/new_library")
exon_mt_3 <- read.table(file = "new_3_ss_exon_mut")
final_3ss$x <- str_replace(final_3ss$X1,pattern = ">",replacement = "")
exon_mt_3_fasta <- merge(final_3ss,exon_mt_3, by.x="x",by.y="x")

##### 2 column into one column fasta file
fasta_3ss_exon <- data.frame(matrix(nrow=3746,ncol=1))
fasta_3ss_exon[seq(1, 3745, by = 2),1] <- exon_mt_3_fasta$X1
fasta_3ss_exon[seq(2, 3746, by = 2),1] <- exon_mt_3_fasta$X2
write.table(fasta_3ss_exon,file = "reordered/reordered_3_ss_exon_splicing.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)

fasta_3ss_intron <- data.frame(matrix(setdiff(fasta_3ss$matrix.nrow...5050..ncol...1.,fasta_3ss_exon$matrix.nrow...3746..ncol...1.),ncol = 1))
write.table(fasta_3ss_intron,file = "reordered/reordered_3_ss_intron_splicing.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)





###############################################################
########### 5ss  ##############################################

setwd("~/Documents/little Hi/new_library")
s5 <- read.table(file = "5ss-1_form_mod.fasta")
length(which(str_detect(s5$V1,">")))

fasta5 <- data.frame(matrix(nrow=2790,ncol=2))

i <- seq(1,16735,by=6)
fasta5[,1] <- s5[i,1]
fasta5[,2] <- paste(s5[i+1,1],s5[i+2,1],s5[i+3,1],s5[i+4,1],s5[i+5,1],sep = "")
fasta5$nchar <- nchar(fasta5[,2])
fasta5$name <- fasta5[,1]
fasta5 <- separate(fasta5,"name", into=c("l","r"),sep="_")  


test5 <- group_by(fasta5, l) %>%  
  arrange(desc(.$r), .by_group=TRUE) %>% 
  mutate(count = n()) %>%
  mutate(str_detect(dplyr::first(X1),"16:")) %>%
  ungroup()
  
max(test5$count)
## 68

###### separate tsc1/2 by chromosome
fasta_5s_tsc1 <- dplyr::filter(test5,`str_detect(dplyr::first(X1), "16:")`==FALSE)
fasta_5s_tsc2 <- dplyr::filter(test5,`str_detect(dplyr::first(X1), "16:")`==TRUE)

fasta_5s_tsc2$r[fasta_5s_tsc2$r==-1] <- Inf
fasta_5s_tsc2 <- group_by(fasta_5s_tsc2, l) %>%  
  arrange(r, .by_group=TRUE) %>% 
  ungroup()

final_5ss <- rbind(fasta_5s_tsc1,fasta_5s_tsc2)
final_5ss$X1 <- as.character(final_5ss$X1) 

##### 2 column into one column fasta file
fasta_5ss <- data.frame(matrix(nrow=5580,ncol=1))
fasta_5ss[seq(1, 5579, by = 2),1] <- final_5ss[,1]
fasta_5ss[seq(2, 5580, by = 2),1] <- final_5ss[,2]

##### build lib bed
lib_bed_5ss <- data.frame(final_5ss$X1,start,final_5ss$nchar)
lib_bed_5ss$final_5ss.X1 <- str_replace(lib_bed_5ss$final_5ss.X1,pattern = ">",replacement = "")

##### check reference for fisher test
ref_5s <- final_5ss[,c(1,4)]
ref_5s_mt <- filter(ref_5s,str_detect(X1,"to") == TRUE)
ref_5s_wt <- setdiff(ref_5s,ref_5s_mt)

write.table(fasta_5ss,file = "reordered_5ss-1_form_mod.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(lib_bed_5ss, quote= FALSE,file = "reordered_5ss-1_new_lib.bed",sep="\t",row.names= FALSE,col.names = FALSE)

##### separate 5'ss exon/ intron fasta for building gtf file
setwd("~/Documents/little Hi/new_library")
exon_mt_5 <- read.table(file = "new_5_ss_exon_mut")
final_5ss$x <- str_replace(final_5ss$X1,pattern = ">",replacement = "")
exon_mt_5_fasta <- merge(final_5ss,exon_mt_5, by.x="x",by.y="x")

##### 2 column into one column fasta file
fasta_5ss_exon <- data.frame(matrix(nrow=4314,ncol=1))
fasta_5ss_exon[seq(1, 4313, by = 2),1] <- exon_mt_5_fasta$X1
fasta_5ss_exon[seq(2, 4314, by = 2),1] <- exon_mt_5_fasta$X2
write.table(fasta_5ss_exon,file = "reordered/reordered_5_ss_exon_splicing.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)

fasta_5ss_intron <- data.frame(matrix(setdiff(fasta_5ss$matrix.nrow...5580..ncol...1.,fasta_5ss_exon$matrix.nrow...4314..ncol...1.),ncol = 1))
write.table(fasta_5ss_intron,file = "reordered/reordered_5_ss_intron_splicing.fasta",row.names = FALSE,col.names = FALSE,quote = FALSE)



