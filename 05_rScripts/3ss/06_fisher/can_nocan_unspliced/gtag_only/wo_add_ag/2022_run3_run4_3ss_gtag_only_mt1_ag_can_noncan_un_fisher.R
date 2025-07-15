library(dplyr)
library(tidyverse)
options(stringsAsFactors = FALSE)
gtable<-read.table("~/Desktop/phD/Rotation/2/TSC alignment/R script/5ss/new_library/3ss_new_reference(fisher)", head=F)
View(gtable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
t1 <- read.csv(file="3_1_correct_can_non_gtag_only.csv",row.names = 1)
t2 <- read.csv(file="3_2_correct_can_non_gtag_only.csv",row.names = 1)
t3 <- read.csv(file="3_3_correct_can_non_gtag_only.csv",row.names = 1)


t5 <- read.csv(file="3_1_most_can_non_gtag_only.csv",row.names = 1)
t6 <- read.csv(file="3_2_most_can_non_gtag_only.csv",row.names = 1)
t7 <- read.csv(file="3_3_most_can_non_gtag_only.csv",row.names = 1)


s1_mt1 <- data.frame(t1[which(!is.na(t1$mt1_.ag_form)),1])        #152
s2_mt1 <- data.frame(t2[which(!is.na(t2$mt1_.ag_form)),1])        #147
s3_mt1 <- data.frame(t3[which(!is.na(t3$mt1_.ag_form)),1])        #152
t_mt1_3 <- data.frame(unique(c(s1_mt1[,1],s2_mt1[,1],s3_mt1[,1])))#190
length(intersect(intersect(s1_mt1[,1],s2_mt1[,1]),s3_mt1[,1]))    #112

s1_mt2 <- data.frame(t1[which(!is.na(t1$barcode_.ag_form)),1])    #26
s2_mt2 <- data.frame(t2[which(!is.na(t2$barcode_.ag_form)),1])    #23
s3_mt2 <- data.frame(t3[which(!is.na(t3$barcode_.ag_form)),1])    #29
t_mt2 <- data.frame(unique(c(s1_mt2[,1],s2_mt2[,1],s3_mt2[,1])))  #32
length(intersect(intersect(s1_mt2[,1],s2_mt2[,1]),s3_mt2[,1]))    #19

s5_mt1 <- data.frame(t5[which(!is.na(t5$mt1_.ag_form)),1])        #148
s6_mt1 <- data.frame(t6[which(!is.na(t6$mt1_.ag_form)),1])        #142
s7_mt1 <- data.frame(t7[which(!is.na(t7$mt1_.ag_form)),1])        #150
t_mt1_7 <- data.frame(unique(c(s5_mt1[,1],s6_mt1[,1],s7_mt1[,1])))#186
length(intersect(intersect(s5_mt1[,1],s6_mt1[,1]),s7_mt1[,1]))    #109

s5_mt2 <- data.frame(t5[which(!is.na(t5$barcode_.ag_form)),1])    #26
s6_mt2 <- data.frame(t6[which(!is.na(t6$barcode_.ag_form)),1])    #23
s7_mt2 <- data.frame(t7[which(!is.na(t7$barcode_.ag_form)),1])    #29
t_mt2_7 <- data.frame(unique(c(s5_mt2[,1],s6_mt2[,1],s7_mt2[,1])))#32
length(intersect(intersect(s5_mt2[,1],s6_mt2[,1]),s7_mt2[,1]))    #19


###canonical vs noncanonical
## correct
#################################################################################
############# t1 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t1[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t1[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t1[wtname, c(25,27)])
  wt_can <- t1[wtname, 25]
  wt_noncan_un <- t1[wtname, 27]
  mtsum <- sum(t1[mtname, c(25,27)])
  mt_can <- t1[mtname, 25]
  mt_noncan_un <- t1[mtname, 27]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_1_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_3_1_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t2 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t2[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t2[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t2[wtname, c(25,27)])
  wt_can <- t2[wtname, 25]
  wt_noncan_un <- t2[wtname, 27]
  mtsum <- sum(t2[mtname, c(25,27)])
  mt_can <- t2[mtname, 25]
  mt_noncan_un <- t2[mtname, 27]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_2_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_3_2_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t3 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t3[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t3[c(wtname,mtname),c(25,27)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t3[wtname, c(25,27)])
  wt_can <- t3[wtname, 25]
  wt_noncan_un <- t3[wtname, 27]
  mtsum <- sum(t3[mtname, c(25,27)])
  mt_can <- t3[mtname, 25]
  mt_noncan_un <- t3[mtname, 27]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_3_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_3_3_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()




############# most ##############################################################
#################################################################################
############# t5 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t5[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t5[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t5[wtname, c(22,24)])
  wt_can <- t5[wtname, 22]
  wt_noncan_un <- t5[wtname, 24]
  mtsum <- sum(t5[mtname, c(22,24)])
  mt_can <- t5[mtname, 22]
  mt_noncan_un <- t5[mtname, 24]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_1_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_most_3_1_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t6 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t6[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t6[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t6[wtname, c(22,24)])
  wt_can <- t6[wtname, 22]
  wt_noncan_un <- t6[wtname, 24]
  mtsum <- sum(t6[mtname, c(22,24)])
  mt_can <- t6[mtname, 22]
  mt_noncan_un <- t6[mtname, 24]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_2_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_most_3_2_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t7 ################################################################
##perform chisquare test
ptable=matrix(nrow=2386,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t7[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {pval<-NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {pval<-NA}
  else{pval<-chisq.test(ssi)$p.value}
  ptable[i,3] <- pval
  print(i)
}


##perform fisher test
for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t7[c(wtname,mtname),c(22,24)]
  ssi[is.na(ssi)] <- 0
  if (sum(ssi[mtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else if (sum(ssi[wtname, c(1,2)]) == 0) {fval<-NA;odd <- NA}
  else{fval<-fisher.test(ssi)$p.value
  odd<-fisher.test(ssi)$estimate}
  ptable[i,4] <- fval
  ptable[i,5] <- odd
  print(i)
}


## count
for (i in 1:nrow(gtable)){
  wtname<-as.character(gtable[i,2])
  mtname<-as.character(gtable[i,3])
  wtsum <- sum(t7[wtname, c(22,24)])
  wt_can <- t7[wtname, 22]
  wt_noncan_un <- t7[wtname, 24]
  mtsum <- sum(t7[mtname, c(22,24)])
  mt_can <- t7[mtname, 22]
  mt_noncan_un <- t7[mtname, 24]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_can
  ptable[i,8] <- wt_noncan_un
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_can
  ptable[i,11] <- mt_noncan_un
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_can","wt_noncan_un","wt_sum",
                      "mt_can","mt_noncan_un","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
write.csv(ptable, file="2022_run3_run4_3_3_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_most_3_3_gtag_only_mt1_ag_pval_fval_cu.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()




setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
f1 <- read.csv(file = "2022_run3_run4_3_1_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv",row.names = 1)
f2 <- read.csv(file = "2022_run3_run4_3_2_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv",row.names = 1)
f3 <- read.csv(file = "2022_run3_run4_3_3_gtag_only_mt1_ag_can_noncan_un_chisq_fisher.csv",row.names = 1)

library(qvalue)
f1_qval <- qvalue(f1$fval)
f1$FDR_qval <- f1_qval$qvalues
f2_qval <- qvalue(f2$fval)
f2$FDR_qval <- f2_qval$qvalues
f3_qval <- qvalue(f3$fval)
f3$FDR_qval <- f3_qval$qvalues


ft <- data.frame(f1[,1:12],f2[,3:12],f3[,3:12])
ft_sort <- ft[order(ft$fval,na.last = FALSE),]
pass_3_fval_strict <- (ft_sort$fval<0.05/2386) + (ft_sort$fval.1<0.05/2386)+(ft_sort$fval.2<0.05/2386)
ft_sort_f <- cbind(pass_3_fval_strict,ft_sort)
write.csv(ft_sort_f, file="2022_run3_run4_3ss_gtag_only_mt1_ag_can_noncan_un_sort_fvals.csv")



setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/can_noncan_unspliced/gtag_only/wo_add_ag/")
f1 <- read.csv(file = "2022_run3_run4_3_1_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv",row.names = 1)
f2 <- read.csv(file = "2022_run3_run4_3_2_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv",row.names = 1)
f3 <- read.csv(file = "2022_run3_run4_3_3_gtag_only_mt1_ag_most_can_noncan_un_chisq_fisher.csv",row.names = 1)


f1_qval <- qvalue(f1$fval)
f1$FDR_qval <- f1_qval$qvalues
f2_qval <- qvalue(f2$fval)
f2$FDR_qval <- f2_qval$qvalues
f3_qval <- qvalue(f3$fval)
f3$FDR_qval <- f3_qval$qvalues


ft <- data.frame(f1[,1:12],f2[,3:12],f3[,3:12])
ft_sort <- ft[order(ft$fval,na.last = FALSE),]
pass_3_fval_strict <- (ft_sort$fval<0.05/2386) + (ft_sort$fval.1<0.05/2386)+(ft_sort$fval.2<0.05/2386)
ft_sort_f <- cbind(pass_3_fval_strict,ft_sort)
write.csv(ft_sort_f, file="2022_run3_run4_3ss_gtag_only_mt1_ag_most_can_noncan_un_sort_fvals.csv")





can1 <- read.csv(file = "2022_run3_run4_3ss_gtag_only_mt1_ag_can_noncan_un_sort_fvals.csv")
can2 <- read.csv(file = "2022_run3_run4_3ss_gtag_only_mt1_ag_most_can_noncan_un_sort_fvals.csv")
can_final <- merge(can1,can2,by.x="X",by.y="X",all=TRUE)
can_final_sort <- can_final[order(can_final$fval.x,na.last = FALSE),]
colnames(can_final_sort) <- str_replace_all(colnames(can_final_sort),"x","can")
colnames(can_final_sort) <- str_replace_all(colnames(can_final_sort),"y","most")
can_final_sort_f <- can_final_sort[,(c(-1))]
write.csv(can_final_sort_f, file="2022_run3_run4_3ss_gtag_only_mt1_ag_final_can_noncan_un_sort_fvals.csv",row.names = can_final_sort$X)

