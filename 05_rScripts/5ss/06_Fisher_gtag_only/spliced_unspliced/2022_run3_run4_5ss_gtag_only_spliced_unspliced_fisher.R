library(qvalue)

gtable<-read.table("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/5ss-1_new_reference(fisher)", head=F)
View(gtable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/")
t1 <- read.csv(file="5_1_correct_can_non_gtag_only.csv",row.names = 1)
t2 <- read.csv(file="5_2_correct_can_non_gtag_only.csv",row.names = 1)
t3 <- read.csv(file="5_3_correct_can_non_gtag_only.csv",row.names = 1)



###spliced vs unspliced

#################################################################################
############# t1 ################################################################
##perform chisquare test
ptable=matrix(nrow=2644,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t1[c(wtname,mtname),c(14,12)]
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
  ssi <- t1[c(wtname,mtname),c(14,12)]
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
  wtsum <- sum(t1[wtname, c(14,12)])
  wt_spliced <- t1[wtname, 14]
  wt_unspliced <- t1[wtname, 12]
  mtsum <- sum(t1[mtname, c(14,12)])
  mt_spliced <- t1[mtname, 14]
  mt_unspliced <- t1[mtname, 12]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_spliced
  ptable[i,8] <- wt_unspliced
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_spliced
  ptable[i,11] <- mt_unspliced
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_spliced","wt_unspliced","wt_sum",
                      "mt_spliced","mt_unspliced","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/spliced_unspliced/gtag_only/")
write.csv(ptable, file="2022_run3_run4_5_1_gtag_only_spliced_unspliced_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_5_1_gtag_only_pval_fval_su.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t2 ################################################################
##perform chisquare test
ptable=matrix(nrow=2644,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t2[c(wtname,mtname),c(14,12)]
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
  ssi <- t2[c(wtname,mtname),c(14,12)]
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
  wtsum <- sum(t2[wtname, c(14,12)])
  wt_spliced <- t2[wtname, 14]
  wt_unspliced <- t2[wtname, 12]
  mtsum <- sum(t2[mtname, c(14,12)])
  mt_spliced <- t2[mtname, 14]
  mt_unspliced <- t2[mtname, 12]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_spliced
  ptable[i,8] <- wt_unspliced
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_spliced
  ptable[i,11] <- mt_unspliced
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_spliced","wt_unspliced","wt_sum",
                      "mt_spliced","mt_unspliced","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/spliced_unspliced/gtag_only/")
write.csv(ptable, file="2022_run3_run4_5_2_gtag_only_spliced_unspliced_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_5_2_gtag_only_pval_fval_su.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()



#################################################################################
############# t3 ################################################################
##perform chisquare test
ptable=matrix(nrow=2644,ncol=11)
ptable=data.frame(ptable)

for (i in 1:nrow(gtable)){
  wtname <- as.character(gtable[i,2])
  mtname <- as.character(gtable[i,3])
  ssi <- t3[c(wtname,mtname),c(14,12)]
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
  ssi <- t3[c(wtname,mtname),c(14,12)]
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
  wtsum <- sum(t3[wtname, c(14,12)])
  wt_spliced <- t3[wtname, 14]
  wt_unspliced <- t3[wtname, 12]
  mtsum <- sum(t3[mtname, c(14,12)])
  mt_spliced <- t3[mtname, 14]
  mt_unspliced <- t3[mtname, 12]
  ptable[i,1] <- wtname
  ptable[i,2] <- mtname
  ptable[i,7] <- wt_spliced
  ptable[i,8] <- wt_unspliced
  ptable[i,9] <- wtsum
  ptable[i,10] <- mt_spliced
  ptable[i,11] <- mt_unspliced
  ptable[i,12] <- mtsum
  print(i); flush.console()
  
}

#qobj<-qvalue(p= !is.na(ptable[,4]), fdr.level=0.1, lambda=0, pi0.method="smoother")
#summary(qobj)


##change column names
colnames(ptable) <- c("wt","mt","pval","fval","odds_ratio","FDR_qval","wt_spliced","wt_unspliced","wt_sum",
                      "mt_spliced","mt_unspliced","mt_sum")
View(ptable)

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/spliced_unspliced/gtag_only/")
write.csv(ptable, file="2022_run3_run4_5_3_gtag_only_spliced_unspliced_chisq_fisher.csv")

##draw histogram of pvalues and fvalues
png(file="hist_5_3_gtag_only_pval_fval_su.png")
par(mfrow=c(1,2))
hist(ptable$pval,breaks = 20)
hist(ptable$fval,breaks = 20)
dev.off()




setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/fisher/output/spliced_unspliced/gtag_only/")
f1 <- read.csv(file = "2022_run3_run4_5_1_gtag_only_spliced_unspliced_chisq_fisher.csv",row.names = 1)
f2 <- read.csv(file = "2022_run3_run4_5_2_gtag_only_spliced_unspliced_chisq_fisher.csv",row.names = 1)
f3 <- read.csv(file = "2022_run3_run4_5_3_gtag_only_spliced_unspliced_chisq_fisher.csv",row.names = 1)


f1_qval <- qvalue(f1$fval)
f1$FDR_qval <- f1_qval$qvalues
f2_qval <- qvalue(f2$fval)
f2$FDR_qval <- f2_qval$qvalues
f3_qval <- qvalue(f3$fval)
f3$FDR_qval <- f3_qval$qvalues


ft <- data.frame(f1[,1:12],f2[,3:12],f3[,3:12])
ft_sort <- ft[order(ft$fval,na.last = FALSE),]
pass_3_fval_strict <- (ft_sort$fval<0.05/2644) + (ft_sort$fval.1<0.05/2644)+(ft_sort$fval.2<0.05/2644)
ft_sort_f <- cbind(pass_3_fval_strict,ft_sort)

write.csv(ft_sort_f, file="2022_run3_run4_5ss_gtag_only_spliced_unspliced_sort_fvals.csv")

