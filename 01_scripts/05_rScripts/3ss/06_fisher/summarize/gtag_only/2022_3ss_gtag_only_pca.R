
### PCA analysis
## option
options(stringsAsFactors = FALSE)

## Library
library(ggplot2)
library(ggfortify)
library(edgeR)

## Load data
setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/")
data.tbl_1 <- read.csv("3_1_most_can_non_gtag_only.csv",row.names = 1)
data.tbl_1 <- data.tbl_1[!rownames(data.tbl_1) %in% c( "2055273_.1" ) , ]
data.tbl_1 <- data.tbl_1[!rownames(data.tbl_1) %in% c( "2055273_chr16:2055287-2055288:TtoC.1" ) , ]
data.tbl_1[is.na(data.tbl_1)] <- 0
data.tbl_1[,10:13] <- data.tbl_1[,10:13] +1
data.tbl_1$t <- data.tbl_1$sum_spliced+data.tbl_1$s1_unspliced

data.tbl_2 <- read.csv("3_2_most_can_non_gtag_only.csv",row.names = 1)
data.tbl_2[is.na(data.tbl_2)] <- 0
data.tbl_2[,10:13] <- data.tbl_2[,10:13] +1
data.tbl_2$t <- data.tbl_2$sum_spliced+data.tbl_2$s2_unspliced

data.tbl_3 <- read.csv("3_3_most_can_non_gtag_only.csv",row.names = 1)
data.tbl_3[is.na(data.tbl_3)] <- 0
data.tbl_3[,10:13] <- data.tbl_3[,10:13] +1
data.tbl_3$t <- data.tbl_3$sum_spliced+data.tbl_3$s3_unspliced

## (canonical) (noncanonical) (unspliced) (spliced)
f1 <- data.frame(data.tbl_1$num_canonical/data.tbl_1$t,data.tbl_1$noncanonical/data.tbl_1$t,data.tbl_1$s1_unspliced/data.tbl_1$t,data.tbl_1$sum_spliced/data.tbl_1$t)
f2 <- data.frame(data.tbl_2$num_canonical/data.tbl_2$t,data.tbl_2$noncanonical/data.tbl_2$t,data.tbl_2$s2_unspliced/data.tbl_2$t,data.tbl_2$sum_spliced/data.tbl_2$t)
f3 <- data.frame(data.tbl_3$num_canonical/data.tbl_3$t,data.tbl_3$noncanonical/data.tbl_3$t,data.tbl_3$s3_unspliced/data.tbl_3$t,data.tbl_3$sum_spliced/data.tbl_3$t)
data.tbl_f <- cbind(f1,f2,f3)   

## PCA analysis
pca.data <- prcomp(t(cpm(data.tbl_f)))

## Visualize
# make group table
grp_tbl <- data.frame(sample=colnames(data.tbl_f),
                      group=rep(c("canonical","noncanonical","unspliced","spliced"),3))


#plot
pca_visual <- autoplot(object = pca.data,data=grp_tbl,colour="group") + theme_bw()

#save
setwd("/Users/angchu/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/3ss/fisher/summarize/")
ggsave("2022_3ss_gtag_only_pca.pdf",plot=pca_visual)

