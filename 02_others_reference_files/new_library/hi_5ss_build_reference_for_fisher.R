setwd("~/Documents/little Hi/new_library")
intron_lib <- read.table("5ss-1_new_lib.bed" ,header=FALSE,stringsAsFactors= FALSE)

intron_lib <- data.frame(intron_lib$V1)
View(intron_lib)
##separate mt and wt
library("dplyr")
mt <- filter(intron_lib, grepl('to',intron_lib$intron_lib.V1))
wt <- setdiff(intron_lib,mt) 
mt$mt2 <- mt$intron_lib.V1
wt$wt2 <- wt$intron_lib.V1

## use replicated columns to combine wt and mt
mt <- separate(mt,"intron_lib.V1",into = c("name1","V"),sep="_")
wt <- separate(wt,"intron_lib.V1",into = c("name2","V"),sep="_")
ref_lib <- merge(mt,wt,by.x="name1",by.y="name2",all.x=TRUE)

## add tag and build reference
t <- 1:2644
ref_lib2 <- data.frame(t,ref_lib$wt2,ref_lib$mt2)
View(ref_lib2)

##write table
write.table(ref_lib2, quote= FALSE,file = "5ss-1_new_reference(fisher)",sep="\t",row.names= FALSE,col.names = FALSE)




