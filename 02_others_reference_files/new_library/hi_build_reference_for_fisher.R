setwd("~/Documents/little Hi/new_library")
intron_lib <- read.table("3ss_new_lib.bed" ,header=FALSE,stringsAsFactors= FALSE)

intron_lib <- data.frame(intron_lib$V1)

##filter out the mutant
library("dplyr")
intron_lib <- filter(intron_lib, grepl('to',intron_lib$intron_lib.V1))

##replicate the mutant, duplicate the column
intron_lib$mt = intron_lib$intron_lib.V1

##rename the wt column
colnames(intron_lib)[colnames(intron_lib) == "intron_lib.V1"] <- "wt"

##create the wt column(wt and mt colnames are the same before symbol"_-")
library("tidyr")
intron_lib <- separate(intron_lib,"wt", into=c("wt","a"),sep="_")
wt_ <- rep("_",times=2386)
intron_lib$wt_ <- wt_
wt <- paste0(intron_lib$wt,intron_lib$wt_)
t <- 1:2386
intron_lib <- data.frame(t,wt,intron_lib$mt)
colnames(intron_lib)[colnames(intron_lib) == "intron_lib.mt"] <- "mt"
View(intron_lib)
##write table
write.table(intron_lib, quote= FALSE,file = "3ss_new_reference(fisher)",sep="\t",row.names= FALSE,col.names = FALSE)




