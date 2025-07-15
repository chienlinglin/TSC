# calculate total read counts of each oligos

setwd("~/Desktop/phD/Rotation/2/TSC alignment/re_aligment/20220318_modified_with with --mp10:4 -k 1 noncaonal_ss 1000/5ss/junction")

tbl_prc <- read.table(file = "2022_5ss_q60_total.count")

colnames(tbl_prc) <- c("chr","start","end","5-1_S12","5-1_S4","5-2_S13","5-2_S5","5-3_S14","5-3_S6")

tbl_prc$"5-1" <- tbl_prc$`5-1_S12` + tbl_prc$`5-1_S4`
tbl_prc$"5-2" <- tbl_prc$`5-2_S13` + tbl_prc$`5-2_S5`
tbl_prc$"5-3" <- tbl_prc$`5-3_S14`+ tbl_prc$`5-3_S6`

#png(filename = "2022_5ss_align_hist_t_count.png")
pdf(file = "2022_5ss_align_hist_t_count.pdf")
par(mfrow=c(3,3))
hist(tbl_prc$`5-1`,breaks = 2790)
hist(tbl_prc$`5-2`,breaks = 2790)
hist(tbl_prc$`5-3`,breaks = 2790)
hist(log10(tbl_prc$`5-1`+1),breaks = 2790)
hist(log10(tbl_prc$`5-2`+1),breaks = 2790)
hist(log10(tbl_prc$`5-3`+1),breaks = 2790)
hist(log10(tbl_prc$`5-1`+1),breaks = 10)
hist(log10(tbl_prc$`5-2`+1),breaks = 10)
hist(log10(tbl_prc$`5-3`+1),breaks = 10)
dev.off()

## Function
#Cor_plot
pairs_upper<-function(x, y){
  par(new=TRUE)
  plot(x,y, pch=20, type = "p", log = "xy", cex = .1, axes = FALSE)
  r <- round(cor(x, y, method = "spearman"), digits=3)
  txt <- paste0("rho : ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(2, 8.2, txt, col = ifelse(r >=0.8, "#38b48b", ifelse(r < 0.6, "#c9171e", "#0094c8")), cex = 0.5/strwidth(txt), font = 2)
}

##options
options(stringsAsFactors = FALSE)

##Pairs plot
##load data
row.names(tbl_prc) <- tbl_prc$chr
tbl <- tbl_prc[,c(10:12)]

View(tbl)
##output
pdf(file= "2022_5ss_align_pairplot.pdf")
pairs(tbl+1, upper.panel = pairs_upper, lower.panel = NULL, main = "2022_5ss_align_pairplot", log = "xy")
dev.off()

######## sessionInfo()
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS 15.4.1
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] zh_TW.UTF-8/zh_TW.UTF-8/zh_TW.UTF-8/C/zh_TW.UTF-8/zh_TW.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] compiler_4.2.0 tools_4.2.0   