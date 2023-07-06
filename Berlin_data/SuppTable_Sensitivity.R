library(dplyr)
library(tidyr)
library(boa)
library(gridExtra)
library(tibble)
library(xtable)
library(gtable)
library(grid)
library(ggplot2)

str95HPD = function(vals) {
  str95HPDI = paste("[", formatC(boa.hpd(vals,0.05)[1],format="e"), ", ", formatC(boa.hpd(vals,0.05)[2],format="e"), "]", sep="")
  return(str95HPDI)
}
strMedian = function(vals) {
  return(formatC(median(vals,na.rm = T), format="e"))
}

calcSummaryStats = function(df, burnin, name1, name2) {
  len = dim(df)[1]
  lenB = burnin*len
  if (burnin != 0.0) {
    df = df[lenB:len,] 
  }
  df = df %>% select(-contains(c("treeLength","Serial", "Sample", "treeLikelihood", "treePriors", "treeHeight")))
  df1 = df %>% summarise(across(everything(), strMedian))
  df2 = df %>% summarise(across(everything(), str95HPD))
  res = rbind(df1,df2)
  row.names(res) = c(name1, name2)
  return(res)
}

setwd("Path/")
Main = read.table("Main.log", header=T)
SFixed = read.table("Sfixed.log", header=T)
DiffD = read.table("DiffD.log", header=T)
WoSus = read.table("WoSus.log", header=T)

FullData = rownames_to_column(as.data.frame(t(calcSummaryStats(DiffD, 0.1, "S1: median", "S1: 95% HPDI"))))
FullData = full_join(FullData, rownames_to_column(as.data.frame(t(calcSummaryStats(SFixed, 0.1, "S2: median", "S2: 95% HPDI")))), by="rowname")
FullData = full_join(FullData, rownames_to_column(as.data.frame(t(calcSummaryStats(WoSus, 0.1, "S3: median", "S3: 95% HPDI")))), by="rowname")
FullData = full_join(FullData, rownames_to_column(as.data.frame(t(calcSummaryStats(Main, 0.1, "median", "95% HPDI")))), by="rowname")

rownames(FullData) = FullData$rowname
FullData$rowname = NULL
rownames(FullData)[1:20] = c(expression(paste("posterior ")), expression(paste("likelihood ")), expression(paste("prior ")), expression(paste("transition/transversion bias ", kappa)), expression(paste("become non-infectious rate ", delta[Delta])), expression(paste("become non-infectious rate ", delta[Omicron])),
                   expression(paste("reproductive number ", R[paste(e,",",1)])), expression(paste("reproductive number ", R[paste(e,",",2)])), expression(paste("lambda ratio ", r[paste(lambda,",",Delta)])), expression(paste("lambda ratio ", r[paste(lambda,",",Omicron)])),
                   expression(paste("sampling proportion ", s[paste(Delta, ",", 1)])), expression(paste("sampling proportion ", s[paste(Delta, ",", 2)])), expression(paste("sampling proportion ", s[paste(Omicron, ",", 1)])), expression(paste("sampling proportion ", s[paste(Omicron, ",", 2)])),
                   expression(paste("frequency parameter ", pi[A])), expression(paste("frequency parameter ", pi[C])), expression(paste("frequency parameter ", pi[G])), expression(paste("frequency parameter ", pi[T])), expression(paste("gamma shape parameter")), expression(paste("clock rate")))

g=tableGrob(xtable(FullData,digits=3),theme = ttheme_default(base_size = 5, colhead=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill="gray100")), core=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill=matrix("gray100", nrow(FullData), ncol(FullData)))), rowhead=list(fg_params = list(parse=TRUE))))
g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 2, b = nrow(g), l = 1, r = ncol(g))
g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 1, l = 1, r = ncol(g))
ggsave(paste("SuppTable_Sensitivity",".pdf",sep=""), g, width=21, heigh=15, units="cm")