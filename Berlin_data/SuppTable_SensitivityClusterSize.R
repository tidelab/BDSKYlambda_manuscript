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
Main = read.table("Min4.log", header=T)
Min10 = read.table("Min10.log", header=T)
Min20 = read.table("Min20.log", header=T)

FullData = rownames_to_column(as.data.frame(t(calcSummaryStats(Min10, 0.1, "min10: median", "min10: 95% HPDI"))))
FullData = full_join(FullData, rownames_to_column(as.data.frame(t(calcSummaryStats(Min20, 0.1, "min20: median", "min20: 95% HPDI")))), by="rowname")
FullData = full_join(FullData, rownames_to_column(as.data.frame(t(calcSummaryStats(Main, 0.1, "min4: median", "min4: 95% HPDI")))), by="rowname")

rownames(FullData) = FullData$rowname
FullData$rowname = NULL
rownames(FullData)[1:19] = c(expression(paste("posterior ")), expression(paste("likelihood ")), expression(paste("prior ")), expression(paste("transition/transversion bias ", kappa)), expression(paste("become non-infectious rate ", delta)),
                   expression(paste("reproductive number ", R[paste(e,",",1)])), expression(paste("reproductive number ", R[paste(e,",",2)])), expression(paste("lambda ratio ", r[paste(lambda,",",Delta)])), expression(paste("lambda ratio ", r[paste(lambda,",",Omicron)])),
                   expression(paste("sampling proportion ", s[paste(Delta, ",", 1)])), expression(paste("sampling proportion ", s[paste(Delta, ",", 2)])), expression(paste("sampling proportion ", s[paste(Omicron, ",", 1)])), expression(paste("sampling proportion ", s[paste(Omicron, ",", 2)])),
                   expression(paste("frequency parameter ", pi[A])), expression(paste("frequency parameter ", pi[C])), expression(paste("frequency parameter ", pi[G])), expression(paste("frequency parameter ", pi[T])), expression(paste("gamma shape parameter")), expression(paste("clock rate")))

g=tableGrob(xtable(FullData,digits=3),theme = ttheme_default(base_size = 5, colhead=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill="gray100")), core=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill=matrix("gray100", nrow(FullData), ncol(FullData)))), rowhead=list(fg_params = list(parse=TRUE))))
g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 2, b = nrow(g), l = 1, r = ncol(g))
g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 1, l = 1, r = ncol(g))
ggsave(paste("SuppTable_SensitivityClusterSize",".pdf",sep=""), g, width=21, heigh=15, units="cm")
