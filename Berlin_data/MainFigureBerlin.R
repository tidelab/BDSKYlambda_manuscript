library(ape)
library(ggtree)
library(ggplot2)
library(gridExtra)
library(treeio)
library(hablar)
require(reshape2)
library(dplyr)
library(seqinr)
library(knitr)
library(phylotools)
library(tidyr)
library(grid)
library(utile.visuals)
library(boa)

makeSuppPlot = function(URLlinFile,URLTreeFolder, URLlogFIle, URLchangeDates, lastSamplingDate, GISLIN = TRUE) {

yeardays = c(20210101:20210131, 20210201:20210229, 20210301:20210331, 20210401:20210430, 20210501:20210531, 20210601:20210630, 20210701:20210731, 20210801:20210831, 20210901:20210930, 20211001:20211031, 20211101:20211130, 20211201:20211231, 20220101:20220131, 20220201:20220228)

latestSamplingDate = function(phylo) {
  yeardays = c(20210101:20210131, 20210201:20210229, 20210301:20210331, 20210401:20210430, 20210501:20210531, 20210601:20210630, 20210701:20210731, 20210801:20210831, 20210901:20210930, 20211001:20211031, 20211101:20211130, 20211201:20211231, 20220101:20220131, 20220201:20220228)
  labs = phylo$tip.label
  maxDate = 0
  for (lab in labs) {
    date = strsplit(lab, "_", fixed=TRUE)[[1]][9]
    date = strtoi(gsub("-", "", date))
    if (date > maxDate) {
      maxDate = date
    }
  }
  ratio = match(maxDate, yeardays)/366
  return(2021+ratio)
}

changeDates = read.table(URLchangeDates, header=FALSE, sep=",")
changeDates = changeDates[,1]
vLinePos = 2021 + match(changeDates, yeardays)/366
vLinePos = data.frame(int = 1:length(vLinePos), pos = vLinePos)
ints = length(changeDates) + 1


### LTT PLOT MCC TREES

setwd(URLTreeFolder)

fs = list.files(pattern = "\\.tree$")
treeb = read.beast(fs[1])
tree = as.phylo(treeb)
xy = ltt.plot.coords(tree)
ld = latestSamplingDate(tree)
xy[,1] = xy[,1] + ld
data_xy = data.frame(time = xy[,1], N=xy[,2])
lin1 = strsplit(strsplit(treeb@file, "cluster")[[1]][2], "tree")[[1]]
data2 = data_xy
data2$linsep = lin1
data2$lin = paste(lin1, " (", length(tree$tip.label), ")", sep="")


fs = fs[-1]
lins = c()
maxD = 0
maxO = 0
for (i in 1:length(fs)) {
  tb = read.beast(fs[i])
  lin = strsplit(strsplit(tb@file, "cluster")[[1]][2], "tree")[[1]]
  if (grepl("BA", lin, fixed=TRUE)) {
    maxO = max(maxO,max(tb@data$height_median, na.rm = TRUE))
  } else {
    maxD = max(maxD,max(tb@data$height_median, na.rm = TRUE))
  }
  lins = c(lins, strsplit(strsplit(tb@file, "cluster")[[1]][2], "tree")[[1]])
  t = as.phylo(tb)
  #ltt.lines(t, col = rainbow(length(fs))[i])
  xy = ltt.plot.coords(t)
  ld = latestSamplingDate(t)
  xy[,1] = xy[,1] + ld
  data_xy_help = data.frame(time = xy[,1], N=xy[,2])
  data_xy_help$linsep = lin
  data_xy_help$lin = paste(lin, " (", length(t$tip.label), ")", sep="")
  data2 = rbind(data2,data_xy_help)
}


lins = read.table(URLlinFile, header=TRUE, sep="\t")[,c(2,14)]
lins = data.frame(lins)
colnames(lins) = c("taxon", "lineage")

linCount = lins %>% group_by(lineage) %>% summarise(n())
names(linCount)[2] = "num"

data2 = data2 %>% mutate(VOC = ifelse(grepl("BA", linsep, fixed=TRUE), "Omicron", "Delta"))
cols = c("Delta" = "seagreen3", "Omicron" = "dodgerblue")
g1 = ggplot() + geom_line(data2, mapping=aes(time,N,col=VOC,group=lin), alpha = 0.7, show.legend = FALSE) + theme_bw() +
  labs(y="Number of tree branches") +
  geom_vline(xintercept = vLinePos$pos, alpha = 0.7, linetype="longdash", size=0.2)


### STEP PLOT R THROUGH TIME

monthRatios = c(2021 + c(0,31,31+29,31+29+31,31+29+31+30,31+29+31+30+31,31+29+31+30+31+30,31+29+31+30+31+30+31,31+29+31+30+31+30+31+31,31+29+31+30+31+30+31+31+30,31+29+31+30+31+30+31+31+30+31,31+29+31+30+31+30+31+31+30+31+30,366)/366)
monthLabs = c("Jan 21", "Feb 21", "Mar 21", "Apr 21", "May 21", "Jun 21", "Jul 21", "Aug 21", "Sep 21", "Oct 21", "Nov 21", "Dez 21", "Jan 22")

dataR = read.table(URLLogFile, header=TRUE)
print(dim(dataR))

lastSamplingInYears = 2021 + match(lastSamplingDate, yeardays)/366
print(lastSamplingInYears)
recent = lastSamplingInYears
print(recent)
changeDatesRelLastSamp = recent - rev(vLinePos$pos)
print(changeDatesRelLastSamp)
valsChange = c(0.0, changeDatesRelLastSamp, maxD + 0.03)
valsChangeO = c(0.0, changeDatesRelLastSamp, maxO)
F_times = rev(recent - valsChange)
F_timesO = rev(recent - valsChangeO)

RvalsD = data.frame(time=F_times, Rmedian=0, R5=0, R95=0)
RvalsO = data.frame(time=F_timesO, Rmedian=0, R5=0, R95=0)

hpds = matrix(ncol=2*ints,nrow=2)
medians = matrix(ncol=2*ints, nrow=1)
hpds2 = matrix(ncol=2*ints,nrow=2)
medians2 = matrix(ncol=2*ints, nrow=1)
for (i in 1:ints) {
  name = paste("reproductive.t.", i, sep="")
  print(head(dataR[[name]]))
  hpds[,i*2-1] = boa.hpd(dataR[[name]], 0.05)
  hpds[,i*2] = boa.hpd(dataR[[name]], 0.05) - 0.001
  medians[i*2-1] = median(dataR[[name]])
  medians[i*2] = median(dataR[[name]])
  RvalsD[i,2] = median(dataR[[name]])
  RvalsD[i,3] = hpds[,i*2-1][1]
  RvalsD[i,4] = hpds[,i*2-1][2]
  
  name = paste("lambda_ratio_o")
  hpds2[,i*2-1] = boa.hpd(dataR[[name]], 0.05)
  hpds2[,i*2] = boa.hpd(dataR[[name]], 0.05) - 0.001
  medians2[i*2-1] = medians[i*2] * median(dataR[[name]]) * 36.5/40.5
  medians2[i*2] = medians[i*2] * median(dataR[[name]])  * 36.5/40.5
  RvalsO[i,2] = medians[i*2] * median(dataR[[name]])  * 36.5/40.5
  RvalsO[i,3] = hpds[,i*2-1][1] * median(dataR[[name]]) * 36.5/40.5
  RvalsO[i,4] = hpds[,i*2-1][2] * median(dataR[[name]]) * 36.5/40.5
}

Rvals = full_join(RvalsD, RvalsO, by="time")
RvalsD["var"] = "Delta"
RvalsD[ints + 1,2:4] = RvalsD[ints,2:4]
RvalsO["var"] = "Omicron"
RvalsO[ints + 1,2:4] = RvalsO[ints,2:4]
Rvals2 = bind_rows(RvalsD,RvalsO)

Rvals[ints + 1,2:7] = Rvals[ints,2:7]

cols2 = c("Delta" = "darkseagreen3", "Omicron" = "deepskyblue3")
g2 = g1 + geom_step(data=Rvals2, mapping=aes(x=time, y=Rmedian*40, color=var), alpha = 1, size = 1.5) +
  geom_stepconfint(data=Rvals2, mapping=aes(x=time, ymin=R5*40, ymax=R95*40, fill=var), alpha=0.3) + 
  scale_color_manual(values = cols, guide="none") + scale_fill_manual(values = cols)  +
  scale_y_continuous(limits=c(0,100), sec.axis = sec_axis(trans=~./40, name=expression(bold("Re")))) +
  scale_x_continuous(name="Month", breaks = monthRatios, label = monthLabs) +
  theme(legend.position="bottom", legend.title = element_blank(), axis.text.y=element_text(face="bold")) +
  annotate(geom="text", x=2021.713, y=56, label="R[e]", parse=T, colour="seagreen3") +
  annotate(geom="text", x=2021.865, y=96.0, label="R[e]", parse=T, colour="dodgerblue")


ggsave(paste(URLTreeFolder, "/", "Berlin_MainFig.pdf", sep=""), g2, width=18.4, height=15, units="cm")
}
