library(s20x)
library(boa)
library(xtable)
library(miscTools)
library(coda)
require(ggplot2)
require(gridExtra)
library(gtable)
library(grid)

colnameindex = function(M , colname0) 
{ 
  colsnames = names(M[1,]); 
  theindex = which(colsnames==colname0); 
  return(theindex); 
}


read = 0
scenario = 1; 
scenarios=c("Scenario9") 
#Scenario9: 10 transmission trees, N=25, five transmission trees with lambda_ratio_low = 0.8-1.2 and five with lambda_ratio_high = 1.4-1.8.

#lambda_ratio and samplingProportion inferred individually for all transmission trees.


L=c(87.6,98.55,109.5,120.45,131.4,153.3,164.25,175.2,186.15,197.1)
for(scenario in 1:1){
  
  print(scenarios[scenario])
  
  if (read==0){
    loglist = read.table(paste("logs",scenarios[scenario],".list",sep=''), as.is=TRUE, header=F)
    #loglist = read.table(paste("logs",".list",sep=''), as.is=TRUE, header=F)
    
    files= c()
    dim = min(100,length(loglist[,1]))
    ess=rep(NA,dim)
    esslimit =200
    samplearray=rep(1000000,dim)
    burninperc = 0.1	

#if scenario1 --> thissample,2:60
    
    for(i in 1:dim){
      assign(paste("log", i, sep=''), read.table(loglist[i,], header=T))
      thissample = length(get(paste("log", i, sep=''))[,1])
      ess1 = effectiveSize(get(paste("log", i, sep=''))[round(burninperc*thissample):thissample,2:60]) # check all likelihood and parameter ESSs 
      ess[i] = min(ess1[which(ess1>0)])
      if (ess[i]>esslimit) {
        files = c(files,i) 
        samplearray[i] = thissample
      }				
    }
  }
  
  
  if (scenario==1){Rb = c(3.0)
  s = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01)
  l=c(0.8,0.9,1.0,1.1,1.2,1.4,1.5,1.6,1.7,1.8)
  pars = 21 #length(names)
  c = colnameindex(log1,'baseReproductiveNumber.t')
  names = names(log1)[c]
  c = colnameindex(log1,'lambda_ratio_87.6')
  names = c(names,names(log1)[c:(c+9)])
  c = colnameindex(log1,'samplingProportion_87.6.t')
  names = c(names,names(log1)[c:(c+9)])
  }
  


  D = c(Rb,l,s)#,s/(d+s))#,m)
  #D = c(b,d,s,m)
  
  
  
  truth = t(matrix(data = D[1:pars], nrow = pars, ncol = dim))
  
  #means:
  bds=matrix(data = NA, nrow = dim, ncol = pars)
  #medians:
  bdsMed=matrix(data = NA, nrow = dim, ncol = pars)
  
  
  #hpd widths
  hpdWidths = matrix(data = NA, nrow = dim, ncol = pars)
  counts = matrix(data = 0, nrow = pars, ncol = 1)
  
  
  for(i in files){
    attach(get(paste("log", i, sep='')))
    samples = samplearray[i]
    burnin = round(burninperc*samples)
    
    for(j in 1:pars){
      current = get(names[j])
      bds[i,j]=mean(current[burnin:samples])
      bdsMed[i,j]=median(current[burnin:samples])
      hpd = boa.hpd(current[burnin:samples], 0.05)[1:2]
      if (truth[1,j] > hpd[1] && hpd[2] > truth[1,j]){
        counts[j] = counts[j] + 1
      }
      
      hpdWidths[i,j]=(hpd[2]-hpd[1])/truth[1,j]
    }
    detach(get(paste("log", i, sep='')))
  }
  
  #bias & error
  bias = (bdsMed - truth) / truth
  error = abs(bdsMed - truth) / truth
  
  
  plotBoolean = 1;
  if(plotBoolean){
    #plotBIAS
    pdf(paste("MedBiasPlot",scenarios[scenario], ".pdf", sep=""))
    par(mar=c(15,4,1,1))
    boxplot(bias, main = "Median relative bias", names=names, las=2)
    abline(0,0, col="green")
    dev.off()
    
    #plotERROR:
    pdf(paste("MedErrPlot", scenarios[scenario], ".pdf", sep = ""))
    par(mar=c(15,4,1,1))
    boxplot(error, main = "Median relative error", names=names, las=2)
    abline(0,0, col="green")
    dev.off()
  }
  
  
  
  print(paste("number of files with ESS>",esslimit,":",sep="")); print(length(files))
 
  
  if (scenario == 1) {
    paramnames=c(expression(R[base]),expression(r[paste(lambda,",",1)]),expression(r[paste(lambda,",",2)]),expression(r[paste(lambda,",",3)]),expression(r[paste(lambda,",",4)]),expression(r[paste(lambda,",",5)]),expression(r[paste(lambda,",",6)]),expression(r[paste(lambda,",",7)]),expression(r[paste(lambda,",",8)]),expression(r[paste(lambda,",",9)]),expression(r[paste(lambda,",",10)]),expression(s[1]),expression(s[2]),expression(s[3]),expression(s[4]),expression(s[5]),expression(s[6]),expression(s[7]),expression(s[8]),expression(s[9]),expression(s[10]))
  }  
  
  results = data.frame(c(truth[1,]), round(colMeans(bdsMed,na.rm=TRUE),4), round(colMeans(error,na.rm=TRUE),4), round(colMeans(bias,na.rm=TRUE), 4), round(colMeans(hpdWidths,na.rm=TRUE),4), round(counts/length(files)*100), row.names=paramnames)
  colnames(results)=c('Truth', 'Median', 'Rel. error', 'Rel. bias', 'Rel. HPD width', '95% HPD accuracy')
  write.csv(file=paste("logresults",scenarios[scenario],".csv",sep=''), results)
  
  print(xtable(results, digits=3),sanitize.rownames.function = identity)
  g=tableGrob(xtable(results, digits=3),theme = ttheme_default(base_size = 7, colhead=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill="gray100")), core=list(fg_params=list(hjust=1, x=0.95), bg_params=list(fill=matrix("gray100", nrow(results), ncol(results)))), rowhead=list(fg_params = list(parse=TRUE))))
  g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 2, b = nrow(g), l = 1, r = ncol(g))
  g=gtable_add_grob(g,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),t = 1, l = 1, r = ncol(g))
  ggsave(paste("Results", scenarios[scenario],".pdf",sep=""), g)
  ggsave(paste("Results", scenarios[scenario],".png",sep=""), g)
}
