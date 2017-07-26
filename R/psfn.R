psfn <-
function(i,p,LS.resp,TrLS1,ngenes,nsubs){
  psv=matrix(NA,nsubs,p[i])  # genes selected for each group and gene set
  pcl=matrix(NA,nsubs,p[i])  # ps for each group and gene set
  for(si in 1:nsubs){
    mu1=colMeans(TrLS1[LS.resp==si,1:ngenes])
    mu2=colMeans(TrLS1[LS.resp!=si,1:ngenes])
    s1=apply(TrLS1[LS.resp==si,1:ngenes],2,sd)
    s2=apply(TrLS1[LS.resp!=si,1:ngenes],2,sd)
    PS=as.vector((mu1-mu2)/(s1+s2))         # PS statistic
    PS=abs(PS)
    psv[si,]=stat.gnames(PS,1:length(PS),crit=p[i])$gnames # selected genes
    pcl[si,]=stat.gnames(PS,1:length(PS),crit=p[i])$t      # corresponding PS for genes
  }
  ordr=order(c(pcl),decreasing=TRUE)    # ordering PS starting with the highest
  fset=c(psv)[ordr]                     # order genes according to their PS
  ffset=fset
  ftab=table(ffset)
  repg=names(ftab)[ftab>1]              # find replicated genes 
  if(length(repg)>0){
    for(rp in 1:length(repg)){
      wf=which(ffset==repg[rp])
      ffset=ffset[-1*wf[(2:(length(wf)))]]
    } 
  }
  ffset[1:p[i]]
}
