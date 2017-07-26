gsetfn <-
function(response,array,fac,narray,np,p,nrmeth,ngenes,nsubs,test.idx)
{
  ##########################################
  # resampling of training & testing data
  dataset<-trainTest(response,array,fac)        # call function to generate training and testing data  
  dataset1<-as.data.frame(dataset$Learning_set)  # training set
  dataset2<-as.data.frame(dataset$Test_set)      # test set
  mat.tr<-as.vector(row.names(dataset1))  # get names of samples in the training set 
  mat.ts<-as.vector(row.names(dataset2))   # get names of samples in the test set 
  nameIdx=sapply(row.names(dataset2),function(x) c(1:narray)[x==names(array)]) # identify the index of test samples
  test.idx[nameIdx]= 1          # samples used in the test data indicated by 1
  LS.resp=dataset1[,1]               # response data for training data
  TS.resp=dataset2[,1]               # response data for test data
  TrLS=dataset1[,-1]                 # exp data for training data
  TrTS=dataset2[,-1]                 # exp data for test data
  LS=t(TrLS)                         # exp training data in the form genes by samples
  TS=t(TrTS)                         # exp test data in the form genes by samples
  TrLS1=data.frame(TrLS,LS.resp)
  TrTS1=data.frame(TrTS,TS.resp)
  LS.resp2=LS.resp
  LS.resp3=LS.resp-1
  Genes=lapply(1:np,function(i) matrix(0,p[i],nrmeth))  # matrix to store the genes selected for each feature selection method
  nrm=1
  
  ######    
  # PAM
  sink(file="undesired_output.txt")
  data=list(x=LS,y=factor(LS.resp2), geneid=as.character(1:nrow(LS)), genenames=paste("g",as.character(1:nrow(LS)),sep=""))
  LS.pam1=pamr.train(data=data)
  pamres=lapply(1:np, function(i) pamfn(i,p,LS.pam1,data) ) 
  for(i in 1:np){ Genes[[i]][,nrm]<-pamres[[i]] }
  nrm=nrm+1
  sink()
  
  ########    
  #### PS 
  psres=lapply(1:np, function(i) psfn(i,p,LS.resp,TrLS1,ngenes,nsubs) ) 
  for(i in 1:np){ Genes[[i]][,nrm]<-psres[[i]] }
  nrm=nrm+1
  
  ##############
  #### BW-ratio 
  LS.bw=stat.bwss(LS,as.integer(LS.resp2))$bw
  for(i in 1:np){ Genes[[i]][,nrm]=stat.gnames(LS.bw,1:length(LS.bw),crit=p[i])$gnames }
  nrm=nrm+1
  
  list(mat.tr=mat.tr,mat.ts=mat.ts,nameIdx=nameIdx,ts.idx=test.idx,res=Genes)
}
