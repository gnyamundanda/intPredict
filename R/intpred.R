# This is a pipeline of different gene selection and class prediction methods to identify a set of genes for predicting known classes
# edata: is the gene expression data with samples on the column.
# labs: is a vector with known classes (samples should be matched to the edata)
# seqp: is a vector with sequence of the number of to genes to selected
# fac: is a proportion how the data should be split for cross validation
# seeed: is a number used to fix results for reproducibility

intpred <- function(edata,labs,seqp,fac,seed=145)
{
  #### Required packages
  installed<-installed.packages()[,1] 
  required<-c("e1071","randomForest","rpart","MASS","siggenes",
              "multtest","samr","pamr","RColorBrewer",'devtools')
  toinstall<-required[!(required %in% installed)]
  if(length(toinstall) != 0){
    source("https://bioconductor.org/biocLite.R")
    biocLite(toinstall) }
  lapply(required, require, character.only = TRUE)
  requirfMM<-c("sma")
  toinstallfMM<-requirfMM[!(requirfMM %in% installed)]
  if(length(toinstallfMM) != 0){ install_github("gnyamundanda/sma") } 
  lapply(requirfMM, require, character.only = TRUE)
  
  #### GeneSelection folder for storing all results
  dir.create("GeneSelection", showWarnings = FALSE)
  setwd("GeneSelection")
  dir.create("indata", showWarnings = FALSE)
  dir.create("result", showWarnings = FALSE)
  dir.create("temp", showWarnings = FALSE)
  dir.create("topgenes", showWarnings = FALSE)
  pad='indata'    # training datasets storage
  padout='result' # results of gene selection and misclassification
  padout2='temp'  # temporary results storage
  
  #### prepare data for the classification problem
  set.seed(seed)                                       # set seed to remove randomness  due to sampling
  labels<-labs; length(labels)
  ord_lables<-order(labels)                            # order labels
  labs<-labels[ord_lables]
  expdata<-edata[,ord_lables]                          # order exp data according to ordered labels
  resp<-as.numeric(labs)
  nsubs=length(names(table(resp)))
  samplez<-paste(c("sample_"),1:length(labs),sep="") # create new naming for samples
  resp1<-as.data.frame(cbind(samplez,resp))
  genames<-paste(c("gene_"),1:nrow(expdata),sep="")    # create new naming for genes
  expd<-expdata
  rownames(expd)<-genames
  colnames(expd)<-samplez
  array<-as.data.frame(expd)                           # array is the expession data (genes on rows)
  
  #### setting up parameters
  nsamp=50                  # simulation size
  nrmethNames=c('pam','ps','bw')
  nrmeth=length(nrmethNames)	                # number of feature seletion methods
  nrmethclNames=c('rf','dlda','svmln','svmrd')
  nrmethcl=length(nrmethclNames)                # number of class prediction methods
  response=resp
  narray=length(response)
  dataset=trainTest(response, array, fac)   # generate training and testing data only to for setting up matrices  
  dataset1=as.data.frame(dataset$Learning_set)
  dataset2=as.data.frame(dataset$Test_set)
  nLarray=nrow(dataset1)                    # number of training samples (its fixed by fac)
  nTarray=nrow(dataset2)                    # number of test samples (its fixed by fac)
  test.idx=rep(0,narray)   # Matrix to store the sample selected as test data  for each of the nsamp runs
  p=seqp                                    # sequence of gene selection   
  np=length(p)                              # nr of gene sets                             
  ngenes=nrow(array)                        # nr of genes 
  unlistfn<-function(x,g,byrow) matrix(unlist(x),ncol=g,byrow=byrow)
  
  ######################
  #	1st. gene selection 
  ######################
  resgen=lapply(1:nsamp,function(k) gsetfn(response,array,fac,narray,np,p,nrmeth,ngenes,nsubs,test.idx))
  mat.train=unlistfn(lapply(resgen[1:nsamp],function(x){x['mat.tr'][[1]]}),nLarray,byrow=TRUE) 
  mat.test=unlistfn(lapply(resgen[1:nsamp],function(x){x['mat.ts'][[1]]}),nTarray,byrow=TRUE) 
  test.indexmatrix=unlistfn(lapply(resgen[1:nsamp],function(x){x['ts.idx'][[1]]}),narray,byrow=TRUE) 
  nameIndex=unlistfn(lapply(resgen[1:nsamp],function(x){x['nameIdx'][[1]]}),nTarray,byrow=TRUE) 
  Results=lapply(resgen[1:nsamp],function(x){x['res'][[1]]})
  
  ############## Storage of gene selection results
  # Storage matrix indicating which sample was selected as test data at each resampling: 1 identifies the test sample
  testIndex=paste(padout2,'/','testIndex','.','txt',sep='')
  write.table(test.indexmatrix, file=testIndex, append=FALSE, quote=F, sep=" ", eol="\n", na="NA", dec = ".", row.names = F, 
              col.names=c(1:narray), qmethod= c("escape"))  
  # Storage matrix with indexes samples selected in test data for each random sampling.
  nameIndex2 <- paste(padout2,'/','nameIndex','.','txt',sep='')
  write.table(nameIndex, file=nameIndex2, append=FALSE, quote=F, sep=" ", eol="\n", na="NA", dec=".", row.names=F, 
              col.names = c(1:nTarray), qmethod = c("escape"))
  # Matrix for storing the names of samples in the training set for each random sampling
  trainMatrix=paste(pad,'/','trainMatrix','.','txt',sep='')
  write.table( mat.train, file =  trainMatrix, append = FALSE, quote = F,sep = " ",eol = "\n", na = "NA", dec = ".", row.names =F , 
               col.names = c(1:nLarray), qmethod = c("escape"))
  # Matrix for storing the names of samples in the test set for each random sampling
  testMatrix=paste(pad,'/','testMatrix','.','txt',sep='')
  write.table( mat.test, file =  testMatrix, append = FALSE, quote = F,sep = " ",eol = "\n", na = "NA", dec = ".", row.names =F , 
               col.names = c(1:nTarray), qmethod = c("escape"))
  # Gene selection output: storing results of each gene selection output
  for(k in 1:nsamp){ 
    for(i in 1:np){ 
      gfile=paste(padout,'/','simulation_',k,'_top_',p[i],'.','txt',sep='')
      write.table(Results[[k]][[i]], file = gfile, append = FALSE, quote = F, sep = " ",
                  eol = "\n", na = "NA", dec = ".", row.names = F, col.names = TRUE, qmethod = c("escape"))}}
  filetemp=paste(padout,'/','feature_selection_worksp.RData',sep='')
  save.image(filetemp)
  #####
  
  #######################
  #	2nd. class prediction 
  #######################
  
  # storing number of errors and total number of samples considered per each simulation, each feature sel method in combn with class pred method
  method=list()
  for(l in 1:nrmeth){                 # running over number of feature selection methods
    method[[l]]=list()
    for(cp in 1:nrmethcl){ method[[l]][[cp]]=matrix(0,nsamp,np*2) } }
  for(k in 1:nsamp){  # run over number of simulations
    LS.resp=resp1[sapply(mat.train[k,], function(x){which(x==resp1[,1])} ),2]  # matching the samples
    TS.resp=resp1[sapply(mat.test[k,], function(x){which(x==resp1[,1])} ),2]
    TrLS=t(array[,sapply(mat.train[k,], function(x){which(x==resp1[,1])})])    # matching the samples to obtain exp data
    TrTS=t(array[,sapply(mat.test[k,], function(x){which(x==resp1[,1])})])
    LS=t(TrLS)
    TS=t(TrTS)
    TrLS1=data.frame(TrLS,LS.resp)
    TrTS1=data.frame(TrTS,TS.resp)
    for(i in 1:np){ # run over the number top p[i] features
      gfile=paste(padout,'/simulation_',k,'_top_',p[i],'.txt',sep='')  # find selected gene set for the top nr of genez (i) and sampling (k)
      Res2=read.table(file = gfile, header = TRUE, quote = "", sep = " ", na = "NA", dec = ".")
      for(l in 1:nrmeth){  # run over each feature selection method
        ## find the gene set selected by each feature selection method 
        LStop=Res2[,l]                              # select gene set for each gene selection method(l), top nr of genez (i) and sampling (k)
        TrLS2=data.frame(TrLS[,LStop],LS.resp)      # select training set associated with the gene set
        TrTS2=data.frame(TrTS[,LStop],TS.resp)      # select testing set
        
        #####
        ## RF
        LS.rf2<-randomForest(as.factor(LS.resp)~.,data=TrLS2,importance=TRUE,proximity=TRUE)
        TS.pred2<-predict(LS.rf2,TrTS2)               # predict the classes using RF  
        e3=table(TS.pred2,TS.resp)
        method[[l]][[1]][k,i]=sum(TS.pred2!=TS.resp)   # number of misclassified samples
        method[[l]][[1]][k,i+np]=sum(e3)                      # total number of samples
        
        ########
        ## DLDA 
        disc.dlda=stat.diag.da(TrLS[,LStop],TrTS[,LStop],cl=as.integer(LS.resp),pool=1)[[1]]
        f3=table(disc.dlda,TS.resp)
        method[[l]][[2]][k,i]=sum(disc.dlda!=TS.resp)
        method[[l]][[2]][k,i+np]=sum(f3)
        
        ##############
        ## SVM-Linear 
        model1=svm(as.factor(LS.resp)~.,data=TrLS2, kernel="linear")
        pred1=predict(model1, newdata=TrTS2)
        g3=table(pred1,TS.resp)
        method[[l]][[3]][k,i]=sum(pred1!=TS.resp)
        method[[l]][[3]][k,i+np]=sum(g3)
        
        ##############
        ## SVM-Radial 
        model1=svm(as.factor(LS.resp)~.,data=TrLS2, kernel="radial")
        pred1=predict(model1, newdata=TrTS2)
        g4=table(pred1,TS.resp)
        method[[l]][[4]][k,i]=sum(pred1!=TS.resp)
        method[[l]][[4]][k,i+np]=sum(g4)
      } # end l
    } # end i
    
    topcolnames<-c(paste0("Top",p),paste0("Top",p,"T"))
    cpm=paste0(nrmethclNames,'file')
    for(l in 1:nrmeth){ 
      for(cp in 1:nrmethcl){
        ffile=paste0(padout2,'/',cpm[cp],'-',nrmethNames[l],'.txt')
        write.table(method[[l]][[cp]], file=ffile, append=FALSE,quote=F, sep=" ",eol="\n",na="NA",
                    dec=".",row.names=F,col.names=topcolnames,qmethod=c("escape"))}}
  } # end k
  
  #####################
  #	Assess performance 
  #####################
  
  # Misclassification Error Rate (MCR)
  MCR=lapply(1:nrmeth,function(l){lapply(1:nrmethcl,function(m){method[[l]][[m]][1:nsamp,1:np]/nTarray })})
  # mean misclassification error rate for ALL data partition
  MCR2=lapply(1:nrmeth,function(l){lapply(1:nrmethcl,function(m){ colMeans(MCR[[l]][[m]]) })})
  # sd misclassification error rate for ALL data partition
  MCR3=lapply(1:nrmeth,function(l){lapply(1:nrmethcl,function(m){ apply(MCR[[l]][[m]],2,sd) })})
  
  ## Rearraing the MCR2 and MCR3 into a data frame of class pred methods, feature sel methods and the top feeatures
  MCRresult=MCRstd=matrix(0,nrmeth*nrmethcl,np+2)
  for(m in 1:nrmethcl){ for(l in 1:nrmeth){ 
    tm=(m-1)*nrmeth+l
    MCRresult[tm,1]=MCRstd[tm,1]=nrmethclNames[m]
    MCRresult[tm,2]=MCRstd[tm,2]=nrmethNames[l]
    MCRresult[tm,3:(np+2)]=round(MCR2[[l]][[m]][1:np],2)
    MCRstd[tm,3:(np+2)]=round(MCR3[[l]][[m]][1:np],2)
  }
  }
  MCRresult=data.frame(MCRresult)
  MCRstd=data.frame(MCRstd)
  filetemp2=paste(padout,'/','MCRmeangenes.txt',sep='')
  filetemp3=paste(padout,'/','MCRstdgenes.txt',sep='')
  testIndex=paste(padout2,'/','testIndex','.','txt',sep='')
  filetemp=paste(padout,'/','allworksp.RData',sep='')
  write.table(MCRresult,file=filetemp2,append=FALSE,quote=F,sep=" ",eol="\n",na="NA",dec=".",row.names=F,qmethod=c("escape"))
  write.table(MCRstd,file=filetemp3,append=FALSE,quote=F,sep=" ",eol="\n",na="NA",dec=".",row.names=F,qmethod=c("escape"))
  write.table(test.indexmatrix,file=testIndex,append=FALSE,quote=F,sep=" ",eol="\n",na="NA",dec=".",row.names=F,col.names=c(1:narray), qmethod=c("escape"))
  save.image(filetemp)
  
  ########
  #	Plots 
  #########
  
  # transform the average and std mcr into np by different classifier 
  means.all<-read.delim("result/MCRmeangenes.txt", header = TRUE, quote = "", sep = "", dec = ".")
  std.all<-read.delim("result/MCRstdgenes.txt", header = TRUE, quote = "", sep = "", dec = ".")
  order.means<-means.all[order(means.all$X2),]
  splitter<-order.means$X2              # split average mcr based on gene selection methods
  parts<-split(order.means,splitter)
  order.std<-std.all[order(std.all$X2),]
  partsd<-split(order.std,splitter)     # split std mcr based on gene selection methods
  cmean<-csd<-mnames<-NULL
  for(i in 1:length(parts)){    
    cmean<-cbind(cmean,matrix(as.numeric(t(parts[[i]])[-c(1,2),]),np,nrmethcl))
    csd<-cbind(csd,matrix(as.numeric(t(partsd[[i]])[-c(1,2),]),np,nrmethcl))
    mnames<-c(mnames,paste0(names(parts)[i],'_',nrmethclNames))}
  ps.mean<-cbind(p,cmean)
  ps.std<-cbind(p,csd)
  colnames(ps.mean)<-colnames(ps.std)<-c("Top",mnames)
  
  # mcr figure
  pdf(file=paste0(Sys.Date(),"_misclassification_error_rate.pdf"),onefile=TRUE,pointsize=16,width=10, height=10)
  par(mar=c(4,4,2,1))
  colg<-brewer.pal(n=12,name="Paired"); colg[11]<-"grey50"
  plotMCR(colg,ps.mean,ps.std,nrmethcl,nrmeth,p,np,nsamp)
  dev.off()
  gfile1=paste0(Sys.Date(),'_MCR.txt')
  gfile2=paste0(Sys.Date(),'_MCR_std.txt')
  write.table(ps.mean,file=gfile1,append=FALSE,quote=F,sep="\t",row.names=F,col.names=TRUE)
  write.table(ps.std,file=gfile2,append=FALSE,quote=F,sep="\t",row.names=F,col.names=TRUE)
  
  #################
  #	Selected genes 
  #################
  
  # final set of genes with the least MER
  # 1st. based on common genes among different gene selection methods at the with the optimal number if genes
  indx=rowSums(ps.mean[,-1]==min(ps.mean[,-1])) 
  idx=((1:nrow(ps.mean))[indx>0])[1]         # identify the smallest top nr of genes with least mcr
  msng<-matrix(NA,(nrow(array)-p[idx]),nsamp*nrmeth) # matrix for genes which were NOT selected by each gene sel metho at each sampling
  cnt=0
  for(k in 1:nsamp){     # each sampling
    gfile3=paste(padout,'/simulation_',k,'_top_',p[idx],'.txt',sep='')  
    for(j in 1:nrmeth){  # each gene sel method
      cnt=cnt+1
      ptop=read.table(file=gfile3, header=TRUE, quote="", sep=" ", na="NA", dec=".")
      msng[,cnt]<-c(1:nrow(array))[-ptop[,j]] # identify the index of genes which were NOT being selected
    }}
  sgenes<-sort(table(unlist(msng)),decreasing=TRUE)  # sort by genes which were frequently NOT selected
  rmgenes=as.numeric(names(sgenes)[1:(nrow(array)-p[idx])])  # select the top (total nr of genes - selected top nr of genes) genes were frequently NOT selected 
  topgenes<-rownames(expdata)[-rmgenes]
  gfile4=paste0('topgenes/',Sys.Date(),'_top_',p[idx],'_genes_based_on_most_commonly_selected_genes.txt')
  write.table(topgenes,file=gfile4,append=FALSE,quote=F,sep="",eol="\n",na="NA",dec=".",row.names=F,col.names=TRUE,qmethod=c("escape"))
  
  # 2nd. based on the best gene selection method
  nmdx=((colnames(ps.mean)[-1])[ps.mean[idx,-1]==min(ps.mean[idx,-1])])[1] # find best gene selection method
  mdx=strsplit(nmdx,'_')[[1]][1]               # name of the best gene selection method
  fsng<-matrix(NA,(nrow(array)-p[idx]),nsamp)  # matrix of genes which were NOT selected by best gene sel metho at each sampling
  for(k in 1:nsamp){ gfile5=paste(padout,'/','simulation_',k,'_top_',p[idx],'.','txt',sep='')  
  fptop=read.table(file=gfile5, header=TRUE, quote="", sep=" ", na="NA", dec=".")
  fsng[,k]<-c(1:nrow(array))[-fptop[,nrmethNames==mdx]] # identify those genes which are NOT being selected
  }
  msgenes<-sort(table(unlist(fsng)),decreasing=TRUE)  # sort by genes which were frequently NOT selected
  mrmgenes=as.numeric(names(msgenes)[1:(nrow(array)-p[idx])])  # select the top (p - best number of genes) genes were frequently NOT selected 
  mtopgenes<-rownames(expdata)[-mrmgenes]
  gfile6=paste0('topgenes/',Sys.Date(),'_top_',p[idx],'_genes_based_on_',mdx,'_gene_selection_method.txt')
  write.table(mtopgenes,file=gfile6,append=FALSE,quote=F,sep="",eol="\n",na="NA",dec=".",row.names=F,col.names=TRUE,qmethod=c("escape"))
  
  # 3rd. final set of genes used in the paper
  idx=19 # top 38 genes
  pmsng<-matrix(NA,(nrow(array)-p[idx]),nsamp*nrmeth)
  cnt=0
  for(k in 1:nsamp){ gfile7=paste(padout,'/','simulation_',k,'_top_',p[idx],'.','txt',sep='')
  for(j in 1:nrmeth){
    cnt=cnt+1
    pptop=read.table(file = gfile7, header = TRUE, quote = "", sep = " ", na = "NA", dec = ".")
    pmsng[,cnt]<-c(1:nrow(array))[-pptop[,j]]
  }}
  psgenes<-sort(table(unlist(pmsng)),decreasing = TRUE)
  prmgenes=as.numeric(names(psgenes)[1:(nrow(array)-p[idx])])
  ptopgenes<-rownames(expdata)[-prmgenes]
  gfile8=paste0('topgenes/',Sys.Date(),'_top_38_genes_based_on_most_commonly_selected_genes_used_in_the_paper.txt')
  write.table(ptopgenes,file=gfile8,append=FALSE,quote=F,sep="",eol="\n",na="NA",dec=".",row.names=F,col.names=TRUE,qmethod=c("escape"))
}
