pamfn <-
function(i,p,LS.pam1,data)
{  # run this over the different number of features
  thres <- cbind(LS.pam1$threshold, LS.pam1$nonzero, LS.pam1$errors)
  thres=thres[(thres[,3]==(min(thres[((thres[,2]<=p[i]*3/2+20)&(thres[,2]>=p[i])),3])))&(thres[,2]<=p[i]*3/2+20)&(thres[,2]>=p[i]),1]  # use ALL criterias to identify diff thresholds
  thres=thres[length(thres)]     # take the highest threshold
  LS.p1=pamr.listgenes(LS.pam1,data=data,threshold=thres)     # pam centroids...genes are ranked by importance
  LS.p1=cbind(as.double(LS.p1[,1]),abs(as.double(LS.p1[,2])))
  LS.p1[1:p[i],1]
}
