plotMCR <-
function(colg,ps.mean,ps.std,nrmethcl,nrmeth,p,np,nsamp){
  plot(ps.mean[,1:2],ylim=c(0,.6),xlim=c(min(p)-5,max(p)+5),xlab="Number of genes selected",
       ylab="Misclassification error rate",type="n",cex.sub=1.2,cex.axis=1.2,cex=1.2,cex.lab=1.2)
  for(k in 1:(nrmethcl*nrmeth)){
    est<-ps.mean[,k+1]
    std<-ps.std[,k+1]
    lower<-upper<-rep(0,np)
    for(i in 1:np){
      lower[i]<-est[i]-1.96*(std[i]/sqrt(nsamp))	
      upper[i]<-est[i]+1.96*(std[i]/sqrt(nsamp))
    }
    if(any(lower<0)){lower[lower<0]=0}
    if(any(upper>1)){lower[lower>1]=1}
    ly=1; lev=1; lw=1.5
    lines(p,est,lwd=lw,lty=ly,col=colg[k])
    for(j in 1:np){
      points(p[j],est[j],col=adjustcolor(colg[k],lev),pch=16)
      exp2<-rep(p[j],2)
      exp3<-c(lower[j],upper[j])
      exp4<-c((p[j]-.35),(p[j]+.35))
      exp5<-c(lower[j],lower[j])
      exp6<-c(upper[j],upper[j])
      lines(exp2,exp3,lwd=lw,lty=ly,col=adjustcolor(colg[k],lev))
      lines(exp4,exp5,lwd=lw,lty=ly,col=adjustcolor(colg[k],lev))
      lines(exp4,exp6,lwd=lw,lty=ly,col=adjustcolor(colg[k],lev))
    }
  }
  legend("topright",bty="n",legend = colnames(ps.mean)[-1],lwd=3,lty=1,col=colg,cex=1.2)
  box(lwd=3) 
}
