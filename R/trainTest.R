trainTest <-
function(response, array,fac)
{
  tarray<-t(array)
  tarray2<-as.data.frame(tarray[order(response),])
  response2<-response[order(response)]
  nspl<-split(tarray2,response2)
  
  resp1<-nspl[[1]]
  resp2<-nspl[[2]]
  resp3<-nspl[[3]]
  resp4<-nspl[[4]]
  resp5<-nspl[[5]]
  
  
  nresp1<-nrow(resp1)
  nresp2<-nrow(resp2)
  nresp3<-nrow(resp3)
  nresp4<-nrow(resp4)
  nresp5<-nrow(resp5)
  
  nr1<-round(nresp1*fac,0)
  nr2<-round(nresp2*fac,0)
  nr3<-round(nresp3*fac,0)
  nr4<-round(nresp4*fac,0)
  nr5<-round(nresp5*fac,0)
  
  r1sample<-sample(c(1:nresp1),nr1)
  r2sample<-sample(c(1:nresp2),nr2)
  r3sample<-sample(c(1:nresp3),nr3)
  r4sample<-sample(c(1:nresp4),nr4)
  r5sample<-sample(c(1:nresp5),nr5)
  
  
  TR1<-resp1[r1sample,]
  TR2<-resp2[r2sample,]
  TR3<-resp3[r3sample,]
  TR4<-resp4[r4sample,]
  TR5<-resp5[r5sample,]
  
  LR1<-resp1[-r1sample,]
  LR2<-resp2[-r2sample,]
  LR3<-resp3[-r3sample,]
  LR4<-resp4[-r4sample,]
  LR5<-resp5[-r5sample,]
  
  TR<-rbind(TR1,TR2,TR3,TR4,TR5)
  LR<-rbind(LR1,LR2,LR3,LR4,LR5)
  cid<-unique(response2)
  cellType<-c(rep(cid[1],nr1),rep(cid[2],nr2),rep(cid[3],nr3),rep(cid[4],nr4),rep(cid[5],nr5))
  cellType2<-c( rep(cid[1],(nresp1-nr1)), rep(cid[2],(nresp2-nr2)), rep(cid[3],(nresp3-nr3)),rep(cid[4],(nresp4-nr4)),rep(cid[5],(nresp5-nr5)))
  TR2<-cbind(cellType,TR)
  LR2<-cbind(cellType2,LR)
  data<-list(TR2,LR2)
  names(data)<-c('Learning_set','Test_set')
  return(data)
}
