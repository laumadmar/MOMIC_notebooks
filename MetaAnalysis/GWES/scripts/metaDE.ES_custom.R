##########################################################################################
###### Modified version of METADE.ES from github metaOmics/MetaDE                 ########
###### This version corrects for all NA in a row when the estimator (var and ES)  ########
###### of a gene is not present in every single datasets combined.                ########
###### This version adds na.rm=TRUE to all calculations of rowSums, the result is ########
###### having also an estimator when the gene is not present in all datasets.     ########
###### Author: Laura Madrid  May 2018                                             ########
##########################################################################################

MetaDE.ES.custom<-function(x,meta.method=c("FEM","REM")){
    meta.method<-match.arg(meta.method)
    K<-ncol(x$ES)
    if(meta.method=="REM"){
        res<-get.REM(x$ES,x$Var,pe=x$perm.ES,pv=x$perm.Var)
        tempFDR<-matrix(res$FDR,ncol=1)
        rownames(tempFDR)<-rownames(x$ES)
        colnames(tempFDR)<-meta.method
        meta.res<-list(n=res$n,mu.hat=res$mu.hat,mu.var=res$mu.var,Qval=res$Qval,Qpval=res$Qpval,tau2=res$tau2,zval=res$zval,pval=res$pval,FDR=tempFDR)	
    }else{
        res<-get.FEM(x$ES,x$Var,x$perm.ES,x$perm.Var)
        tempFDR<-matrix(res$FDR,ncol=1)
        rownames(tempFDR)<-rownames(x$ES)
        colnames(tempFDR)<-meta.method
        meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,zval=res$zval,pval=res$pval,FDR=tempFDR)	
    }
    attr(meta.res,"nstudy")<-K
    attr(meta.res,"meta.method")<-meta.method 
    class(meta.res)<-"MetaDE.ES"
    return(meta.res)
}

get.REM<-function(em,vm,pe=NULL,pv=NULL){
  k<-ncol(em)
  n<-k-rowSums(is.na(em))
  Q.val<-get.Q(em,vm)
  tau2<-get.tau2(Q.val,vm,k)
  temp.wt<-1/(vm+tau2)
  mu.hat<-rowSums(temp.wt*em,na.rm=TRUE)/rowSums(temp.wt,na.rm=TRUE)
  mu.var<-1/rowSums(temp.wt,na.rm=TRUE)
  Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
  z.score<-get.REM2(em,vm)$zval
  if(!is.null(pe)&!is.null(pv)){
    rnum<-which(apply(em,1,function(x) !any(is.na(x))))
    Z0<-matrix(get.REM2(pe,pv)$zval,nrow(em),nrow(pe)/nrow(em))
    z.p<-rep(NA,nrow(em))
    z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
  }else{
    z.p<-2*(1-pnorm(abs(z.score)))
  }
  qval<-p.adjust(z.p,method="BH")
  res<-list(n=n,mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,zval=z.score,pval=z.p,FDR=qval)
  return(res)
}

get.REM2<-function(em,vm){
  k<-ncol(em)
  Q.val<-get.Q(em,vm)
  tau2<-get.tau2(Q.val,vm,k)
  temp.wt<-1/(vm+tau2)	
  mu.hat<-rowSums(temp.wt*em,na.rm=TRUE)/rowSums(temp.wt,na.rm=TRUE)
  mu.var<-1/rowSums(temp.wt,na.rm=TRUE)
  Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
  z.score<-mu.hat/sqrt(mu.var)
  z.p<-2*(1-pnorm(abs(z.score)))
  qval<-p.adjust(z.p,method="BH")
  res<-list(mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,zval=z.score,pval=z.p,FDR=qval)
  return(res)
}

get.Q<-function(em,vm){
  wt <- 1/vm
  temp1 <- wt * em
  mu.hat <- rowSums(temp1,na.rm=TRUE)/rowSums(wt,na.rm=TRUE)
  Q <- rowSums(wt * (em - mu.hat)^2,na.rm=TRUE)
  return(Q)
}

get.tau2<-function(Q,vm,k){
  wt<-1/vm
  s1 <- rowSums(wt,na.rm=TRUE)
  s2 <- rowSums(wt^2,na.rm=TRUE)
  temp<- (Q - (k - 1))/(s1 - (s2/s1))
  tau2<-pmax(temp,0)
  return(tau2)
}
