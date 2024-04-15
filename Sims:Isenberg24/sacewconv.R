#estimator and variance of estimator using random intercept 

#required source, must correctly set file path
source("/home/danepi/saceglm.R")
#adding a small sample correction 

sacecluster<-function(data,trt="A",surv="S",out="Y",clustid="Id",indv="X",set1=T,set2=F,conf=.95,boot=F,logform=T,
                      partial=T,nagq=10,iters=200){
  #functions required for analytic method numerical integration
  if(boot==F){
    f1<-function(x){
      fin<-(-x^2/2+sum(dfint$S*x-log(1+exp(dfcov%*%betahat+x))))
      return(fin)
    }
    
    f2<-function(x){
      fin<-f1(x)+log(x^2)
      return(fin)
    }
    
    f3<-function(x){
      fin<-f1(x)+log(x^4)
      return(fin)
    }
    
    #vector of functions
    #indexed by
    f4<-function(x,a){
      g<-f1(x)
      v<-exp(dfcov%*%betahat+x)/(1+exp(dfcov%*%betahat+x))
      ncov<-ncol(dfcov)
      vrep<-matrix(rep(v,ncov),ncol=ncov)
      vfin<-apply(dfcov*vrep,2,sum)
      fin<-g+sign(vfin)*log(abs(vfin))
      if(vfin[a]==0){
        fin[a]<-0
      }
      return(fin[a])
    }
    
    f5<-function(x,a){
      g<-f1(x)
      v<-x^2*exp(dfcov%*%betahat+x)/(1+exp(dfcov%*%betahat+x))
      ncov<-ncol(dfcov)
      vrep<-matrix(rep(v,ncov),ncol=ncov)
      vfin<-apply(dfcov*vrep,2,sum)
      fin<-g+sign(vfin)*log(abs(vfin))
      if(vfin[a]==0){
        fin[a]<-0
      }
      return(fin[a])
    }
    
    #matrix of functions
    #indexed by a,b
    f6<-function(x,a,b){
      g<-f1(x)
      v<-exp(dfcov%*%betahat+x)/((1+exp(dfcov%*%betahat+x))^2)
      ncov<-ncol(dfcov)
      vrep<-matrix(rep(v,ncov),ncol=ncov)
      vprod<-dfcov*vrep
      vfin<-matrix(rep(0,ncov^2),nrow=ncov)
      for(i in 1:nrow(vprod)){
        vfin<-vfin+vprod[i,]%*%t(vprod[i,])
      }
      fin<-g+sign(vfin)*log(abs(vfin))
      if(vfin[a,b]==0){
        fin[a,b]<-0
      }
      return(fin[a,b])
    }
    
    f7<-function(x,a,b){
      g<-f1(x)
      v<-exp(dfcov%*%betahat+x)/(1+exp(dfcov%*%betahat+x))
      ncov<-ncol(dfcov)
      vrep<-matrix(rep(v,ncov),ncol=ncov)
      vprod<-dfcov*vrep
      vsum1<-apply(vprod,2,sum)
      vfin<-vsum1%*%t(vsum1)
      fin<-g+sign(vfin)*log(abs(vfin))
      if(vfin[a,b]==0){
        fin[a,b]<-0
      }
      return(fin[a,b])
    }
    
    
    #adaptive Gaussian-Hermite quadrature for arbitrary function
    #unspecified arguments to allow for indexing
    quadfun<-function(x,w,mode,d2,f,...){
      v<-sapply(mode+sqrt(2)*sqrt(1/abs(d2))*x,FUN=f,...)
      sqrt(2)*sqrt(1/abs(d2))*sum(w*exp(x^2+v))
    }
    
    
    #allows for iterations of above
    #am=0, bm=0, scalar
    #am=k, bm=0, k dim vector, k>1
    #am=k, bm=k, matrix k by k
    aghqvect<-function(f,am=0,bm=0,h=10^(-5)){
      #provides weights and evaluation points
      #nagq is number of quad points determined by wrapper function
      rule<-fastGHQuad::gaussHermiteData(nagq)
      x<-rule$x
      w<-rule$w
      if(am==0 & bm==0){
        #finds max/mode value
        k<-optimize(f,lower=-10,upper=10,maximum = T)$maximum
        #approximate for second derivative
        d2<-(f(k+h)-2*f(k)+f(k-h))/(h^2)
        #approximated value using AGQ
        val<-quadfun(x=x,w=w,mode=k,d2=d2,f=f)
      }
      #repeat above for vectors of integrals
      if(am!=0 & bm==0){
        k<-rep(0,am)
        #out<-rep(0,am)
        d2<-rep(0,am)
        val<-rep(0,am)
        for(i in 1:am){
          #prevent integrals that are exactly 0 from being computed within a cluster
          #occurs when there is a binary cluster-level variable
          #over multiple scenarios .5 seemed like a consistenly good value to check
          if(f(.5,a=i)==0){
            val[i]<-0
            next
          }else{
            k[i]<-optimize(f,lower=-10,upper=10,a=i,maximum = T)$maximum
            #out[i]<--f(k[i],a=i)
            d2[i]<-(f(k[i]+h,a=i)-2*f(k[i],a=i)+f(k[i]-h,a=i))/(h^2)
            #val[i]<-sqrt(2*pi/abs(d2[i]))*exp(out[i])
            val[i]<-quadfun(x=x,w=w,mode=k[i],d2=d2[i],f=f,a=i)}}}
      #repeat above over double index
      if(am!=0 & bm!=0){
        k<-matrix(rep(0,am*bm),nrow=am)
        out<-matrix(rep(0,am*bm),nrow=am)
        d2<-matrix(rep(0,am*bm),nrow=am)
        val<-matrix(rep(0,am*bm),nrow=am)
        for(i in 1:am){
          for(j in 1:bm){
            if(f(.5,a=i,b=j)==0){
              val[i,j]<-0
            }else{
              k[i,j]<-optimize(f,lower=-10,upper=10,a=i,b=j,maximum = T)$maximum
              out[i,j]<--f(k[i,j],a=i,b=j)
              d2[i,j]<-(f(k[i,j]+h,a=i,b=j)-2*f(k[i,j],a=i,b=j)+f(k[i,j]-h,a=i,b=j))/(h^2)
              #val[i,j]<-sqrt(2*pi/abs(d2[i,j]))*exp(out[i,j])
              val[i,j]<-quadfun(x=x,w=w,mode=k[i,j],d2=d2[i,j],f=f,a=i,b=j)}}}
      }
      return(val)
    }}
  
  names<-c(trt,surv,out,clustid,indv)
  #sace point estimates
  saceestim<-function(data,trt,surv,out,clustid,indv,set1,set2){
    #allows for one to include cluster size as a covariate when ics=T
    #potential informative cluster size
    df1<-dplyr::select(data,any_of(names))
    #number of clusters
    nc<-length(unique(df1[,4]))
    #renaming columns for notation
    colnames(df1)[1:4]<-c("A","S","Y","Id")
    labs<-c(colnames(df1[,-c(2,3,4)]),"(1|Id)")
    #glmm logistic model fit with a random intercept
    glmmfit<-lme4::glmer(paste("S~ ", paste(labs, collapse = "+")),family="binomial",data=df1)
    #fixed effects coefficients
    betafix<-lme4::fixef(glmmfit)
    #variance of random intercept
    sigma2<-(as.data.frame(lme4::VarCorr(glmmfit))$sdcor[1])^2
    #adding other treatment to df
    dfnew<-cbind(df1[,-1],A=1-df1$A,Aold=df1$A)
    #estimated predicted probability under other treatment
    valst<-predict(glmmfit,newdata=dfnew,type="response")
    valso<-predict(glmmfit,newdata=df1,type="response")
    dfnewv<-cbind(dfnew,pred=valst,predog=valso)
    #set1 estimator
    if(set1==T){
      exp1h<-sum(dfnewv$pred*dfnewv$Y*dfnewv$S*dfnewv$Aold)/sum(dfnewv$pred*dfnewv$S*dfnewv$Aold)
      exp0h<-sum(dfnewv$pred*dfnewv$Y*dfnewv$S*(1-dfnewv$Aold))/sum(dfnewv$pred*dfnewv$S*(1-dfnewv$Aold))
      set1est<-exp1h-exp0h}
    
    #set2 estimator
    if(set2==T){
      #prediction under treatment received
      exp1j<-sum(dfnewv$pred*dfnewv$Y*dfnewv$S*dfnewv$Aold/(dfnewv$predog))/sum(dfnewv$pred*dfnewv$Aold*dfnewv$S/(dfnewv$predog))
      exp0j<-sum(dfnewv$Y*dfnewv$S*(1-dfnewv$Aold))/sum(dfnewv$S*(1-dfnewv$Aold))
      set2est<-exp1j-exp0j}
    
    #alternatives for outputs given information
    #could futher subset by what's necessary for analytic vs boot, but slows down
    if(set1==T & set2==T){
      return(list(df1,nc,valst,valso,betafix,sigma2,c(exp1h,exp0h,exp1j,exp0j),set1est,set2est))
    }
    if(set1==T & set2==F){
      return(list(df1,nc,valst,valso,betafix,sigma2,c(exp1h,exp0h),set1est))
    }
    if(set1==F & set2==T){
      return(list(df1,nc,valst,valso,betafix,sigma2,c(exp1j,exp0j),set2est))
    }
    if(set1==F & set2==F){
      stop("You must specify an estimator")
    }
  }
  
  
  #output of function
  results<-saceestim(data,trt,surv,out,clustid,indv,set1,set2)
  
  df<-results[[1]]
  nc<-results[[2]]
  valst<-results[[3]]
  valso<-results[[4]]
  betahat<-results[[5]]
  sigma2<-results[[6]]
  
  
  #for analytic variances, needs for matrices
  if(set1==T & set2==T){
    exp1h<-results[[7]][[1]]
    exp0h<-results[[7]][[2]]
    exp1j<-results[[7]][[3]]
    exp0j<-results[[7]][[4]]
  }
  if(set1==T & set2==F){
    exp1h<-results[[7]][[1]]
    exp0h<-results[[7]][[2]]
  }
  if(set1==F & set2==T){
    exp1j<-results[[7]][[1]]
    exp0j<-results[[7]][[2]]
  }
  
  #obtaining estimates
  estimators<-unlist(results[-c(1:7)])
  
  newnames<-c("int","A",indv)
  #analytic variance
  if(boot==F){
    #add to real function but give a different warning 
    #boundary condition 
    if(sigma2<5*10^(-4)){
      #add penalty
      final<-saceglm(data,trt,surv,out,clustid,indv,crobust=T,varpen=T,set1,set2,conf,boot,iters) 
      return(c(final))
      stop()
    }else{
      #number of covariates, adding 1 for intercept
      ncov<-ncol(df[,-c(2,3,4)])+1
      #empty m-estimator inner matrix, h=set1, j=set2
      phimath<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
      phimatj<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
      #empty m-estimator outer matrix, before inversion
      Ah<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
      Aj<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
      #add intercept do data frame
      dffullint<-data.frame(int=rep(1,nrow(df)),df)
      #for each cluster compute that terms as described in accompanying paper
      for(i in 1:nc){
        clust1<-which(df$Id==i,)
        dfclust1<-dffullint[clust1,]
        dfint<-dfclust1
        #covariate only matrices under treatment received
        dfcov<-dplyr::select(dfclust1,any_of(newnames))
        dfcovop<-dfcov
        #under alternative treatment
        dfcovop$A<-1-dfcovop$A
        dfcov<-as.matrix(dfcov)
        dfcovop<-as.matrix(dfcovop)
        #estimated values under each treatment received vs alternative
        p0hat<-valst[clust1]
        p1hat<-valso[clust1]
        nvar<-ncol(dfcov)
        
        if(logform==F){
          #entries of meat matrix common to both estimators
          phi1<-apply(dfcov*dfint$S,2,sum)-aghqvect(f4,am=nvar)/aghqvect(f1)
          phi2<--1/(2*sigma2)+1/(2*sigma2^2)*aghqvect(f2)/aghqvect(f1)
          #phi2<--1/(2*.21^2)+1/(2*.21^4)*aghqvect(f2)/aghqvect(f1)
          
          #entries of outer matrix common to both estimators before inversion
          A11<--(aghqvect(f6,am=nvar,bm=nvar)+aghqvect(f7,am=nvar,bm=nvar))/aghqvect(f1)+aghqvect(f4,am=nvar)%*%t(aghqvect(f4,am=nvar))/(aghqvect(f1))^2
          A12<--1/(2*sigma2^2)*(aghqvect(f5,am=nvar)/aghqvect(f1)-aghqvect(f4,am=nvar)/(aghqvect(f1))^2*aghqvect(f2))
          A22<-1/(2*sigma2^2)-1/(sigma2^3)*aghqvect(f2)/aghqvect(f1)+1/(4*sigma2^4)*(aghqvect(f3)/aghqvect(f1)-(aghqvect(f2))^2/(aghqvect(f1))^2)}
        
        #applies a partial log form the above when computation is unfeasible, recommended for cluster sizes >100
        if(logform==T & partial==T){
          phi1<-apply(dfcov*dfint$S,2,sum)-aghqvect(f4,am=nvar)/aghqvect(f1)
          phi2<--1/(2*sigma2)+1/(2*sigma2^2)*aghqvect(f2)/aghqvect(f1)
          A11<--(aghqvect(f6,am=nvar,bm=nvar)+aghqvect(f7,am=nvar,bm=nvar))/aghqvect(f1)+exp(log(aghqvect(f4,am=nvar)%*%t(aghqvect(f4,am=nvar)))-2*log(abs(aghqvect(f1))))
          A12<--1/(2*sigma2^2)*(aghqvect(f5,am=nvar)/aghqvect(f1)-
                                  exp(log(aghqvect(f4,am=nvar))-2*log(abs(aghqvect(f1)))+sign(aghqvect(f2))*log(abs(aghqvect(f2)))))
          
          A22<-1/(2*sigma2^2)-1/(sigma2^3)*aghqvect(f2)/aghqvect(f1)+1/(4*sigma2^4)*(aghqvect(f3)/aghqvect(f1)-exp(2*log(abs(aghqvect(f2)))-2*log(abs(aghqvect(f1)))))
        }
        
        #full log form of above
        if(logform==T & partial==F){
          phi1<-apply(dfcov*dfint$S,2,sum)-ifelse(aghqvect(f4,am=nvar)!=0,exp(sign(aghqvect(f4,am=nvar))*log(abs(aghqvect(f4,am=nvar)))-sign(aghqvect(f1))*log(abs(aghqvect(f1)))),0)
          phi2<--1/(2*sigma2)+1/(2*sigma2^2)*
            ifelse(aghqvect(f2)!=0,exp(sign(aghqvect(f2))*log(abs(aghqvect(f2)))-sign(aghqvect(f1))*log(abs(aghqvect(f1)))),0)
          A11<-ifelse((aghqvect(f6,am=nvar,bm=nvar)+aghqvect(f7,am=nvar,bm=nvar))!=0,-exp(sign(aghqvect(f6,am=nvar,bm=nvar)+aghqvect(f7,am=nvar,bm=nvar))*log(abs(aghqvect(f6,am=nvar,bm=nvar)+aghqvect(f7,am=nvar,bm=nvar)))-sign(aghqvect(f1))*log(abs(aghqvect(f1)))),0)+ifelse(aghqvect(f4,am=nvar)%*%t(aghqvect(f4,am=nvar))!=0,exp(log(aghqvect(f4,am=nvar)%*%t(aghqvect(f4,am=nvar)))-2*log(abs(aghqvect(f1)))),0)
          A12<--1/(2*sigma2^2)*(ifelse(aghqvect(f5,am=nvar)!=0,exp(sign(aghqvect(f5,am=nvar))*log(abs(aghqvect(f5,am=nvar)))-sign(aghqvect(f1))*log(abs(aghqvect(f1)))),0)-ifelse(aghqvect(f4,am=nvar)!=0,exp(log(aghqvect(f4,am=nvar))-2*log(abs(aghqvect(f1)))+sign(aghqvect(f2))*log(abs(aghqvect(f2)))),0))
          A22<-1/(2*sigma2^2)-1/(sigma2^3)*ifelse(aghqvect(f2)!=0,exp(sign(aghqvect(f2))*log(abs(aghqvect(f2)))-sign(aghqvect(f1))*log(abs(aghqvect(f1)))),0)+1/(4*sigma2^4)*(ifelse(aghqvect(f3)!=0,exp(sign(aghqvect(f3))*log(aghqvect(f3))-sign(aghqvect(f1))*log(aghqvect(f1))),0)-ifelse(aghqvect(f2)!=0,exp(2*log(abs(aghqvect(f2)))-2*log(abs(aghqvect(f1)))),0))}
        
        #set1 only entries to matrices
        if(set1==T){
          #inner
          phi3h<-sum(dfint$Y*dfint$A*dfint$S*p0hat)-exp1h*sum(dfint$A*dfint$S*p0hat)
          phi4h<-sum(dfint$Y*(1-dfint$A)*dfint$S*p1hat)-exp0h*sum((1-dfint$A)*dfint$S*p1hat)
          #outer before inversion
          A31h<-apply(dfcovop*dfint$Y*dfint$A*dfint$S*p0hat*(1-p0hat),2,sum)-exp1h*apply(dfcovop*dfint$A*dfint$S*p0hat*(1-p0hat),2,sum)
          A41h<-apply(dfcovop*dfint$Y*(1-dfint$A)*dfint$S*p1hat*(1-p1hat),2,sum)-exp0h*apply(dfcovop*(1-dfint$A)*dfint$S*p1hat*(1-p1hat),2,sum)
          A33h<-sum(dfint$A*dfint$S*p0hat)
          A44h<-sum((1-dfint$A)*dfint$S*p1hat)
          
          #filling the matrix as one cycles through cluster indices
          phimatih<-c(phi1,phi2,phi3h,phi4h)%*%t(c(phi1,phi2,phi3h,phi4h))
          phimath<-phimath+phimatih
          
          Aih<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
          
          Aih[1:ncov,1:ncov]<-A11
          Aih[1:ncov,ncov+1]<-A12
          Aih[ncov+1,1:ncov]<-A12
          Aih[ncov+1,ncov+1]<-A22
          Aih[ncov+2,1:ncov]<-A31h
          Aih[ncov+3,1:ncov]<-A41h
          Aih[ncov+2,ncov+2]<-A33h
          Aih[ncov+3,ncov+3]<-A44h
          
          
          Ah<-Aih+Ah
        }
        
        #set2 only entries to matrices
        if(set2==T){
          #inner
          phi3j<-sum(dfint$Y*dfint$A*dfint$S*p0hat/p1hat)-exp1j*sum(dfint$A*dfint$S*p0hat/p1hat)
          phi4j<-sum(dfint$Y*(1-dfint$A)*dfint$S)-exp0j*sum((1-dfint$A)*dfint$S)
          #outer
          A31j<-apply(dfcovop*dfint$Y*dfint$A*dfint$S*p0hat*(1-p0hat)/p1hat-dfcov*dfint$Y*dfint$A*dfint$S*p0hat*(1-p1hat)/p1hat,2,sum)-exp1j*apply(dfcovop*dfint$A*dfint$S*p0hat*(1-p0hat)/p1hat-dfcov*dfint$A*dfint$S*p0hat*(1-p1hat)/p1hat,2,sum)
          A41j<-rep(0,ncov)
          A33j<-sum(dfint$A*dfint$S*p0hat/p1hat)
          A44j<-sum((1-dfint$A)*dfint$S)
          
          #filling the matrix as one cycles through cluster indices
          phimatij<-c(phi1,phi2,phi3j,phi4j)%*%t(c(phi1,phi2,phi3j,phi4j))
          phimatj<-phimatj+phimatij
          
          Aij<-matrix(rep(0,(ncov+3)*(ncov+3)),nrow=ncov+3)
          
          Aij[1:ncov,1:ncov]<-A11
          Aij[1:ncov,ncov+1]<-A12
          Aij[ncov+1,1:ncov]<-A12
          Aij[ncov+1,ncov+1]<-A22
          Aij[ncov+2,1:ncov]<-A31j
          Aij[ncov+3,1:ncov]<-A41j
          Aij[ncov+2,ncov+2]<-A33j
          Aij[ncov+3,ncov+3]<-A44j
          
          Aj<-Aij+Aj}
      }
    }
    
    #determining which category of output, decided by what type of estimator chosen
    alpha2<-(1-conf)/2
    zl<-qnorm(alpha2)
    zu<--zl
    if(set1==T & set2==F){
      Ainh<-solve(Ah)
      sam<-Ainh%*%phimath%*%t(Ainh)
      v1max<-max(abs(sam))
      varsest<-nc/(nc-nrow(Ah))*t(c(rep(0,ncov+1),1,-1))%*%Ainh%*%phimath%*%t(Ainh)%*%c(rep(0,ncov+1),1,-1)
      lb<-estimators+zl*sqrt(as.numeric(varsest))
      ub<-estimators+zu*sqrt(as.numeric(varsest))
      bounds<-c(lb,ub)
      convg<-ifelse(v1max>=nc,0,1)
      finnames<-c("EstimateSet1","VarEstSet1","LBS1","UBS1","Conv")}
    
    if(set1==F & set2==T){
      Ainj<-solve(Aj)
      sam<-Ainj%*%phimatj%*%t(Ainj)
      v1max<-max(abs(sam))
      varsest<-nc/(nc-nrow(Aj))*t(c(rep(0,ncov+1),1,-1))%*%Ainj%*%phimatj%*%t(Ainj)%*%c(rep(0,ncov+1),1,-1)
      lb<-estimators+zl*sqrt(as.numeric(varsest))
      ub<-estimators+zu*sqrt(as.numeric(varsest))
      bounds<-c(lb,ub)
      convg<-ifelse(v1max>=nc,0,1)
      finnames<-c("EstimateSet2","VarEstSet2","LBS2","UBS2","Conv")}
    
    if(set1==T & set2==T){
      Ainh<-solve(Ah)
      #varsesth<-t(c(rep(0,ncov+1),1,-1))%*%Ainh%*%phimath%*%t(Ainh)%*%c(rep(0,ncov+1),1,-1)
      samh<-Ainh%*%phimath%*%t(Ainh)
      #v11h<-samh[1,1]
      v1maxh<-max(abs(samh))
      varsesth<-nc/(nc-nrow(Ah))*t(c(rep(0,ncov+1),1,-1))%*%Ainh%*%phimath%*%t(Ainh)%*%c(rep(0,ncov+1),1,-1)
      lb1<-estimators[1]+zl*sqrt(as.numeric(varsesth))
      ub1<-estimators[1]+zu*sqrt(as.numeric(varsesth))
      Ainj<-solve(Aj)
      #don't need both
      #samj<-Ainj%*%phimatj%*%t(Ainj)
      #v11j<-samj[1,1]
      #v1maxj<-max(abs(samj))
      #varsestj<-t(c(rep(0,ncov+1),1,-1))%*%Ainj%*%phimatj%*%t(Ainj)%*%c(rep(0,ncov+1),1,-1)
      varsestj<-nc/(nc-nrow(Aj))*t(c(rep(0,ncov+1),1,-1))%*%Ainj%*%phimatj%*%t(Ainj)%*%c(rep(0,ncov+1),1,-1)
      lb2<-estimators[2]+zl*sqrt(as.numeric(varsestj))
      ub2<-estimators[2]+zu*sqrt(as.numeric(varsestj))
      varsest<-c(varsesth,varsestj)
      bounds<-c(lb1,ub1,lb2,ub2)
      #convg<-v1maxh
      convg<-ifelse(v1maxh>=nc,0,1)
      #finnames<-c("EstimateSet1","EstimateSet2","VarEstSet1","VarEstSet2","LBS1","UBS1","LBS2","UBS2","Conv11h",
                  #"ConvMaxh","Conv11j","ConvMaxj")
      finnames<-c("EstimateSet1","EstimateSet2","VarEstSet1","VarEstSet2","LBS1","UBS1","LBS2","UBS2","Conv")}
  }
  
  #non parametric bootstrap
  if(boot==T){
    #number of types of estimators specified
    colboot<-sum(c(set1,set2))
    #empty boot matrix by number of bootstrap iterations
    bootsample<-matrix(rep(0,iters*colboot),ncol=colboot)
    #non-parametric bootstrap, sampling with replacement entire clusters
    for(i in 1:iters){
      bootind<-sample(1:nc,replace=T)
      dflist<-list()
      for(j in 1:nc){
        dflist[[j]]<-df[df$Id==bootind[j],]
      }
      #boostrap data with potentially repeated clusters
      dfboot<-do.call(rbind,dflist)
      names<-names(df)
      #reestimate parameters on boot data
      resultsb<-suppressMessages(saceestim(data=dfboot,trt,surv,out,clustid,indv,set1,set2))
      #generate bootstrap estimates
      bootsample[i,]<-unlist(resultsb[-c(1:7)])
    }
    if(set1==T & set2==F){
      finnames<-c("EstimateSet1","VarEstSet1","LBS1","UBS1")}
    
    if(set1==F & set2==T){
      finnames<-c("EstimateSet2","VarEstSet2","LBS2","UBS2")}
    
    if(set1==T & set2==T){
      finnames<-c("EstimateSet1","EstimateSet2","VarEstSet1","VarEstSet2","LBS1","UBS1","LBS2","UBS2")}
    
    #find variance of non-parametric bootstrap distribution
    varsest<-apply(bootsample,2,var)
    alpha2<-(1-conf)/2
    percl<-alpha2
    percu<-1-alpha2
    bounds<-apply(bootsample, 2, quantile, probs = c(percl,percu))}
  
  if(boot==F){
  final<-c(estimators,varsest,bounds,convg)
  names(final)<-finnames
  return(c(final[1:(length(finnames)-1)],Degen=0,final[length(finnames)]))}
  if(boot==T){
  final<-c(estimators,varsest,bounds) 
  names(final)<-finnames
  return(c(final[1:(length(finnames))],Degen=0))
  }
}