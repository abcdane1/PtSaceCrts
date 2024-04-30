#simulations for bootstrap performance and time, parallelized for mac or linux 
#analytic variance only, boot set to True 
#creates values as in Table S4

#names of sourced files, for running simulations  
#must set your own file paths
source("/home/danepi/sacewconv.R")
source("/home/danepi/saceglm.R")

#simulation output names
names<-c("TruthAvg","Set1EstAvg","Set2EstAvg", "Set1BiasEst", "Set2BiasEst",
         "Set1MCVarEst","Set2MCVarEst","Set1VarAvg","Set2VarAvg","Set1Cov","Set2Cov","Sigma2Degen",
         "Always-Survivors","Protected","Never-Survivors", "Harmed","EmpCorrSurv","Time")

#wrapper for simulation
newwrapparamest_sim<-function(nc,a,b,alphamon,alphasurv,iccy,iccs,crobust,boot,iters,nsim,ncores){
  RNGkind("L'Ecuyer-CMRG")
  set.seed(01051990)
  newparamest_sim<-function(nc,a,b,alphamon,alphasurv,iccy,iccs,crobust,boot,iters,nsim){
    newdata_sim<-function(nc,a,b,alphamon,alphasurv,iccy,iccs,nsim){
      #cluster sizes
      cs<-sample(a:b,size=nc,replace=T)
      #number of individuals
      n<-sum(cs)
      
      #cluster id
      id<-rep(1:nc,cs)
      
      #intercept 
      #alphas<-0 old
      
      #new
      alphas<-.75
      
      #alphay<-0.5 old
      #new
      alphay<-1
      
      #individual variables
      x1<-rep(0,n)
      x2<-rep(0,n)
      
      #effect on survival
      
      #old
      #beta1s<-0.5
      #beta2s<--0.25
      
      #new 
      beta1s<-.1
      beta2s<--0.05
      
      
      #effect on outcome
      #old
      #beta1y<-0.1
      #beta2y<-0.25
      
      #new
      beta1y<-0.25
      beta2y<-0.125
      
      
      #cluster variable
      c1<-rep(0,n)
      
      #effect on survival 
      
      #old
      #beta3s<-0.5
      
      #new
      beta3s<-0.1
      
      
      #effect on outcome
      beta3y<-0
      
      #survival vectors
      s1<-rep(0,n) #treatment 1
      s0<-rep(0,n) #treatment 0
      s<-rep(0,n)  #observed treatment
      
      tv<-rep(0,n) #observed treatment of group
      trt<-rep(0,nc) #treatment for individuals
      
      #effect on survival (i.e. monotonicity parameter) #alphamon
      #effect on outcome 
      betaty<-1
      
      #nonmortal outcome vectors
      y1<-rep(0,n)
      y0<-rep(0,n)
      y<-rep(0,n)
      
      #initiate for loop
      c<-1
      
      #for random intercept
      #ICC
      #residual variance 
      #varsg<-pi^2/3
      #sig2s<-(varsg*rho)/(1-rho)
      
      vary<-1
      sig2y<-(vary*iccy)/(1-iccy)
      
      #scale such that iccs is as stated 
      gam<-pi/sqrt(3)*sqrt((1-iccy)/(vary*iccy))*sqrt(iccs/(1-iccs))
      
      #random intercept 
      #bi<-rnorm(nc,mean=0,sd=sqrt(sig2s))
      biy<-rnorm(nc,mean=0,sd=sqrt(sig2y))
      bi<-gam*biy #what happens if drawn independently
      #bi<-gam*rnorm(nc,mean=0,sd=sqrt(sig2y))
      for(i in 1:nc){
        #treatment assignment
        t<-rbinom(1,1,.5)
        #treatment by group
        trt[i]<-t
        #treatment by individual
        tv[c:(c+cs[i]-1)]<-rep(t,cs[i])
        
        #population info
        
        #individual covariates
        mean1<-rep(2,cs[i])
        sigma21<-.5*diag(cs[i])
        mean2<-rep(.5,cs[i])
        sigma22<-.25*diag(cs[i])
        val1<-MASS::mvrnorm(n=1,mu=mean1,Sigma=sigma21)
        val2<-MASS::mvrnorm(n=1,mu=mean2,Sigma=sigma22)
        
        #cluster covariates
        val3<-rep(rbinom(1,1,.3),cs[i])
        
        #filling in vectors 
        x1[c:(c+cs[i]-1)]<-val1
        x2[c:(c+cs[i]-1)]<-val2
        c1[c:(c+cs[i]-1)]<-val3
        
        #pop survival under treatment 1
        #alphamon, impact of treatment on survival
        sp1<-1/(1+exp(-alphas-alphamon-beta1s*val1-beta2s*val2-beta3s*val3-bi[i]))
        s1val<-rbinom(n=cs[i],size=1,p=sp1)
        s1[c:(c+cs[i]-1)]<-s1val
        #pop survival under treatment 0
        #alphasurv, effect of survival under trt=1 on survival under trt=0
        sp0<-1/(1+exp(-alphas-alphasurv*s1val-beta1s*val1-beta2s*val2-beta3s*val3-bi[i]))
        s0val<-rbinom(n=cs[i],size=1,p=sp0)
        s0[c:(c+cs[i]-1)]<-s0val
        
        #error terms for pop non-mortal outcomes 
        eiy1<-rnorm(cs[i])
        eiy0<-rnorm(cs[i])
        
        #uder treatment 1
        y1val<-(betaty+1)*(alphay+beta1y*val1+beta2y*val2+beta3y*val3)+biy[i]+eiy1
        y1[c:(c+cs[i]-1)]<-y1val
        #under treatment 0
        y0val<-alphay+beta1y*val1+beta2y*val2+beta3y*val3+biy[i]+eiy0
        y0[c:(c+cs[i]-1)]<-y0val
        
        #observed values 
        if(t==1){
          s[c:(c+cs[i]-1)]<-s1val
          y[c:(c+cs[i]-1)]<-y1val
        }else{
          s[c:(c+cs[i]-1)]<-s0val
          y[c:(c+cs[i]-1)]<-y0val
        }
        c<-c+cs[i]
      }
      #observed df
      dfsim<-data.frame(Id=id,Y=y,S=s,X1=x1,X2=x2,C1=c1,A=tv)
      #pop df 
      dftruth<-data.frame(Id=id,Y1=y1,Y0=y0,S1=s1,S0=s0,X1=x1,X2=x2,C1=c1,A=tv)
      return(list(dfsim,dftruth))
    }
    
    dfsim<-newdata_sim(nc,a,b,alphamon,alphasurv,iccy,iccs,nsim)[[1]]
    #fitting and timing for SACE GLMM boot 
    start_timeglmm <- Sys.time()
    out<-sacecluster(data=dfsim,indv=c("X1","X2","C1"),set1=T,set2=T,conf=.95,boot=boot,iters=iters)
    end_timeglmm <- Sys.time()
    
    #degen replace 
    #degen<-function(df,y,k=0){
    #if(y[10]==0 & is.na(y[10])==F){
    #dfsim<-newdata_sim(nc,a,b,alphamon,alphasurv,iccy,iccs,nsim)[[1]]
    #out<-sacecluster(data=dfsim,indv=c("X1","X2","C1"),set1=T,set2=T,conf=.95,boot=boot,iters=iters)
    #k<-k+1
    #degen(df=dfsim,y=out,k)
    #  }else{
    #    return(list(y,k))}}
    
    #convergence check, ignored since boot=F
    if(boot==F){
    degen<-function(df,y){
      if(y[10]==0 & is.na(y[10])==F){
        out<-rep(NA,length(y))
        return(list(out,1))
      }else{
        return(list(y,0))}}
    
    ldegen<-degen(df=dfsim,y=out)
    out<-ldegen[[1]]
    convgind<-ldegen[[2]]}
    
    #GLMM output 
    set1est<-out[1]
    set2est<-out[2]
    set1var<-out[3]
    set2var<-out[4]
    lbs1<-out[5]
    ubs1<-out[6]
    lbs2<-out[7]
    ubs2<-out[8]
    degen<-out[9]
    #GLMM time 
    timewith<-end_timeglmm-start_timeglmm
    
    #fitting and timing of SACE GLM boot 
    start_timeglm <- Sys.time()
    outwo<-saceglm(data=dfsim,indv=c("X1","X2","C1"),crobust=crobust,varpen=F,set1=T,set2=T,conf=.95,boot=boot,iters=iters)
    end_timeglm <- Sys.time()
    
    #GLM output 
    set1estwo<-outwo[1]
    set2estwo<-outwo[2]
    set1varwo<-outwo[3]
    set2varwo<-outwo[4]
    lbs1wo<-outwo[5]
    ubs1wo<-outwo[6]
    lbs2wo<-outwo[7]
    ubs2wo<-outwo[8]
    degenwo<-outwo[9]
    #GLM time
    timewo<-end_timeglm-start_timeglm
    
    #make more clusters
    dftrue<-newdata_sim(nc=1000,a,b,alphamon,alphasurv,iccy,iccs,nsim)[[2]]
    y1t<-dftrue$Y1
    y0t<-dftrue$Y0
    s1t<-dftrue$S1
    s0t<-dftrue$S0
    truth<-(sum(s1t*s0t*y1t)-sum(s1t*s0t*y0t))/sum(s1t*s0t)
    
    #coverage
    covset1<-ifelse(truth<ubs1 & truth>lbs1,1,0)
    covset2<-ifelse(truth<ubs2 & truth>lbs2,1,0)
    covset1wo<-ifelse(truth<ubs1wo & truth>lbs1wo,1,0)
    covset2wo<-ifelse(truth<ubs2wo & truth>lbs2wo,1,0)
    
    #observed stratum 
    s1<-dftrue$S1
    s0<-dftrue$S0
    np<-nrow(dftrue)
    alwayssurv<-length(which(s1==1 & s0==1))/np
    protected<-length(which(s1==1 & s0==0))/np
    neversurv<-length(which(s1==0 & s0==0))/np
    #montonicity violations
    harmed<-1-alwayssurv-protected-neversurv
    #marginal cor
    pcor<-cor(s1,s0)
    
    with<-c(truth,set1est,set2est,set1est-truth,set2est-truth,set1var,set2var,covset1,covset2,degen,alwayssurv,protected,neversurv,harmed,pcor,timewith)
    without<-c(truth,set1estwo,set2estwo,set1estwo-truth,set2estwo-truth,set1varwo,set2varwo,
               covset1wo,covset2wo,degenwo,alwayssurv,protected,neversurv,harmed,pcor,timewo)
    c(with,without)
    
  }
  
  l<-parallel::mcmapply(newparamest_sim,nsim=1:nsim,MoreArgs=list(nc=nc,a=a,b=b,alphamon=alphamon,alphasurv=alphasurv,iccy
                                                                  =iccy,iccs=iccs,crobust=crobust,
                                                                  boot=T,iters=250),mc.cores=ncores)
  lw<-l[1:16,] #performance GLM
  lwo<-l[17:32,] #performance GLMM
  
  #averaging results for SACE GLMM, includes time 
  resultswomcw<-apply(lw,1,mean,na.rm=T)
  #mc est of variance 
  mcestvarset1w<-1/((nsim-1))*sum((lw[2,]-resultswomcw[2])^2,na.rm=T)
  mcestvarset2w<-1/((nsim-1))*sum((lw[3,]-resultswomcw[3])^2,na.rm=T)
  resultsw<-c(resultswomcw[1:5],mcestvarset1w,mcestvarset2w,resultswomcw[6:16])
  
  #repeat for SACE GLM
  resultswomcwo<-apply(lwo,1,mean,na.rm=T)

  mcestvarset1wo<-1/((nsim-1))*sum((lwo[2,]-resultswomcwo[2])^2,na.rm=T)
  mcestvarset2wo<-1/((nsim-1))*sum((lwo[3,]-resultswomcwo[3])^2,na.rm=T)
  resultswo<-c(resultswomcwo[1:5],mcestvarset1wo,mcestvarset2wo,resultswomcwo[6:16])
  return(rbind(resultsw,resultswo))
}

#parameter scenarios 
#random intercept
ri<-c(T,F)
#number of clusters 
#ncscen<-c(30,60,90)
ncscen<-c(60)
#large 
#ncscen<-100
#ascen<-c(10,50)
ascen<-c(25) #min csize
#bscen<-c(40,150)
bscen<-c(50) #max csize 
#alphamonscen<-c(log(1.25),0,-log(1.25),log(5))
alphamonscen<-c(0,log(1.25),log(5)) #delta 
#alphamonscen<-log(1.25)
alphasurvscen<-0
iccyscen<-.1
#iccscen cant be identically 0
#iccsscen<-c(.3,.15,.1, 0.05)
iccsscen<-c(.1,.3)
#iccsscen<-c(.05)
#add additional later for near 0

#nscens<-length(ncscen)*nrow(csizescen)*length(alphamonscen)*length(alphasurvscen)*
#length(iccyscen)*length(iccsscen)*length(ri)

#constructing scenario matrix 
scenwosize<-expand.grid(Rint=ri,Nc=ncscen,Delta=alphamonscen,ICCy=iccyscen,ICCs=iccsscen)
scenwosize<-cbind(Rint=scenwosize[,1],scenwosize[,-1])
scenariomat<-do.call("rbind", replicate(length(ascen), scenwosize, simplify = FALSE))
acol<-rep(ascen,each=nrow(scenariomat)/2)
bcol<-rep(bscen,each=nrow(scenariomat)/2)

finscenariomat<-cbind(scenariomat[,1:2],a=acol,b=bcol,scenariomat[,3:ncol(scenariomat)])
simoutput<-matrix(rep(0,nrow(finscenariomat)*length(names)),ncol=length(names))

#for loop to produce output
for(i in 1:(nrow(finscenariomat)/2)){
  simoutput[(2*i-1):(2*i),]<-newwrapparamest_sim(nc=finscenariomat[(2*i-1),2],a=finscenariomat[(2*i-1),3],b=finscenariomat[(2*i-1),4],
                                                 alphamon=finscenariomat[(2*i-1),5],alphasurv=alphasurvscen,
                                                 iccy=finscenariomat[(2*i-1),6],iccs=finscenariomat[(2*i-1),7],crobust=T,
                                                 boot=T,iters=250,nsim=1000,ncores=40)
  
  if(i %% 6 == 0){ 
    filename<-paste("new5sacesceboot",i,".csv",sep="")
    colnames(simoutput)<-names
    finalvals<-cbind(finscenariomat,simoutput)
    write.csv(finalvals,filename)
  }
}
