#' Estimation of SACE with Sandwich Variance
#'
#' This function provides point estimates of SACE based on estimators derived from two
#' sets of identification assumptions as described in (cite hayden, ding/lu, jiang, own, papers).
#' All survival models are fit using GLM for binary data.
#' The function also provides estimates of the variance for these estimators. User has the option
#' to choose how variances are estimated and corresponding confidence intervals
#' are constructed. By default, the function provides the variance of the asymptotic distribution of
#' these estimators (cite wy and stefanski) that is cluster-robust in accordance with a GEE model with
#' a working correlation matrix for independence.
#'
#' @param data Data frame containing all data required for estimation with user-specified names below.
#' @param trt A named `character` specifying treatment variable. Default is "A".
#' @param surv A named `character` specifying survival status, where survival through study is indicated by 1 and death by 0. Default is "S".
#' @param out A named `character` specifying non-mortal outcome. Default is "Y".
#' @param clustid A named `character` specifying non-mortal cluster membership. Default is "Id".
#' @param indv A named `character` vector for covariates to be treated as fixed effects. Group-level variables can be included but they must be defined
#' for each individual in the group. Default is "X".
#' @param crobust A `logical` specifying whether cluster-robust variance estimate is provided (submatrix equal to `vcov(gee`). If set to `F`, the clustering structure is completely ignored
#' in the sandwich variance expression (submatrix equal to `sandwich(glm())`). Default is `T`.
#' @param set1 A `logical` argument for whether identified estimator uses Set 1 Assumptions. Default is `T`.
#' @param set2 A `logical` argument for whether identified estimator uses Set 2 Assumptions. Default is `F`.
#' If `set2=T` and `set2=T`, function will provide results for both estimators.
#' @param conf A `numeric` argument in the interval (0,1) for % confidence interval. Default is `.95`.
#' @param boot A `logical` argument for variance estimation. If `boot=F`, variance is estimated
#' @param iters A `double` for number of bootstrap samples to be taken when `boot=T`. Default is `iters=200`. This argument is ignored when `boot=F`.
#'
#' @return A named `double` including point estimates, estimates of variance, and confidence intervals.
#'
#'
#' @export


#sace glm
saceglm<-function(data,trt="A",surv="S",out="Y",clustid="Id",indv="X",crobust=T,set1=T,set2=T,conf=.95,boot=F,iters=200){
  names<-c(trt,surv,out,clustid,indv)
  #sace estimators
  saceestimglm<-function(data,trt,surv,out,clustid,indv,set1,set2){
    #allows for one to include cluster size as a covariate when ics=T
    #potential informative cluster size
    df1<-dplyr::select(data,any_of(names))
    #number of clusters
    nc<-length(unique(df1[,4]))
    #renaming columns for notation
    colnames(df1)[1:4]<-c("A","S","Y","Id")
    df1$Id<-as.factor(df1$Id)
    labs<-c(colnames(df1[,-c(2,3,4)]))
    #glmm logistic model fit with a random intercept
    glmfit<-glm(paste("S~ ", paste(labs, collapse = "+")),family="binomial",data=df1)
    #fixed effects coefficients
    betahat<-coef(glmfit)
    #variance of random intercept
    #adding other treatment to df
    dfnew<-cbind(df1[,-1],A=1-df1$A,Aold=df1$A)
    #estimated predicted probability under other treatment
    valst<-predict(glmfit,newdata=dfnew,type="response")
    valso<-predict(glmfit,newdata=df1,type="response")
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
      return(list(df1,nc,valst,valso,betahat,c(exp1h,exp0h,exp1j,exp0j),set1est,set2est))
    }
    if(set1==T & set2==F){
      return(list(df1,nc,valst,valso,betahat,c(exp1h,exp0h),set1est))
    }
    if(set1==F & set2==T){
      return(list(df1,nc,valst,valso,betahat,c(exp1j,exp0j),set2est))
    }
    if(set1==F & set2==F){
      stop("You must specify an estimator")
    }
  }


  #output of function
  results<-saceestimglm(data,trt,surv,out,clustid,indv,set1,set2)

  df<-results[[1]]
  nc<-results[[2]]
  valstno<-results[[3]]
  valsono<-results[[4]]
  betahat<-results[[5]]

  #for analytic variances, needs for matrices
  if(set1==T & set2==T){
    exp1h<-results[[6]][[1]]
    exp0h<-results[[6]][[2]]
    exp1j<-results[[6]][[3]]
    exp0j<-results[[6]][[4]]
  }
  if(set1==T & set2==F){
    exp1h<-results[[6]][[1]]
    exp0h<-results[[6]][[2]]
  }
  if(set1==F & set2==T){
    exp1j<-results[[6]][[1]]
    exp0j<-results[[6]][[2]]
  }


  #obtaining estimates
  estimators<-unlist(results[-c(1:6)])

  newnames<-c("int","A",indv)
  if(boot==F){
    if(crobust==F){
      #number of covariates, adding 1 for intercept
      ncov<-ncol(df[,-c(2,3,4)])+1
      #empty m-estimator inner matrix, h=set1, j=set2
      phimatnoh<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      phimatnoj<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      #empty m-estimator outer matrix, before inversion
      Anoh<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      Anoj<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2) #here
      #add intercept do data frame
      dffullint<-data.frame(int=rep(1,nrow(df)),df)

      dfcovfull<-dplyr::select(dffullint,newnames)
      dfcovfullop<-dfcovfull
      #under alternative treatment
      dfcovfullop$A<-1-dfcovfullop$A
      dfcovfull<-as.matrix(dfcovfull)
      dfcovfullop<-as.matrix(dfcovfullop)

      for(i in 1:nrow(dffullint)){
        Yi<-dffullint$Y[i]
        Ai<-dffullint$A[i]
        Si<-dffullint$S[i]
        valstnoi<-valstno[i]
        valsonoi<-valsono[i]
        dfcovfulli<-dfcovfull[i,]
        dfcovfullopi<-dfcovfullop[i,]

        phi1noi<-dfcovfulli*(Si-valsonoi)
        dui<-dfcovfulli*valsonoi*(1-valsonoi)
        A11noi<-dfcovfulli%*%t(dui)

        if(set1==T){
          phi2noih<-Yi*Ai*Si*valstnoi-exp1h*Ai*Si*valstnoi
          phi3noih<-Yi*(1-Ai)*Si*valstnoi-exp0h*(1-Ai)*Si*valstnoi
          A21noih<-dfcovfullopi*Yi*Ai*Si*valstnoi*(1-valstnoi)-exp1h*dfcovfullopi*Ai*Si*valstnoi*(1-valstnoi)
          A31noih<-dfcovfullopi*Yi*(1-Ai)*Si*valstnoi*(1-valstnoi)-exp0h*dfcovfullopi*(1-Ai)*Si*valstnoi*(1-valstnoi)
          A22noih<-Ai*Si*valstnoi
          A33noih<-(1-Ai)*Si*valstnoi

          phimatnoih<-c(phi1noi,phi2noih,phi3noih)%*%t(c(phi1noi,phi2noih,phi3noih))
          phimatnoh<-phimatnoh+phimatnoih

          Anoih<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)

          Anoih[1:ncov,1:ncov]<-A11noi
          Anoih[ncov+1,1:ncov]<-A21noih
          Anoih[ncov+2,1:ncov]<-A31noih
          Anoih[ncov+1,ncov+1]<-A22noih
          Anoih[ncov+2,ncov+2]<-A33noih

          Anoh<-Anoh+Anoih}

        if(set2==T){
          phi2noij<-Yi*Ai*Si*valstnoi/valsonoi-exp1j*Ai*Si*valstnoi/valsonoi
          phi3noij<-Yi*(1-Ai)*Si-exp0j*(1-Ai)*Si
          A21noij<-dfcovfullopi*Yi*Ai*Si*valstnoi*(1-valstnoi)/valsonoi-dfcovfulli*Yi*Ai*Si*valstnoi*(1-valsonoi)/valsonoi
          -exp1j*dfcovfullopi*Ai*Si*valstnoi*(1-valstnoi)/valsonoi-dfcovfulli*Ai*Si*valstnoi*(1-valsonoi)/valsonoi
          A31noij<-rep(0,ncov)
          A22noij<-Ai*Si*valstnoi/valsonoi
          A33noij<-(1-Ai)*Si

          phimatnoij<-c(phi1noi,phi2noij,phi3noij)%*%t(c(phi1noi,phi2noij,phi3noij))
          phimatnoj<-phimatnoj+phimatnoij

          Anoij<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)

          Anoij[1:ncov,1:ncov]<-A11noi
          Anoij[ncov+1,1:ncov]<-A21noij
          Anoij[ncov+2,1:ncov]<-A31noij
          Anoij[ncov+1,ncov+1]<-A22noij
          Anoij[ncov+2,ncov+2]<-A33noij

          Anoj<-Anoj+Anoij}
      }}

    if(crobust==T){
      ncov<-ncol(df[,-c(2,3,4)])+1
      #empty m-estimator inner matrix, h=set1, j=set2
      phimatnoh<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      phimatnoj<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      #empty m-estimator outer matrix, before inversion
      Anoh<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)
      Anoj<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)

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
        #estimated values under each treamtent received vs alternative
        valstnoi<-valstno[clust1]
        valsonoi<-valsono[clust1]
        nvar<-ncol(dfcov)
        Yi<-dfint$Y
        Ai<-dfint$A
        Si<-dfint$S
        phi1noi<-t(t(dfcov)%*%(Si-valsonoi))
        A11noi<-t(dfcov*valsonoi*(1-valsonoi))%*%dfcov

        if(set1==T){
          phi2noih<-sum(Yi*Ai*Si*valstnoi)-exp1h*sum(Ai*Si*valstnoi)
          phi3noih<-sum(Yi*(1-Ai)*Si*valstnoi)-exp0h*sum((1-Ai)*Si*valstnoi)
          A21noih<-apply(dfcovop*Yi*Ai*Si*valstnoi*(1-valstnoi),2,sum)-exp1h*apply(dfcovop*Ai*Si*valstnoi*(1-valstnoi),2,sum)
          A31noih<-apply(dfcovop*Yi*(1-Ai)*Si*valstnoi*(1-valstnoi),2,sum)-exp0h*apply(dfcovop*(1-Ai)*Si*valstnoi*(1-valstnoi),2,sum)
          A22noih<-sum(Ai*Si*valstnoi)
          A33noih<-sum((1-Ai)*Si*valstnoi)

          phimatnoih<-c(phi1noi,phi2noih,phi3noih)%*%t(c(phi1noi,phi2noih,phi3noih))
          phimatnoh<-phimatnoh+phimatnoih

          Anoih<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)

          Anoih[1:ncov,1:ncov]<-A11noi
          Anoih[ncov+1,1:ncov]<-A21noih
          Anoih[ncov+2,1:ncov]<-A31noih
          Anoih[ncov+1,ncov+1]<-A22noih
          Anoih[ncov+2,ncov+2]<-A33noih

          Anoh<-Anoh+Anoih}

        if(set2==T){
          phi2noij<-sum(Yi*Ai*Si*valstnoi/valsonoi)-exp1j*sum(Ai*Si*valstnoi/valsonoi)
          phi3noij<-sum(Yi*(1-Ai)*Si)-exp0j*sum((1-Ai)*Si)
          A21noij<-apply(dfcovop*Yi*Ai*Si*valstnoi*(1-valstnoi)/valsonoi-dfcov*Yi*Ai*Si*valstnoi*(1-valsonoi)/valsonoi,2,sum)-exp1j*apply(dfcovop*Ai*Si*valstnoi*(1-valstnoi)/valsonoi-dfcov*Ai*Si*valstnoi*(1-valsonoi)/valsonoi,2,sum)
          A31noij<-rep(0,ncov)
          A22noij<-sum(Ai*Si*valstnoi/valsonoi)
          A33noij<-sum((1-Ai)*Si)

          phimatnoij<-c(phi1noi,phi2noij,phi3noij)%*%t(c(phi1noi,phi2noij,phi3noij))
          phimatnoj<-phimatnoj+phimatnoij

          Anoij<-matrix(rep(0,(ncov+2)*(ncov+2)),nrow=ncov+2)

          Anoij[1:ncov,1:ncov]<-A11noi
          Anoij[ncov+1,1:ncov]<-A21noij
          Anoij[ncov+2,1:ncov]<-A31noij
          Anoij[ncov+1,ncov+1]<-A22noij
          Anoij[ncov+2,ncov+2]<-A33noij

          Anoj<-Anoj+Anoij}
      }}

    alpha2<-(1-conf)/2
    zl<-qnorm(alpha2)
    zu<--zl
    if(set1==T & set2==F){
      Ainh<-solve(Anoh)
      varsest<-t(c(rep(0,ncov),1,-1))%*%Ainh%*%phimatnoh%*%t(Ainh)%*%c(rep(0,ncov),1,-1)
      lb<-estimators+zl*sqrt(as.numeric(varsest))
      ub<-estimators+zu*sqrt(as.numeric(varsest))
      bounds<-c(lb,ub)
      finnames<-c("EstimateSet1","VarEstSet1","LBS1","UBS1")}

    if(set1==F & set2==T){
      Ainj<-solve(Anoj)
      varsest<-t(c(rep(0,ncov),1,-1))%*%Ainj%*%phimatnoj%*%t(Ainj)%*%c(rep(0,ncov),1,-1)
      lb<-estimators+zl*sqrt(as.numeric(varsest))
      ub<-estimators+zu*sqrt(as.numeric(varsest))
      bounds<-c(lb,ub)
      finnames<-c("EstimateSet2","VarEstSet2","LBS2","UBS2")}

    if(set1==T & set2==T){
      Ainh<-solve(Anoh)
      varsesth<-t(c(rep(0,ncov),1,-1))%*%Ainh%*%phimatnoh%*%t(Ainh)%*%c(rep(0,ncov),1,-1)
      lb1<-estimators[1]+zl*sqrt(as.numeric(varsesth))
      ub1<-estimators[1]+zu*sqrt(as.numeric(varsesth))
      Ainj<-solve(Anoj)
      varsestj<-t(c(rep(0,ncov),1,-1))%*%Ainj%*%phimatnoj%*%t(Ainj)%*%c(rep(0,ncov),1,-1)
      lb2<-estimators[2]+zl*sqrt(as.numeric(varsestj))
      ub2<-estimators[2]+zu*sqrt(as.numeric(varsestj))
      varsest<-c(varsesth,varsestj)
      bounds<-c(lb1,ub1,lb2,ub2)
      finnames<-c("EstimateSet1","EstimateSet2","VarEstSet1","VarEstSet2","LBS1","UBS1","LBS2","UBS2")}
  }

  if(boot==T){
    colboot<-sum(c(set1,set2))
    #empty boot matrix by number of bootstrap iterations
    bootsample<-matrix(rep(0,iters*colboot),ncol=colboot)

    if(crobust==F){
      for(i in 1:iters){
        bootind<-sample(1:nrow(df),replace=T)
        #necessary matrix values for estimation
        dfboot<-dfsim[bootind,]
        resultsb<-suppressMessages(saceestimglm(data=dfboot,trt,surv,out,clustid,indv,set1,set2))
        bootsample[i,]<-unlist(resultsb[-c(1:6)])
      }}

    if(crobust==T){
      for(i in 1:iters){
        bootind<-sample(1:nc,replace=T)
        dflist<-list()
        for(j in 1:nc){
          dflist[[j]]<-df[df$Id==bootind[j],]
        }
        #boostrap data with potentially repeated clusters
        dfboot<-do.call(rbind,dflist)
        #names<-names(df) #why do i need?
        #reestimate parameters on boot data
        resultsb<-suppressMessages(saceestimglm(data=dfboot,trt,surv,out,clustid,indv,set1,set2))
        #generate bootstrap estimates
        bootsample[i,]<-unlist(resultsb[-c(1:6)])
      }
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

  final<-c(estimators,varsest,bounds)
  names(final)<-finnames
  #final<-list(Ainj,phimatnoj)
  return(final)
}
