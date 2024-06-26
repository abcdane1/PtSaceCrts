% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saceestimwvar.R
\name{sacecluster}
\alias{sacecluster}
\title{Estimation of SACE in Cluster Randomized Trials}
\usage{
sacecluster(
  df,
  trt = "A",
  surv = "S",
  out = "Y",
  clustid = "Id",
  csize = "Csize",
  clustv = "C",
  indv = "X",
  set1 = T,
  set2 = F,
  ics = F,
  boot = F,
  nagq = 10,
  iters = 200
)
}
\arguments{
\item{df}{Data frame containing all data required for estimation with names listed below.}

\item{trt}{A named \code{character} specifying treatment variable. Default is "A".}

\item{surv}{A named \code{character} specifying survival status. Default is "S".}

\item{out}{A named \code{character} specifying non-mortal outcome. Default is "Y".}

\item{clustid}{A named \code{character} specifying non-mortal cluster membership. Default is "Id".}

\item{csize}{A named \code{character} specifying cluster size. Default is "Csize". This need not be specified if \code{ics=F}.}

\item{clustv}{A named \code{character} vector for cluster-level variables. Default is "C"}

\item{indv}{A named \code{character} vector for individual-level varibles. Default is "X".}

\item{set1}{A \code{logical} argument for whether identified estimator uses Set 1 Assumptions. Default is \code{T}.}

\item{set2}{A \code{logical} argument for whether identified estimator uses Set 2 Assumptions. Default is \code{F}.
If \code{set2=T} and \code{set2=T}, function will provide out for both estimators}

\item{ics}{A \code{logical} argument for whether estimator should account for cluster size. Default is \code{F}.}

\item{boot}{A \code{logical} argument for variance estimation. If \code{boot=F}, variance is estimated
using the asymptotic variance of estimator. If \code{boot=T}, variance of estimator is
computed using nonparametric boostrap. Default is \code{F}.}

\item{nagq}{A \code{double} for number of quadrature points used in the AGQ approximation of integrals
involved in computing estimate of asymptotic variance. Default is \code{nagq=10}. It is not recommended
to use in excess of 20 quadrature points. If \code{boot=T}, this argument is unused.}

\item{iters}{A \code{double} for number of bootstrap samples to be taken. Default is \code{iters=200}.}
}
\value{
A named \code{double} including point estimates, estimates of variance, and confidence intervals.
To add: names, CIs.
}
\description{
This function provides you point estimates of SACE in cluster randomized
trials based on estimators derived from two sets of identification assumptions
as described in (cite paper). All principal score based models are fit using GLMM for binary data.
The function also provides estimates of the variance for these estimators. User has the option
to choose how variance is estimated and corresponding confidence intervals
are computed. By default, the function provides the variance of the asymptotic distribution of
these estimators (cite wy), which relies on numeric approximation of integrals via Adaptive Gauss-Hermite
Quadrature (AGQ) (cite). The user may also select to estimate variance using non-parametric
cluster bootstrap option (cite welsh).
}
\references{
\enumerate{
\item{add.}
\item{add.}
\item{add.}
}
}
