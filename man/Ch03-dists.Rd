\name{dists}
\alias{dists}
\docType{data}
\title{
  The Distribution family object
}
\description{
  The \code{pwrFDR} package currently incorporates 3 distribution types,
  normal, t and F. The first two of these are strictly for statistics formed 
  from two group comparison while the third is for statistics formed from the 
  omnibus test of any difference among an arbitrary number of groups >=2. The
  structure is general and user expandable. One must specify the density, 
  CDF and quantile function for a given distribution and its parameters under
  the null and under the alternative. These parameters must be expressions
  to be evaluated inside the kernel of the power program, functions of the
  arguments \code{n.sample}, \code{groups} and \code{effect.size}. This is 
  not used directly by the user at all unless she (he) wants to add a 
  distribution type. 
}
\format{
  A data frame with 3 observations on the following 6 variables.
  \describe{
    \item{\code{pars0}}{a list vector having components 'c(nd, p1, p2, ...)'
                        where 'nd' is the distribution number starting with 0,
                        and p1, p2, ..., are paramters of the distribution,
                        which are functions of 'n.sample', 'groups' and 
                        'effect.size' as mentioned above.  These must be 
                        expressed as a call e.g.  
                        as.call(expression(c, nd, p1, p2, ...))  etc.
                        'pars0' are the parameters under the null.}
    \item{\code{pars1}}{a list vector. See directly above. Parameters under the 
                        alternative.}
    \item{\code{minv}}{a list vector with components given the values -Inf or 0,
                       which will be used to decide if the two sided corrections
                       are used or not.}
    \item{\code{ddist}}{a list vector with components set to functions, each
                        one computing the probability density function corresponding
                        to the particular distribution. A function of arguments 'x' 
                        and 'par'. See details below.}
    \item{\code{pdist}}{a list vector with components set to the functions, each
                        one computing the cumulative distribution function corresponding
                        to the particular distribution. A function of arguments 'x' 
                        and 'par'. See details below.}
    \item{\code{qdist}}{a list vector with components set to the functions, each
                        one computing the quantile function (inverse cumulative 
                        distribution function) corresponding to the particular 
                        distribution. A function of arguments 'x' and 'par'. 
                        See details below.}
  }
}
\details{
   The instance shipped with the current version of the package was created using
   the following commands. Read this and you'll figure out how to add a distribution
   type yourself.
\cr
"dists" <-\cr
as.data.frame(rbind(\cr
                    \verb{### Normal with 2 groups ###}\cr
                      c(pars0=as.call(expression(c,0,ncp=0,sd=1)),\cr
			pars1=as.call(expression(c,0,ncp=(n.sample/2)^0.5*effect.size, sd=1)),\cr
                        minv=-Inf,\cr
                        ddist=function(x, par) dnorm(x, mean=par[2], sd=par[3]),\cr
                        pdist=function(x, par) pnorm(x, mean=par[2], sd=par[3]),\cr
                        qdist=function(x, par) qnorm(x, mean=par[2], sd=par[3])),\cr\cr
		      \verb{### t with 2 groups ###}\cr
                      c(pars0=as.call(expression(c,1,ncp=0, ndf=2*n.sample - 2)),\cr
                        pars1=as.call(expression(c,1,ncp=(n.sample/2)^0.5*effect.size, ndf=2*n.sample - 2)),\cr
                        minv=-Inf,\cr
                        ddist=function(x, par) dt(x, ncp=par[2], df=par[3]),\cr
                        pdist=function(x, par) pt(x, ncp=par[2], df=par[3]),\cr
                        qdist=function(x, par) qt(x, ncp=par[2], df=par[3])),\cr\cr
		      \verb{### F with 'groups' groups, effect.size=theta*c(0, 0.5, 0.5, ..., 0.5, 1) ###}\cr
                      c(pars0=as.call(expression(c,2,ncp=0, ndf1=groups-1, ndf2=groups*(n.sample-1))),\cr
                        pars1=as.call(expression(c,2,ncp=n.sample*effect.size^2/2, ndf1=groups-1, ndf2=groups*(n.sample-1))),\cr
                        minv=0,\cr
                        ddist=function(x, par) df(x, ncp=par[2], df1=par[3], df2=par[4]),\cr
                        pdist=function(x, par) pf(x, ncp=par[2], df1=par[3], df2=par[4]),\cr
                        qdist=function(x, par) qf(x, ncp=par[2], df1=par[3], df2=par[4]))))
}
\source{
  This isn't 'data' data, its a kind of a 'family' object.

  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. arXiv:1801.03989
}
