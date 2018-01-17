\name{CDF.Pval}
\alias{CDF.Pval}
\title{
  CDF of pooled (H0 and HA) population p-values
}
\description{
  Computes the CDF of the pooled population p-values under the
  mixture model, e.g. the p-values are i.i.d. with CDF 
  a mixture between a uniform (CDF in the null distributed
  population) and a concave function (CDF in the non-null 
  distributed population).
}
\usage{
  CDF.Pval(u, groups, r.1, effect.size, n.sample, control)
}
\arguments{
  \item{u}{
    Argument of the CDF. Result will be Pr( P_i <= u )
  }
  \item{groups}{
    The number of experimental groups to compare. Default value is 2.
  }
  \item{r.1}{
    The proportion of all test statistics that are distributed under HA.
  }
  \item{effect.size}{
    The effect size (mean over standard deviation) for test statistics 
    having non-zero means. Assumed to be a constant (in magnitude) over
    non-zero mean test statistics.
  }
  \item{n.sample}{
    The number of experimental replicates.
  }
  \item{control}{
    Optionally, a list with components with the following components:
    'groups', used when distop=3 (F-dist), specifying number of groups. 
    'version', used only in the 'JL' method, choice 0 gives the 'JL'
    version as published, whereas choice 1 replaces the FDR with
    r.0*FDR resulting in the infinite simultaneous tests limiting
    average power, which is the 'Iz' version, but this is redundant
    because you can specify the 'Iz' method to use this option.
    'tol' is a convergence criterion used in iterative methods
    which is set to 1e-8 by default
    'max.iter' is an iteration limit, set to 1000 by default
    'distop', specifying the distribution family of the central and
      non-centrally located sub-populations. =1 gives normal (2 groups)
      =2 gives t- (2 groups) and =3 gives F- (2+ groups)    
    'CS', correlation structure, for use only with 'method="simulation"'
      which will simulate m simulatenous tests with correlations 'rho'
      in blocks of size 'n.WC'. Specify as list
      CS = list(rho=0.80, n.WC=50) for example
  }
}
\details{
  Computes the CDF of the pooled population p-values under the
  mixture model, e.g. the p-values are i.i.d. with CDF 
  a mixture between a uniform (CDF in the null distributed
  population) and a concave function (CDF in the non-null 
  distributed population). If Fc_0 is the cCDF of a test 
  statistic under H0 and Fc_A is the cCDF of a test statistic
  under HA then the CDF of the P-values is
  
  G(u) = (1-r) u  + r Fc_A(Fc_0^{-1}(u))

  The limiting positve call fraction, lim_m J_m/m = gamma (a.s.)
  is the solution to the equation 

  G( gamma f) = gamma 

  where f is the FDR
}
\value{
  A list with components
  \item{call}{The call which produced the result}
  \item{u}{The argument that was passed to the function}
  \item{CDF.Pval}{The value of the CDF}
}
\references{
  Genovese, C. and L. Wasserman. (2004) A stochastic process approach to
  false discovery control. Annals of Statistics. 32 (3), 1035-1061.

  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. arXiv:1801.03989
}
\author{
  Grant Izmirlian <izmirlian at nih dot gov>
}
\seealso{
  \code{\link{CDF.Pval.HA}}
}
\examples{
  ## First calculate an average power for a given set of parameters
  rslt.avgp <- pwrFDR(effect.size=0.79, n.sample=46, r.1=2000/54675, FDR=0.15)

  ## Now verify that G( gamma f ) = gamma

  gamma <- rslt.avgp$gamma
  f <- rslt.avgp$call$FDR

  G.gma.f <- CDF.Pval(u=gamma*f, r.1=2000/54675, effect.size=0.79, n.sample=46)

  c(G.of.gamma.f=G.gma.f$CDF.Pval, gamma=gamma) 
}
\keyword{FDR}
\keyword{Benjamini}
\keyword{Hochberg}
\keyword{microarrays}
\keyword{Multiple.Testing}
\keyword{average.power}
\keyword{k.power}
\keyword{lambda.power}
