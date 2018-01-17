\name{var.J.o.rtm}
\alias{var.J.o.rtm}
\title{
  Calculates asymptotic variance of significant call fraction, J_m/m
  in the BH-FDR procedure on m=N.tests simultaneous tests.
}
\description{
  A function which computes the asymptotic variance of the proportion
  of significant calls, J_N/N.  
}
\usage{
  var.J.o.rtm(x, groups, effect.size, n.sample, r.1, FDR, N.tests, control)
}
\arguments{
  \item{x}{
    Calls to this function can be made either specifying the single
    argument, x, which is an object of class "pwr" returned from the
    function, \code{pwrFDR} or by a full call specification, by
    specifying all the arguments required by the function \code{pwrFDR}
    (see below)
  }
  \item{groups}{
    The number of experimental groups to compare. Default value is 2.
  }  
  \item{effect.size}{
    The effect size (mean over standard deviation) for test statistics 
    having non-zero means. Assumed to be a constant (in magnitude) over
    non-zero mean test statistics.
  }
  \item{n.sample}{
    The number of experimental replicates. Required for calculation
    of average.power or L.power.
  }
  \item{r.1}{
    The proportion of simultaneous tests that are non-centrally located 
  } 
  \item{FDR}{
    the false discovery rate.
  }
  \item{N.tests}{
    The number of simultaneous tests. 
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
  The proportion of significant calls, J_m/m, i.e., the proportion
  of all hypothesis tests declared significant by the Benjamini-Hochberg
  procedure is shown in the cited publication to be root-m consistent
  and asymptotically normal i.e.
  
  m^{1/2} ( J_m/m - gamma )  --D--> N(0, tau^2)
  
  which is to say that the distribution of J_m/m is spiked about its
  mean, \code{gamma}, and the width dies off as 1/m^{1/2}. This is of
  interest in its own right, and is also needed as an intermediate in
  the calculation of the asymptotic variance of the empirical average
  power, S_m/M_m.  The normal approximation for this latter quantity
  can be used to power multiple testing experiments on a lower quantile
  of the empirical mean rather than on its expected value. The width of
  this distribution is non-negligible for as many as 200 simultaneous
  tests, where the average power could be 80\% but the lower quantile
  of the empirical average power could be as low as 50\%.   
}
\value{
  Returns a value of class \code{vvv}, containing components
  \item{var.J.o.rtm}{The computed asymptotic variance}
  \item{power}{The average power}
  \item{gamma}{The expected proportion of significant calls}
  \item{c.g}{The 'q-value', which is the value of the criterion on the
    scale of the statistic (t of given number of 2 n - 2 degrees of
    freedom) which can be used as a per test criterion resulting in
    the equivalent Benjamini-Hochberg procedure}
  \item{call}{The call which produced the result}
}
\references{
  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. arXiv:1801.03989
}
\author{
  Grant Izmirlian <izmirlig at mail dot nih dot gov>
}
\seealso{
  \code{\link{var.rtm.SoM}}
  \code{\link{var.rtm.ToJ}}
}
\examples{
  ## call using result of pwrFDR 
  rslt.Iz <- pwrFDR(effect.size=0.79, n.sample=46, r.1=2000/54675, FDR=0.15)
  vJ <- var.J.o.rtm(rslt.Iz)

  ## call via argument list specification
  vJ <- var.J.o.rtm(effect.size=0.79, n.sample=46, r.1=2000/54675, FDR=0.15)
}
\keyword{Benjamini-Hochberg}
\keyword{Multiple.Testing}
\keyword{k.power}
\keyword{lambda.power}
