\name{var.rtm.ToJ}
\alias{var.rtm.ToJ}
\title{
  Calculates asymptotic variance of false discovery fraction, T_N/J_N.
}
\description{
  A function which computes the asymptotic variance of the false
  discovery fraction, T_N/J_N.
}
\usage{
  var.rtm.ToJ(x, groups, effect.size, n.sample, r.1, FDR, N.tests, control)
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
    of power
  }
  \item{r.1}{
    The proportion of simultaneous tests that are non-centrally located 
  } 
  \item{FDR}{
    The false discovery rate.
  }
  \item{N.tests}{
    Number of simultaneous tests.
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
  The false discovery fraction (FDF), T_m/J_m, i.e., is the proportion
  of all centrally located statistical tests that are declared
  significant by the Benjamini-Hochberg procedure. It is shown, in the
  cited publication, to be root-m consistent and asymptotically normal i.e.
  
  m^{1/2} ( T_m/J_m - (1-r.1) f )  --D--> N(0, sigma^2)
  
  which is to say that the distribution of T_m/J_m is spiked about its
  mean, the false discovery rate, and the width dies off as
  1/m^{1/2}. This noraml approximation is of interest in its own right,
  allowing tighter control over the false discovery fraction that
  afforded by controlling its mean, the false discovery rate.
  One can determine the FDR required to bound the FDF by a given
  value with a probability no less than a given value.
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
  \code{\link{var.J.o.rtm}}
  \code{\link{var.rtm.SoM}}
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
\keyword{average.power}
\keyword{k.power}
\keyword{lambda.power}
