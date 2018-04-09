\name{cCDF.ToJ}
\alias{cCDF.ToJ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Complementary CDF for the False Discovery Fraction
}
\description{
  Computes the complementary CDF for the False Discovery Fraction, T_m/J_m, via
  approximation based upon the asymptotic distribution.
}
\usage{
  cCDF.ToJ(lambda, x=NULL, groups, effect.size, n.sample, r.1, FDR, N.tests, 
           control)
}
\arguments{
  \item{lambda}{
    Argument of the complementary CDF. Result will be Pr( S_m/M_m > lambda )
  }
  \item{x}{
    In the abbreviated call sequence, the user only needs to specify \code{lambda},
    \code{N.tests}, and \code{x}, an object of class \code{pwr} returned by the 
    function \code{pwrFDR}. 
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
    The number of experimental replicates.
  }
  \item{r.1}{
    The proportion of all test statistics that are distributed under HA.
  }
  \item{FDR}{
    The false discovery rate.
  }
  \item{N.tests}{
    The number of simultaneous hypothesis tests. 
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
  The complementary CDF for the False Discovery Fraction, T_m/J_m, 
  is approximated using its asymptotic distribution. Since

  m^0.5 ( T_m/J_m - (1 - r.1) FDR ) --D--> N(0, tau^2)

  then 

  P( T_m/J_m > lambda) 

    ~  1 - Phi( m^0.5 (lambda - (1 - r.1) FDR)/tau )

  The approximation is reasonable as long as m tau^2 is large enough
  A formula for the asymptotic variance is given in the cited manuscript.
  There is a user level function, var.rtm.ToJ, which computes the asymptotic 
  variance.  
}
\value{
  An object of class "vvv" which is a list having components
    \item{cCDF.ToJ}{The result}
    \item{average.power}{The average power at the supplied arguments}
    \item{c.g}{The per test threshold that is equivalent to the BH-FDR, on the 
               test statistic scale}
    \item{gamma}{The limiting proportion of tests that were called significant}
    \item{objective}{Result of optimization producing the average power, should
                     be close to zero.}
    \item{err.III}{The probability mass on side of the oppositely signed alternative
                   in two sided tests}
    \item{sigma.rtm.ToJ}{The square root of the asymptotic variance of the root-m 
                         scaled false discovery fraction, m^0.5 * T_m/J_m}
    \item{call}{The call which produced the result}
}
\references{
  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. arXiv:1801.03989
}
\author{
  Grant Izmirlian <izmirlian at nih dot gov>
}
\seealso{
  \code{\link{controlFDF}}
  \code{\link{cCDF.SoM}}
}
\examples{
## Example 1: Explicit call 
   ccdf <- cCDF.ToJ(lambda=(1-2000/54675)*0.15, effect.size=0.79, n.sample=46, r.1=2000/54675, 
                    FDR=0.15, N.tests=1000)
   ccdf

## Example 2: Abbreviated call using result of pwrFDR
   rslt.avgp <- pwrFDR(effect.size=0.79, n.sample=46, r.1=2000/54675, FDR=0.15)
   ccdf <- cCDF.ToJ(lambda=(1-2000/54675)*0.15, x=rslt.avgp, N.tests=1000)
   ccdf
}
\keyword{FDR}
\keyword{Benjamini}
\keyword{Hochberg}
\keyword{microarrays}
\keyword{Multiple.Testing}
\keyword{average.power}
\keyword{k.power}
\keyword{lambda.power}
