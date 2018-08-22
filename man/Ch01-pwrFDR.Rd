\name{pwrFDR}
\alias{pwrFDR}
\title{
  Power in the BH-FDR procedure.
}
\description{
  Calculates average power or \eqn{\lambda}-power in the Benjamini-Hochberg
  FDR procedure, or sample size required for given average- or \eqn{\lambda}- power.
}
\usage{
  pwrFDR(groups = 2, effect.size, n.sample, r.1, FDR, use.LFD = FALSE, 
         N.tests, average.power, L.power, lambda, method = c("approximate", 
         "simulation", "JL", "Iz"), control = list(version = 0, 
         tol = 1e-08, max.iter = 1000, distopt = 1, CS = list(NULL), 
         verb = FALSE), n.sim = 1000, temp.file) 
}
\arguments{
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
  \item{use.LFD}{
    Under the default value, 'FALSE', the procedure is controlled
    under the specified FDR.  Specifying "use.LFD=TRUE" will
    instead first find an FDR, FDR.star, so that the probability that
    the FDF exceeds FDR.star is less than or equal to the originally
    specified FDR. Power (either average- or lambda-) or sample size
    under the average- or lambda- power is then calculated when the
    procedure is controlled at FDR.star. This is recommended when the
    distribution of the FDF is wide enough to matter but still well
    approximated by its asymptotic distribution, e.g. between
    several hundred and several thousand simultaneous tests. The
    argument 'N.tests' is required when "use.LFD=TRUE".
  }
  \item{N.tests}{
    The number of simultaneous hypothesis tests. 
  }
  \item{average.power}{
    The desired average power. Sample size calculation requires
    specification of either 'average.power' or 'L.power'.
  }
  \item{L.power}{
    The desired \eqn{\lambda}-power (see details for explanation). Sample
    size calculation requires specification of either 'average.power'
    or 'L.power'.
  }
  \item{lambda}{
    The \eqn{\lambda}-power threshold, required when calculating the
    \eqn{\lambda}-power (see details for explanation) or when calculating
    the sample size required for \eqn{\lambda}-power.    
  }
  \item{method}{
    Specify the method whereby the average power is calculated. The
    selection 'approximate' uses the approximation technique of
    Izmirlian (2016). To use direct simulation, specify 'simulation'.
    Specification of 'JL' uses the method of Jung (2005) and Liu
    (2007). Specification of 'Iz' uses the form which was shown
    (Izmirlian) to be the large number of tests limit of
    the average power. Valid selections are 'Iz' (default), 'JL',
    and any substring of 'approximate' or 'simulation'.
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
  \item{n.sim}{
    If 'simulation' method is chosen you may specify number of
    simulations. Default is 1000.
  }
  \item{temp.file}{
    If 'simulation' method is chosen you may specify a tempfile
    where the current simulation replicate is updated. Very usefull
    for batch runs. You can use the included utility 'gentempfilenm'
  }
}
\details{
  Power for the BH-FDR procedure on 'm' simultaneous tests, at a given
  FDR, f, is computed under the following model. A-priori, a proportion, 
  'r.1', of the 'm' test statistics are given a non-centrally located
  distribution, the remaining proportion are given a centrally
  distribution. Each of the non-centrally located test statistics is
  given the same location parameter. Suppose that a total of 'M'
  are chosen to have a non-central distribution, that a total of 'J'
  of all 'm' statistics are declared significant, and that of these,
  'S' come from the non-centrally located population. This results in
  the following table. 
  \tabular{rrccc}{
    1. \tab             \tab rej H0  \tab acc H0         \tab row Total \cr
    2. \tab H0 is FALSE \tab   S     \tab M-S            \tab M         \cr
    3. \tab H0 is TRUE  \tab J-S     \tab (m-M)-(J-S)    \tab m-M.1     \cr
    4. \tab col Total   \tab   J     \tab m-J            \tab m         \cr
  }
  
  By default, when 'n.sample' is specified, the function computes the 
  the average power,

  AVERAGE POWER:
  \eqn{E[ S/M ]}

  This is approximated using the infinite tests limit.

  When 'n.sample' and '\eqn{\lambda}' are both specified, the function computes
  the \eqn{\lambda}-power, which is the probability:

  \eqn{\lambda}-POWER:
  \eqn{P\{ S/M > \lambda \}}
  
  that the true positive fraction S/M exceeds the threshold, '\eqn{\lambda}'.
  This is done using an asymptotic approximation. If the user 
  desires instead the sample size required for given average- or \eqn{\lambda}-
  power, then leave the argument 'n.sample' unspecified and specify
  instead either the 'average.power' or the '\eqn{\lambda}', the threshold for
  the \eqn{\lambda}-power. By default, the model uses the assumption that test
  statistics are distributed as t-distributed as in the two group
  comparison. If the user specifies a value of 'groups' larger than 2
  then the model assumes the test statistics are F-distributed as in
  the omnibus F-test for any difference. The above approximations
  are done under the default method. If the user wants to compare with
  simulated answers to obtain simulation estimates of the 'exact'
  values, one need only specify 'method="sim"'. 
}
\value{
  An object of class "pwr" with with components including:
  \item{call}{The call which produced the result}
  \item{average.power}{Resulting average power.}
  \item{L.power}{When 'lambda' is specified, the lambda-power is also computed}
  \item{L.eq}{The lambda at which the lambda-power and average-power are equal.}
  \item{n.sample}{If 'n.sample' is missing from the argument list, then
    the sample size required for the specified average- or lambda-
    power.}
  \item{FDR.star}{If "use.LFD=TRUE" was specified, the FDR at which
    the probability that the FDF exceeds FDR.star is less than or equal
    to the originally specified FDR.}
  \item{c.g}{The BH-FDR threshold on the scale of the test statistics.}
  \item{gamma}{The proportion of all 'm' tests declared significant.}
  \item{objective}{Result of optimization yielding the 'average.power'.}
  \item{err.III}{Mass on the wrong side of the threshold.}
  \item{sigma.rtm.SoM}{Asymptotic variance of the true positive fraction.}
} 
\references{
  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. arXiv:1801.03989

  Jung S-H. (2005) Sample size for FDR-control in microarray data
    analysis. Bioinformatics; 21:3097-3104.
  
  Liu P. and Hwang J-T. G. (2007) Quick calculation for sample size while
    controlling false discovery rate with application to microarray
    analysis. Bioinformatics; 23:739-746.
}
\author{
  Grant Izmirlian <izmirlian at nih dot gov>
}
\seealso{
  \code{\link{controlFDF}}
}
\examples{
## Example 1a: average power
   rslt.avgp <- pwrFDR(effect.size=0.79, n.sample=46, r.1=2000/54675, FDR=0.15)
   rslt.avgp

## Example 1b: lambda-power
   rslt.lpwr <- pwrFDR(effect.size=0.79, n.sample=46, r.1=2000/54675,
                       FDR=0.15, lambda=0.80, N.tests=54675)
   rslt.lpwr

## Example 1c: sample size required for given average power
   rslt.ss.avgp <- pwrFDR(effect.size=0.79, average.power=0.82,
                          r.1=2000/54675, FDR=0.15)
   rslt.ss.avgp

## Example 1d: sample size required for given lambda-power
   rslt.ss.lpwr <- pwrFDR(effect.size=0.79, L.power=0.82, lambda=0.80,
                          r.1=2000/54675, FDR=0.15, N.tests=54675)
   rslt.ss.lpwr

## Example 1e: simulation
   rslt.sim <- update(rslt.avgp, method="sim", n.sim=500, N.tests=1000)
   rslt.sim

## Example 2: methods for adding, subtracting, multiplying, dividing, exp, log,
## logit and inverse logit
   rslt.avgp - rslt.sim
   logit(rslt.avgp)       ## etc
   
## Example 3: Compare the asymptotic distribution of S/M with kernel
##            density estimate from simulated data 
   pdf <- with(detail(rslt.sim)$reps, density(S/M1))

   med <- with(detail(rslt.sim)$reps, median(S/M1))
   avg <- rslt.sim$average.power
   sd <- rslt.sim$v.SoM.emp^0.5

   rng.x <- range(pdf$x)
   rng.y <- range(c(pdf$y, dnorm(pdf$x, mean=avg, sd=sd)))

   plot(rng.x, rng.y, xlab="u", ylab="PDF for S/M", type="n")
   with(pdf, lines(x, y))
   lines(rep(rslt.sim$average.power, 2), rng.y, lty=2)
   lines(pdf$x, dnorm(pdf$x, mean=avg, sd=sd), lty=3)
}
\keyword{FDR}
\keyword{Benjamini}
\keyword{Hochberg}
\keyword{microarrays}
\keyword{Multiple.Testing}
\keyword{average.power}
\keyword{k.power}
\keyword{lambda.power}
