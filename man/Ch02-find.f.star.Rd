\name{find.f.star}
\alias{find.f.star}
\title{
  Calculates a reduced FDR required to control the FDF
}
\description{
  Calculates a reduced FDR required to bound the the false 
  discovery rate in probability. 
}
\usage{
  find.f.star(groups, FDR, r.1, N.tests, effect.size, n.sample)
}
\arguments{
  \item{groups}{
    The number of experimental groups to compare. Default value is 2.
  }
  \item{FDR}{
    The false discovery rate.
  }
  \item{r.1}{
    The proportion of simultaneous tests that are non-centrally located 
  }
  \item{N.tests}{
    The number of simultaneous hypothesis tests. 
  }
  \item{effect.size}{
    The effect size (mean over standard deviation) for test statistics 
    having non-zero means. Assumed to be a constant (in magnitude) over
    non-zero mean test statistics.
  }
  \item{n.sample}{
    The number of experimental replicates.
  }
}
\details{
  Calculates a reduced FDR required to bound the the false 
  discovery rate in probability...e.g. finds f* so that
  when the BH-FDR procedure at FDR, f*, is used, we ensure
  that

  Pr( T/J > (1-r) f ) < (1-r) f*

  where 'f' is the original false discovery rate and
  'r' is the proportion of non-null distributed test
  statistics.
}
\value{
  \item{f.star}{The reduced FDR required to bound the FDF in probability}
  \item{obj}{Result of optimization yielding 'f.star' should be close to 0}
  \item{L.star}{The bound on the FDF, should be (1-r) f.  See above.}
  \item{P.star}{The probability that the FDF is greater than L.star. See above.}
  \item{average.power}{Resulting average power.}
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
}
\author{
  Grant Izmirlian <izmirlian at nih dot gov>
}
\seealso{
  \code{\link{pwrFDR}}
}
\examples{
## at FDR=0.15 and other parameters, it takes n.sample=46 replicates for 
## average power > 80%
pwr.46.15 <- pwrFDR(FDR=0.15, r.1=0.03, N.tests=1000, effect.size=0.79, n.sample=46)

## when there are 'only' N.tests=1000 simultaneous tests, the distribution of the
## false discovery fraction, FDF, is not so highly spiked at the FDR=0.15
## You need to set the FDR down to FDR=0.0657 to ensure that  Pr( T/J > 0.145 ) < 0.0657
fstr <- find.f.star(FDR=0.15, r.1=0.03, N.tests=1000, effect.size=0.8, n.sample=46)

## at all the above settings, with FDR=0.0657 at an n.sample of 46, we only have 69% 
## average power.
pwr.46.0657 <- pwrFDR(FDR=0.065747, r.1=0.03, N.tests=1000, effect.size=0.79, n.sample=46)

## it'll cost 7 more replicates to get the average power up over 80%.
pwr.53.0657 <- pwrFDR(FDR=0.065747, r.1=0.03, N.tests=1000, effect.size=0.8, n.sample=53)

## it costs only 8.75% more to get it right!
}
\keyword{FDR}
\keyword{Benjamini}
\keyword{Hochberg}
\keyword{microarrays}
\keyword{Multiple.Testing}
\keyword{average.power}
\keyword{k.power}
\keyword{lambda.power}
