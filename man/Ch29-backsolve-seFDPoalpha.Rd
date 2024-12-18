\name{backsolve.seFDPoalpha}
\alias{backsolve.seFDPoalpha}
\alias{backsolve.seTPPoavgpwr}
\title{
Find missing argument giving required se[FDP]/alpha (or se[TPP]/average.power)
}
\description{
  backsolve.seFDPoalpha finds the missing argument, one of 'N.tests',
  'r.1', 'n.sample' or 'effect size' giving the specified value of
  se[FDP]/alpha under the BH-FDR procedure.
  
  backsolve.seTPPoavgpwr finds the missing argument, one of 'N.tests',
  'r.1', 'n.sample' or 'effect size' giving the specified value of
  se[TPP]/average.power under the BH-FDR procedure.
}
\usage{
backsolve.seFDPoalpha(seFDPoalpha, effect.size, n.sample, r.1, alpha, groups = 2, N.tests,
                      type = "balanced", grpj.per.grp1 = 1, distopt = 1, rho, k.bs)

backsolve.seTPPoavgpwr(seTPPoavgpwr, effect.size, n.sample, r.1, alpha, groups = 2,
                       N.tests, type = "balanced", grpj.per.grp1 = 1, distopt = 1, rho,
                       k.bs)
}
\arguments{
  \item{seFDPoalpha}{
    In backsolve.seFDPoalpha, the user specified value of se[FDP]/alpha
  }
  \item{seTPPoavgpwr}{
    In backsolve.seTPPoavgpwr, the user specified value of se[TPP]/average.power
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
  \item{alpha}{
    The false discovery rate (in the BH case) or the upper bound on the
    probability that the FDP exceeds delta (BHFDX and Romano case)
  }
  \item{groups}{
    The number of experimental groups to compare. Must be integral and
    >=1. The default value is 2.
  }
  \item{N.tests}{
    The number of simultaneous hypothesis tests. 
  }
  \item{type}{
    A character string specifying, in the groups=2 case, whether the
    test is 'paired', 'balanced', or 'unbalanced' and in the case when
    groups >=3, whether the test is 'balanced' or 'unbalanced'. The
    default in all cases is 'balanced'. Left unspecified in the one
    sample (groups=1) case. 
  }
  \item{grpj.per.grp1}{
    Required when \code{type}="unbalanced", specifies the group 0 to
    group 1 ratio in the two group case, and in the case of 3 or more
    groups, the group j to group 1 ratio, where group 1 is the group
    with the largest effect under the alternative hypothesis. 
  }
  \item{distopt}{
    Test statistic distribution in among null and alternatively
    distributed sub-populations. distopt=0 gives normal (2 groups),
    distop=1 gives t- (2 groups) and distopt=2 gives F- (2+ groups)
  }
  \item{rho}{
    This can be done under the assumption of tests that are correlated
    identically in pair within blocks of given size.
  }
  \item{k.bs}{
    When 'rho' is specified, the common block-size for correlated test
    statistics. 
  }
}
\value{
A numeric vector having components
  \item{<missing argument>}{Value of missing argument giving required
  se[FDP]/alpha (backsolve.seFDPoalpha) or se[TPP]/average.power
    (backsolve.seTPPoavgpwr).}
  \item{average.power}{The average power at the given set of conditions}
  \item{se.VoR/se.ToM}{The standard error of the FDP
    (backsolve.seFDPoalpha) or standard error of the TPP 
    (backsolve.seTPPoavgpwr).}
  \item{value}{Value returned by the solver. Should be near zero if a
    solution was found.}
}
\references{
  Izmirlian G. (2020) Strong consistency and asymptotic normality for
    quantities related to the Benjamini-Hochberg false discovery rate
    procedure. Statistics and Probability Letters; 108713,
    <doi:10.1016/j.spl.2020.108713>

  Izmirlian G. (2017) Average Power and \eqn{\lambda}-power in
    Multiple Testing Scenarios when the Benjamini-Hochberg False
    Discovery Rate Procedure is Used. <arXiv:1801.03989>

  Jung S-H. (2005) Sample size for FDR-control in microarray data
    analysis. Bioinformatics; 21:3097-3104.
  
  Kluger D. M., Owen A. B. (2023) A central limit theorem for the
    Benjamini-Hochberg false discovery proportion under a factor model.
    Bernoulli; xx:xxx-xxx.

  Liu P. and Hwang J-T. G. (2007) Quick calculation for sample size while
    controlling false discovery rate with application to microarray
    analysis. Bioinformatics; 23:739-746.

  Lehmann E. L., Romano J. P.. Generalizations of the familywise error
    rate. Ann. Stat.. 2005;33(3):1138-1154.

  Romano Joseph P., Shaikh Azeem M.. Stepup procedures for control of
    generalizations of the familywise error rate. Ann. Stat..
    2006;34(4):1850-1873. 
}
\author{
  Grant Izmirlian Jr <izmirlian at nih dot gov>
}
\examples{
backsolve.seFDPoalpha(seFDPoalpha=0.50, n.sample=50, alpha=0.05, effect.size=0.8,
                      r.1=0.20)

backsolve.seTPPoavgpwr(seTPPoavgpwr=0.20, n.sample=30, alpha=0.05, effect.size=0.8,
                       r.1=0.20)
}
\keyword{FDR}
\keyword{Benjamini}
\keyword{Hochberg}
\keyword{FDX}
\keyword{microarrays}
\keyword{Multiple.Testing}
\keyword{average.power}
\keyword{TPX}
\keyword{k.power}
\keyword{TPX.power}
