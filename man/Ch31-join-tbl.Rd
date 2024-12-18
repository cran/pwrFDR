\name{join.tbl}
\alias{join.tbl}
\title{
  Combine pwrFDR Results
}
\description{
  Joins pwrFDR objects into a single table.
}
\usage{
join.tbl(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
        obj1, obj2, ... each being the result of a call to pwrFDR.
        See the example below.
}
}
\value{
  The table of joined pwrFDR objects as a data.frame 
}
\author{
  Grant Izmirlian
}
\examples{
  rslt.avgp.r15 <- pwrFDR(effect.size = 0.79, average.power=0.80, r.1 = 0.15, alpha = 0.15)
  rslt.avgp.r10 <- update(rslt.avgp.r15, r.1 = 0.10)
  rslt.avgp.r05 <- update(rslt.avgp.r15, r.1 = 0.05)

  join.tbl(rslt.avgp.r15, rslt.avgp.r10, rslt.avgp.r05)
}
