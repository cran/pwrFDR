\name{detail}
\alias{detail}
\title{
  The detail extraction function for simulated power objects
}
\description{
  Objects created by the \code{pwrFDR} function with option
  \code{method}=="simulation" are returned with an attribute named
  \code{detail}. This is its extractor function
}
\usage{
detail(obj)
}
\arguments{
  \item{obj}{
    An object created by the \code{pwrFDR} function with option
    \code{method}=="simulation".
  }
}
\value{
  A list with components
  \item{reps}{A data frame of \code{n.sim} rows containing the results of the
    simulations as columns: \code{M1}, the \# of non-null distributed
    statistics, \code{J}, the \# of statistical tests rejected under the BH-FDR
    procedure, and \code{S}, the \# of true positives.
  }
  \item{CCDF}{A data frame of dimension \code{n.T}, the number of test
    statistics with non-zero mean, by 2. The first column is \code{s}, 
    taking values 1 to \code{n.T}, and the second column is the monte 
    carlo estimated complementary CDF for the distribution of the 
    number of true positives, \code{S}}
  \item{X}{A single simulation replicate of the \code{m} raw test
    statistics, included as a sanity check}
}
\author{
  Grant Izmirlian izmirlig at mail dot nih dot gov
}

