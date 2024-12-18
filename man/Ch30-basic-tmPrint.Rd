\name{basic.tmPrint}
\alias{basic.tmPrint}
\title{
Wrapper to Print a Basic Nicely Formatted Table
}
\description{
Creates a generic call to print.TableMonster which in turn calls xtable
}
\usage{
basic.tmPrint(x, special = NULL, simple = FALSE, dbg = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    Any data.frame object. Here, the result of a call to pwrFDR.
}
  \item{special}{
    Special arguments to print.TableMonster. See package documentation.
}
  \item{simple}{
    The simplest use case
}
  \item{dbg}{
    Set to a value >= 1 for debugging
}
  \item{\dots}{
    Other arguments
}
}
\value{
  The value returned is an invisible version of the argument `x'.
}
\author{
  Grant Izmirlian
}
