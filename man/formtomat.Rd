\name{formtomat}
\alias{formtomat}
\title{Analyze formula string}
\description{
  Analyzes formula string that defines the regression models for
  taxon-environment relationships.  Returns a numerical representation
  of the formula and the derivative of the formula that can be used by
  \code{mlsolve}. 
}
\usage{
formtomat(a, xvar)
}
\arguments{
  \item{a}{Character string defining regression formula.}
  \item{xvar}{Character vector with variable names.}
}
\details{
  Internal function to \code{mlsolve}.  
}
\value{
  List of transformations.
}
\author{Lester L. Yuan}
\keyword{models}

