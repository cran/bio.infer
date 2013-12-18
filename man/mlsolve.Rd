\name{mlsolve}
\alias{mlsolve}
\title{Maximum likelihood inferences}

\description{\code{mlsolve} uses maximum likelihood estimation to infer
  environmental conditions from biological observations.}

\usage{
mlsolve(ss, coef0, site.sel = "all", bruteforce = FALSE)
}
\arguments{
  \item{ss}{A site-species matrix, with sample identifiers in the first
    column, and abundances.  Typically, the result of \code{\link{makess}}.
  }
  \item{site.sel}{Character vector of sample id's at which inferences
    should be computed.  Default selection of 'all' forces script to
    compute inferences for all samples in the matrix.}
  \item{coef0}{Coefficient file with model specifications and regression
    coefficients defining 
    taxon-environment relationships. See \code{\link{taxon.env}} for
    more details.}
  \item{bruteforce}{Logical flag.  Select \code{TRUE} to compute solution by
    brute force rather than by conjugate gradients.  Provides a means of
    examining likelihood contours when needed.}
}

\details{
  \code{mlsolve} uses maximum likelihood estimation to infer environmental
  conditions from biological observations.  Maximum likelihood
  estimation in this context is a constrained optimization problem, in
  which we wish find the point at which the likelihood function is
  maximized, constrained by the range of the environmental variables in
  the calibration data.  \code{mlsolve} formulates the likelihood
  function and calls \code{\link{optim}} to solve the optimization
  problem.
  
  The option \code{bruteforce} can be used to check the performance of the
  iterative solver.  When \code{bruteforce} is set to be \code{TRUE},
  the script also 
  computes likelihood values for a uniformly-spaced grid that spans the
  ranges of the environmental variables.  This is a time-consuming
  calculation, and therefore should only be attempted on a few samples.
  Once the grid is computed, though, contours of the likelihood surface
  are plotted.

  \code{mlsolve} requires that the taxa included in \code{ss} each have
  associated taxon-environment information in \code{coef0}.  The best
  way to ensure that this condition is satisfied is to make sure that
  \code{\link{get.otu}} is run with the same \code{coef0} as used in
  \code{mlsolve}, and that the results from \code{get.otu} are used to
  generate \code{ss} (see \code{\link{makess}}). 
  }

  \examples{
data(coef.west.wt)
data(ss.OR)
inf.out <- mlsolve(ss.OR, coef.west.wt)
}
\value{ 
The script returns a dataframe containing inferences for each sample as
well as a logical flag indicating whether solutions computing with four
different starting locations are consistent with one another.
}
\author{
Lester L. Yuan
}
\keyword{models}
