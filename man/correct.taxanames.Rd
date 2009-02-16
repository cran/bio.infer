\name{correct.taxanames}
\alias{correct.taxanames}
\title{Correct unrecognized taxon names}
\description{
Prompts user to correct unrecognized taxon names
}
\usage{
correct.taxanames(tname.old, get.tax.env)
}
\arguments{
  \item{tname.old}{Character vector of unrecognized names}
  \item{get.tax.env}{Environment for ITIS taxon table}
}
\details{
  Prompts user to enter corrections to names via a tcltk dialog box.
}
\value{
  Character vector of corrected names
}
\author{Lester L. Yuan}
\keyword{models}
