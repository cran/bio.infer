\name{get.dupe.sel}
\alias{get.dupe.sel}
\title{Get duplicate selection}
\description{
Prompt user to select correct version from duplicated entries
}
\usage{
get.dupe.sel(sumstr)
}
\arguments{
  \item{sumstr}{Character vector summarizing duplicated taxon names}
}
\details{
  Prompts the user via a tcltk window to select appropriate name from
  duplicate taxon names.
}
\value{
  Returns integer indicating which names are selected.
}
\author{Lester L. Yuan}
\keyword{models}
