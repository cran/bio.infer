\name{remove.dupes}
\alias{remove.dupes}
\title{Remove unselected duplicates from taxon table}
\description{
Remove unselected duplicate taxon names from taxon table.
}
\usage{
remove.dupes(fulltab, dupe.list, dupe.sel)
}
\arguments{
  \item{fulltab}{Full taxon table for matched taxon names.}
  \item{dupe.list}{List of duplicated entries.}
  \item{dupe.sel}{Character vector of taxa selected by user.}
}
\details{
  Deletes rows from fulltab that are not selected by user
}
\value{
  Returns revised fulltab.
}
\author{Lester L. Yuan}
\keyword{models}
