\name{locate.dupes}
\alias{locate.dupes}
\title{Locate duplicated taxon names}
\description{
  Locate taxon names that match multiple records in ITIS.
}
\usage{
locate.dupes(fulltab)
}
\arguments{
  \item{fulltab}{Full taxonomic table for matched taxon names}
}
\details{
  Finds duplicate names in fulltab
}
\value{
  Returns string that summarizes duplicates and their row numbers.
}
\author{Lester L. Yuan}
\keyword{models}
