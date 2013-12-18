\name{make.species}
\alias{make.species}
\title{Generate species names}
\description{
  Generate species names from possibilities in the parsed taxa name table.
}
\usage{
make.species(df.parse, fulltab)
}
\arguments{
  \item{df.parse}{Matrix with valid taxon names in the second column}
  \item{fulltab}{Full taxonomic table for matched taxon names}
}
\details{
  Generates likely species names from matched taxon names that are
  identified as Genus level names.
}
\value{
  Returns full taxonomic table with SPECIES field appended.
}
\author{Lester L. Yuan}
\keyword{models}
