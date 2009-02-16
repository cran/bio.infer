\name{resolve.mult}
\alias{resolve.mult}
\title{Resolve cases in which two taxa are found in ITIS}
\description{
  Resolves cases in which two taxa listed on the same entry are both
  found in ITIS. 
}
\usage{
resolve.mult(parse.list, get.tax.env)
}
\arguments{
  \item{parse.list}{List of taxa that are found in ITIS (position 1)and
  taxa that are still unresolved (position 2).}
  \item{get.tax.env}{Environment variable where the full ITIS taxon
    table is available}
}
\details{
  The most likely reason for two valid taxa being listed in the same
  entry is that the taxonomist could only identify a specimen down to
  two possibilities.  In these cases, the script downgrades the
  identification to a coarser level that is the same for both taxa.
}
\value{
  Updates parse.list so that resolved multiple listings are included in
  position 1.
}
\author{Lester L. Yuan}
\keyword{models}
