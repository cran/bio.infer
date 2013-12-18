\name{load.revised.species}
\alias{load.revised.species}

\title{Loads a user-revised list of species names}

\description{
  If manual correction of the species names assigned by \code{get.taxonomic} is desired, these
  corrections should be made on the \code{sum.tax.table.txt} file.
 The corrections should be saved as
  tab-delimited text, and then \code{load.revised.species} run.   
}

\usage{
load.revised.species(bcnt.tax, fname)
}
\arguments{
  \item{bcnt.tax}{benthic count data frame that includes the full
    taxonomic hierarchy for each taxon.}
  \item{fname}{Character string containing file name for the revised,
    tab-delimited, file in the format of \code{sum.tax.table.txt}.  See \code{\link{get.taxonomic}}.}
  }

\details{
  The script prompts the user to enter the name of the text file
  containing the revised species names.  The original "SPECIES" column
  in bcnt.tax is deleted and the new, revised version is merged into
  bcnt.tax. 
}

\value{
  The script returns a new benthic count file that incorporates the
  corrections to the species names
}
\seealso{\code{\link{get.taxonomic}}}
\author{
Lester L. Yuan
}
\keyword{models}
