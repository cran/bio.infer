\name{get.taxonomic}
\alias{get.taxonomic}
\title{Merge benthic count data with the ITIS taxonomy table}

\description{
  \code{get.taxonomic} merges the taxon names provided in a benthic
  count data file with a standardized taxonomy table.  The
  purpose of this function is to provide the complete taxonomic
  hierarchy for each taxon.
} 

\usage{
get.taxonomic(bcnt)
}
\arguments{
  \item{bcnt }{A benthic count data frame with the following three
    fields: a sample identifier, taxon name, and abundance.  The three
    fields must be ordered as listed.  Sample identifier is typically a
    character string that is a unique identifier for a given
    sample. Taxon name is a character string that specifies the name of
    the taxon observed in the sample, and abundance is the number of
    individuals of that taxon that were observed.
  }
}

\details{
  \code{get.taxonomic} attempts to match taxon names provided in a
  benthic count data file with taxa listed in ITIS.  The script
  automatically corrects for the most common conventions in taxonomic
  naming.  For example, taxon names are converted to all capital
  letters, and trailing strings such as "SP." and "SPP." are deleted.
  
  
  The script lists taxa remaining that are not
  matched to entries in the ITIS database using \code{\link{fix}}.  The
  user can then enter a 
  corrected spelling for the taxon name.  
  
  The script provides a final summary of the unmatched taxon names, and
  the user can continue revising names, or complete the program.
}

\value{
  A final taxonomic table is returned as tab-delimited text in the file
  \code{sum.tax.table.txt}.  This file is best viewed using a
  spreadsheet.   Within 
  R, the script returns a benthic count file with the original taxon
  name associated with a full taxonomic hierarchy for that taxon.
  The assigned species names can be modified by hand in
  \code{sum.tax.table.txt} and the revised information loaded with
  \code{\link{load.revised.species}}. 
}
\examples{
data(bcnt.OR)
data(itis.ttable)
\dontrun{
bcnt.tax.OR <- get.taxonomic(bcnt.OR)
}
}
\seealso{\code{\link{itis.ttable}}}
\author{
Lester L. Yuan
}
\keyword{models}
