\name{get.otu}
\alias{get.otu}

\title{Defines operational taxonomic units (OTUs)}

\description{
  \code{get.otu} generates OTUs for benthic count data using a simple
  algorithm that is designed to maximize the amount of taxonomic
  information retained across the entire data set, while eliminating
  redundant taxonomic identifications. 
} 
\usage{
get.otu(bcnt, optlist = NULL, ndc = TRUE, outputFile = "sum.otu.txt")
}

\arguments{
  \item{bcnt }{benthic count data frame that includes the full taxonomic
    hierarchy for each taxon.  Typically, \code{bcnt} is the output from
    \code{\link{get.taxonomic}}. 
  }
  \item{optlist }{vector of taxon names for which trait or tolerance value
  is available.  Specify \code{optlist = NULL} for cases in which OTU are
  computed without regard for any additional taxon specific information,
  e.g. for developing a RIVPACS model.}
\item{ndc }{Logical flag indicating whether OTU should be computed that
  eliminate double-counting of taxa.}
\item{outputFile}{File name for summary output file.}
}

\details{
  \code{get.otu} makes decisions out OTU according to the following rule.
  The number samples in which a given taxonomic group (e.g. Baetidae) is
  observed, is compared to the number of samples in which members of the
  that same group are identified to a higher taxonomic level (e.g.,
  Baetis, Fallceon, Acentrella, etc.).  If more samples are observed
  with the coarser identification, then all identifications are
  downgraded to the coarse identification.  If more samples are observed
  at the finer level of identification, then those sample are retained,
  and all identifications at the coarser level are omitted. 
  Before deciding on OTU assignments, the script first reviews the
  species names that are listed in the benthic count file and matches
  them with the species names provided in \code{optlist}.   
}

\value{
  Final OTU decisions are returned as tab-delimited text file,
  \code{sum.otu.txt}, which is best viewed using a spreadsheet.  A new benthic count
  file is also returned in R, including the original benthic count
  information and a new field specifying the OTU designation for each
  taxon.  Manual revisions to the assigned OTUs can be recorded in
  \code{sum.otu.txt} and loaded with \code{\link{load.revised.otu}}
}
\examples{
data(bcnt.tax.OR)
data(coef.west.wt)
bcnt.otu.OR <- get.otu(bcnt.tax.OR, coef.west.wt)
}

\author{
  Lester L. Yuan
}
\seealso{\code{\link{get.taxonomic}}}
\keyword{models}

