\name{bio.infer-package}
\alias{bio.infer-package}
\alias{bio.infer}
\docType{package}
\title{
Scripts for computing biological inferences
}
\description{
Reformats raw benthic count files, combines them with existing
taxon-environment relationships, and computes inferences of
environmental conditions.
}
\details{
\tabular{ll}{
Package: \tab bio.infer\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-10-31\cr
License: \tab GPL version 2 or newer\cr
}
Use this package to import benthic count data, reformat the data, and
compute inferences of environmental conditions using taxon-environment
relationships.

Typically, the user will start by loading a local benthic count
data file (see  \code{\link{bcnt.OR}}
for an example of the required format).
This benthic count file will be merged with standardized
taxonomic hierarchy information using \code{\link{get.taxonomic}}.
Then, operational taxonomic units will be assigned, based on the number
of occurrences of different taxa and based on the identity of taxa for
which environmental preferences have been quantified. (See
\code{\link{get.otu}} for the OTU assignment script, and
\code{\link{coef.west.wt}} for an example file of environmental
preferences.)

Once OTUs have been assigned, the benthic count file is reformatted as a
site-taxon matrix (\code{\link{makess}}).  This file can then be used to
compute biological inferences (\code{\link{mlsolve}}).

This package was developed in Windows but should be fully functional in other
operating systems.

}
\author{
Lester L. Yuan
}
\keyword{models}
