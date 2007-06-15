\name{taxon.env}
\alias{taxon.env}

\title{Estimate taxon-environment relationships}

\description{
  \code{taxon.env} estimates taxon-environment relationships using
  logistic regression.  
}

\usage{
taxon.env(form, bcnt, envdata, bcnt.siteid, bcnt.abndid,
env.siteid,tlevs = "all", dumpdata = FALSE)
}

\arguments{
  \item{form}{Regression
    formula for model.
    Format is identical to that used for \code{glm}. }
  \item{bcnt}{A benthic count file with a full taxonomic hierarchy  that
    results from running \code{\link{get.taxonomic}}.}
  \item{envdata}{Environmental data file}
  \item{bcnt.siteid}{Field name for \code{bcnt} unique site identifier.}
  \item{bcnt.abndid}{Field name for \code{bcnt} taxon abundance.}
  \item{env.siteid}{Field name for \code{envdata} unique site
    identifier.}
  \item{tlevs}{Taxonomic levels to compute taxon-environment
    relationships.}
  \item{dumpdata}{Logical field specifying whether data used to estimate
    taxon-environment relationships should be included in the output.}
}

\details{
  This script generates logistic regression models for taxa at various levels of
  taxonomic resolution (as specified by the user).  The script
  computes the degrees of freedom specified by the model ($df$) and
  selects taxa that occur in at least $10 \times df$ and is absent from at
  least $10 \times df$ sites.  Then, logistic regression models based on
  the formula specified by the user are
  computed for each taxon using \code{\link{glm}}.
  The user can select the taxonomic levels (e.g., CLASS, ORDER, SPECIES,
  etc.) at which taxon-environment relationships are computed using the
  parameter \code{tlevs}.  The default value for \code{tlevs} computes
  relationships at all levels.
  }
\value{
  Regression coefficient file that can be used with \code{mlsolve} to
  infer environmental conditions.
  \item{tnames}{Taxon names for which models are calculated}
  \item{csave}{Matrix of regression coefficients}
  \item{xvar}{Character vector with names of explanatory variables}
  \item{xlims}{Range of each explanatory variable}
  \item{form}{Model formula}
  \item{roc}{Area under receiver operator curve for model for each
    taxon.}
  \item{raw.data}{Data used to estimate relationships.  Only included if
    \code{dumpdata = T}.}
  
}
\examples{
data(envdata.OR)
data(bcnt.tax.OR)
coef <- taxon.env(~sed + sed^2, bcnt.tax.OR, envdata.OR, "SVN",
"CountValue", "STRM.ID")
}
\seealso{\code{\link{get.taxonomic}}}
\author{
Lester L. Yuan
}
\keyword{models}
