\name{envdata.OR}
\alias{envdata.OR}
\docType{data}
\title{Environmental data from western OR}
\description{
  Environmental data for stream sampled in western Oregon.
}
\usage{data(envdata.OR)}
\format{
  A data frame with 245 observations on the following 13 variables.
  \describe{
    \item{\code{STRM.ID}}{Unique sample code}
    \item{\code{jday}}{Sampling day}
    \item{\code{lon}}{Longitude in decimal degrees}
    \item{\code{lat}}{Latitude in decimal degrees}
    \item{\code{sed}}{Percent sands and fines}
    \item{\code{elev}}{Log-transformed elevation}
    \item{\code{slope}}{Slope}
    \item{\code{area}}{Log-transformed catchment area}
    \item{\code{sed.log}}{Log-transformed percent sands and fines}
    \item{\code{temp}}{Grab stream temperature}
    \item{\code{temp.avg}}{Seven day average maximum temperature}
    \item{\code{elev.ut}}{Elevation}
    \item{\code{sp}}{Estimate of unit stream power}
  }
}
\source{
  Oregon Department of Environmental Quality
}
\examples{
data(envdata.OR)
}
\keyword{datasets}
