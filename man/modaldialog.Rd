\name{modalDialog}
\alias{modalDialog}

\title{Modal dialog box using TclTK}

\description{
  \code{modalDialog} generates a dialog box prompting the user to
  enter a list of text before the process continues.  Based on scripts
  provided in http://bioinf.wehi.edu.au/~wettenhall/RTclTkExamples/
}

\usage{
modalDialog(title,itemlist,entryInit,entryWidth=20,
                        returnValOnCancel="ID_CANCEL")
}

\arguments{

  \item{title}{Title for the dialog box}
  \item{itemlist}{Character vector of names that the user must correct}
  \item{entryInit}{Default string to fill in the entries}
  \item{entryWidth}{Width of entry string}
  \item{returnValOnCancel}{Value for script to return if canceled}
}

\value{
  Character vector of user entries.
}

\author{
Lester L. Yuan
}
\keyword{models}
