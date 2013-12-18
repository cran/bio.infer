\name{tklist.modal}
\alias{tklist.modal}

\title{Modal select list dialog using TclTk}

\description{
  \code{tklist.modal} generates a dialog box that prompts the user to
  select from a list of options.  This script uses TclTk tools that
  should be operable on different operating systems.
}

\usage{
tklist.modal(title, elements0, returnValOnCancel="ID_CANCEL",selectmode =
"single")
}

\arguments{

  \item{title}{Title for the list box.}
  \item{elements0}{Character vector of elements listing different
    choices}
  \item{returnValOnCancel}{Character string returned by the script if
    the dialog box is canceled}
  \item{selectmode}{Set to multiple to allow for multiple selections.}

}

\value{
  A TclTk dialog box.
}
\examples{
\dontrun{
tklist.modal("Select a color", c("red", "green", "blue"))
}
}
\author{
Lester L. Yuan
}
\keyword{models}
