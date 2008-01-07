# 1.7.2009
# Select coefficient file from those available

sel.coeffile <- function() {

  path <- R.home()
  flist <- system(paste("ls ", path, "//library//bio.infer//data", sep = ""),
                  intern = TRUE)
  incvec <- substring(flist, 1, 4) == "coef"

  flist <- flist[incvec]

  fsel <- tk_select.list(flist, title = "Available coefficient files")

  return(fsel)
}

