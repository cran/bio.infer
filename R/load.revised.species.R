"load.revised.species" <-
  function(bcnt.tax, fname) {

    tname <- names(bcnt.tax)[2]

    dftemp <- read.delim(fname)
    
    imatch <- match("SPECIES", names(bcnt.tax))
    
    names0 <- names(bcnt.tax)
    
    bcnt.tax <- bcnt.tax[, -imatch]
    
    dftemp <- dftemp[, c(tname, "SPECIES")]

    bcnt.tax <- merge(bcnt.tax, dftemp, by = tname)
    bcnt.tax$SPECIES <- levels(bcnt.tax$SPECIES)[bcnt.tax$SPECIES]
    bcnt.tax <- bcnt.tax[, names0]

    return(bcnt.tax)
  }



