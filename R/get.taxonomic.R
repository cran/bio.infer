"get.taxonomic" <-
function(bcnt, itis.ttable = NULL, exlocal = character(0),
         outputFile = "sum.tax.table.txt", gui = FALSE, txt = NULL) {

  if (is.null(itis.ttable)) {
    data(itis.ttable, envir=environment(NULL))
  }
  names0 <- names(bcnt)
  f.tname <- names0[2]

  tlevs <- names(itis.ttable)
  imatch <- match("TAXON", tlevs)
  tlevs <- tlevs[-imatch]
  

  if (is.factor(bcnt[,2])) {
    f1 <- sort(unique(levels(bcnt[,2])[bcnt[,2]]))
  }
  else {
    if (is.character(bcnt[,2])) {
      f1 <- sort(unique(bcnt[,2]))
    }
    else {
        tkmessageBox(message = "2nd field is neither factor nor character",
                     icon = "error", type = "ok")
    }
  }

  f2 <- toupper(f1)
  dfref <- data.frame(I(f1), I(f2))

  itis.taxa <- itis.ttable$TAXON
  
#  df1 <- merge(dfref, itis.ttable, by.x = "f2", by.y = "TAXON", all.x = TRUE)
#  df1$nomatch <- is.na(df1[, tlevs[1]])
#  dfref <- unique.data.frame(df1[, c("f1", "f2", "nomatch")])

  substr <- as.list(rep(NA, times = 1))
  i <- 1
  tmiss0 <- dfref$f2

  # get rid of anything inside parentheses
  w1 <- regexpr("\\(", tmiss0)
  w2 <- regexpr("\\)", tmiss0)
  incvec <- (w1 != -1) & (w2 != -1)
  tmiss0[incvec] <- paste(substring(tmiss0[incvec], 1, w1[incvec]-1),
                          substring(tmiss0[incvec], w2[incvec]+1,
                                    nchar(tmiss0[incvec])))

  # parse all continuous character strings into separate storage lists
  repeat {
    w <- regexpr("[A-Z]+", tmiss0)
    if (sum(w != -1) == 0) break
    substr[[i]] <- substring(tmiss0, w, w+attributes(w)$match.length-1)
    w3 <- w + attributes(w)$match.length
    tmiss0 <- substring(tmiss0, w3, nchar(tmiss0))
    if (sum(tmiss0 != "") == 0) break
    i <- i + 1
  }

  # First guess at species names
  exlist <- c("DUPLICATE", "SETAE", "CODE", "GROUP", "TYPE", "GENUS",
              "PANEL", "SAND", "TURRET", "CASE", "LARVAE", "ADULT",
              toupper(exlocal))
  sp.name <- rep("", times = length(substr[[1]]))
  if (length(substr) > 1) {
    for (i in 1:length(substr[[1]])) {
      for (j in 2:length(substr)) {
        if ((nchar(substr[[j]][i]) > 3) & (! (substr[[j]][i] %in% exlist))) {
          if (sp.name[i] == "") {
            sp.name[i] <- substr[[j]][i]
          }
          else{
            sp.name[i] <- paste(sp.name[i], substr[[j]][i], sep = "/")
          }
        }
      }
    }
  }

  dfref$f2 <- substr[[1]]
  dfref$sp.name <- sp.name

  #Check for any compound taxa
  #Only allow compound taxa up to the same family
  if (length(substr) > 1) {
    imatch <- match("SUBCLASS", toupper(tlevs))
    tlevs.loc <- rev(tlevs[length(tlevs):imatch])
    for (i in 1:nrow(dfref)) {
      if (nchar(substr[[2]][i]) > 3) {
        imatch1 <- match(substr[[1]][i], itis.taxa)
        imatch2 <- match(substr[[2]][i], itis.taxa)
        if (! is.na(imatch1) & ! is.na(imatch2)) {
          comp1 <- itis.ttable[imatch1, tlevs.loc]
          comp2 <- itis.ttable[imatch2, tlevs.loc]
        tlev.sav <- ""
          for (j in 1:length(comp1)){
            if (! is.na(comp1[j]) & ! is.na(comp2[j])) {
              if ((comp1[j] == comp2[j]) & (comp1[j] != "")) {
                tlev.sav <- tlevs.loc[j]
              }
            }
          }
          if (tlev.sav != "") {
            dfref$f2[i] <- comp1[,tlev.sav]
          dfref$sp.name[i] <- ""
          }
        # default here if both strings match to something
        # but a common higher level is not found is genus.species
        }
        else {
          if (!is.na(imatch1)) {
          # check whether first string is a genus
            if (is.na(itis.ttable[imatch1, "GENUS"])) {
              dfref$sp.name[i] <- ""
            }
            else {
              if (itis.ttable[imatch1, "GENUS"] == "") {
                dfref$sp.name[i] <- ""
              }
            }
          }
          else {
            if (! is.na(imatch2)) {
              dfref$f2[i] <- substr[[2]][i]
              dfref$sp.name[i] <- ""
            }
          }
        }
      }
    }
  }

  tmiss0 <- character(0)
  for (i in 1:length(dfref$f2)) {
    if (! (dfref$f2[i] %in% itis.taxa)) {
      tmiss0 <- c(tmiss0, dfref$f2[i])
    }
  }
  
  tmiss0 <- sort(unique(tmiss0))

  if (length(tmiss0) > 0) {

    tkmessageBox(message= "Please correct unrecognized taxon names.",
                 icon = "info", type = "ok" )

    repeat {
      
      repeat {
        cstr <- rep("", times = length(tmiss0))
        dfcorrect <- data.frame(I(tmiss0), I(cstr))
        names(dfcorrect) <- c("TAXANAME", "CORRECTION")
#        dfcorrect <- edit(dfcorrect)
        dfcorrect$CORRECTION <- modalDialog("Correct misspellings",
                                            tmiss0, "")
        changetname <- FALSE
        for (i in 1:nrow(dfcorrect)) {
          changetname <- tmiss0[i] != dfcorrect$TAXANAME[i]
          if (changetname) {
            tkmessageBox(message = "Do not change TAXANAME\n",
                         icon = "error", type = "ok")
            flush.console()
            break
          }
        }
        if (!changetname) break
      }
      
      dfcorrect$CORRECTION <- toupper(dfcorrect$CORRECTION)
      if (gui) {
        tkinsert(txt, "end", "\nCorrections entered:\n")
        tkinsert(txt, "end", "Taxon name")
        tkinsert(txt, "end", "\t")
        tkinsert(txt, "end", "Corrected name")
        tkinsert(txt, "end", "\n")
        for (iii in 1:nrow(dfcorrect)) {
          tkinsert(txt, "end", dfcorrect[iii,1])
          tkinsert(txt, "end", "\t")
          tkinsert(txt, "end", dfcorrect[iii,2])
          tkinsert(txt, "end", "\n")
        }
        tkinsert(txt, "end", "\n")
        tkfocus(txt)
      }
      else {
        cat("Corrections entered.\n")
        print(dfcorrect)
        cat("\n")
      }
      
      for (i in 1:nrow(dfcorrect)) {
        incvec <- dfref$f2 == dfcorrect$TAXANAME[i]
        if (! is.na(dfcorrect$CORRECTION[i])) {
          imatch <- match(toupper(dfcorrect$CORRECTION[i]), itis.taxa)
          
          if (is.na(imatch)) {
            if (dfcorrect$CORRECTION[i] != "") {
              if (gui) {
                tkmessageBox(message = paste(dfcorrect$CORRECTION[i], "is not found in ITIS"), icon = "error", type = "ok")
              }
              else {
                cat(dfcorrect$CORRECTION[i], " is not found in ITIS.\n")
              }
            }
          }
          else {
            dfref$f2[incvec] <- dfcorrect$CORRECTION[i]
            specpres <- dfref$sp.name != ""
            if (sum(specpres & incvec) > 0) {
              isel <- specpres & incvec
              if (itis.ttable$GENUS[imatch] == "") {
                dfref$sp.name[isel] <- ""
              }
            }
          }
        }
      }
      
      df1 <- merge(dfref, itis.ttable, by.x = "f2", by.y = "TAXON", all.x = TRUE)
      nomatch <- is.na(df1[, tlevs[1]])
      
      df2 <- merge(bcnt, dfref, by.x = f.tname, by.y = "f1", all.x = TRUE)
      df2 <- merge(df2, itis.ttable, by.x = "f2", by.y = "TAXON", all.x = TRUE)
      
      incvec <- is.na(df2[, tlevs[1]])
      if (sum(incvec) > 0) {
        tmiss <- sort(unique(df2[incvec, "f2"]))
        sumnumocc <- rep(NA, times = length(tmiss))
        for (i in 1:length(tmiss)) {
          sumnumocc[i] <- sum(df2[, "f2"] == tmiss[i])
        }
        dfsumnumocc <- data.frame(I(tmiss), sumnumocc)
        names(dfsumnumocc) <- c("TAXANAME", "NUMBER OF OCCURRENCES")

        if (gui) {
          tkinsert(txt, "end", "\nSummary of taxa without matches.\n")
          tkinsert(txt, "end", "Taxon name")
          tkinsert(txt, "end", "\t")
          tkinsert(txt, "end", "Number of occurrences")
          tkinsert(txt, "end", "\n")
          for (iii in 1:nrow(dfsumnumocc)) {
            tkinsert(txt, "end", dfsumnumocc[iii,1])
            tkinsert(txt, "end", "\t")
            tkinsert(txt, "end", dfsumnumocc[iii,2])
            tkinsert(txt, "end", "\n")
          }
          
          tkinsert(txt, "end", "\n")
          tkfocus(txt)
        }
        else {
          cat("Summary of taxa without matches: \n")
          format(dfsumnumocc, justify = "centre")
          print(dfsumnumocc)
        }
      }
      flush.console()
      
      done0 <- tkmessageBox(message="Are you done editing?",
                            icon = "question", type = "yesno",
                            default = "yes")
      if (as.character(done0) == "yes") break
    }
  }

  # Final check for duplications in itis table
  dftemp <- data.frame(I(unique(dfref$f2)))
  names(dftemp) <- c("f2")
  dftemp2 <- merge(itis.ttable, dftemp, by.x = "TAXON", by.y = "f2",
                   all.y = TRUE)
  reps <- unique(dftemp2$TAXON[duplicated(dftemp2$TAXON)])
  if (length(reps) > 0) {
    tkmessageBox(message = "The following taxa match with multiple ITIS records.Please select the appropriate taxon from list", icon = "info", type = "ok")
    for (i in 1:length(reps)) {
      isav <- (1:nrow(dftemp2))[reps[i] == dftemp2$TAXON]
      print(dftemp2[isav, c("TAXON", "PHYLUM", "CLASS", "ORDER", "FAMILY")])
      sumstr <- rep("", times = length(isav))
      for (j in 1:length(isav)) {
        sumstr[j] <- paste(dftemp2[isav[j],
                                   c("PHYLUM", "CLASS", "ORDER", "FAMILY")],
                           collapse = "-")
      }
      a <- tklist.modal(paste("Select appropriate taxon for ", reps[i]),
                              sumstr)
      isel <- match(a, sumstr)
      isav <- isav[-isel]
      dftemp2 <- dftemp2[-isav,]
    }
  }
  itis.ttable.loc <- dftemp2[, names(itis.ttable)]

  # Eliminate fields with no entries
  iomit <- numeric(0)
  for (i in 1:length(tlevs)) {
    if (sum(itis.ttable.loc[, tlevs[i]] != "", na.rm = TRUE) == 0) {
      iomit <- c(iomit,i)
    }
  }
  if (length(iomit) > 0) {
    tlevs <- tlevs[-iomit]
  }

  # Only allow species names for valid genera
  df1 <- merge(dfref, itis.ttable.loc, by.x = "f2", by.y = "TAXON",
               all.x = TRUE)
  df1$SPECIES <- paste(df1$GENUS, df1$sp.name, sep = ".")
  incvec <- (nchar(df1$sp.name) > 0) & (df1$GENUS != "")
  incvec[is.na(incvec)] <- FALSE
  df1$SPECIES[! incvec] <- ""
  dfref <- df1[, c("f1", "f2", "SPECIES")]

  # delete site id and count fields to generate summary taxatable
  df1 <- df1[, c(tlevs, "SPECIES", "f1")]
  df1 <- df1[do.call(order, df1),]
  names(df1) <- c(tlevs, "SPECIES", f.tname)

  if (is.character(outputFile)) {
    write.table(df1, sep = "\t", file = outputFile, row.names = FALSE)
    tkmessageBox(message = paste("Check final taxa name assignments in", outputFile), icon = "info", type = "ok")
  }
  df2 <- merge(bcnt, dfref, by.x = f.tname, by.y = "f1", all.x = TRUE)

  df2 <- merge(df2, itis.ttable.loc, by.x = "f2", by.y = "TAXON", all.x = TRUE)

  varlist <- c(tlevs, "SPECIES")

  for (i in 1:length(varlist)) {
    incvec <- df2[, varlist[i]] == ""
    incvec[is.na(incvec)] <- FALSE
    df2[incvec, varlist[i]] <- NA
  }

  return(df2[, c(names0, tlevs, "SPECIES")])
  
}

