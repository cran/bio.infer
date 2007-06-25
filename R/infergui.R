"infergui" <-
function() {

  getfile <- function() {
    name <- tclvalue(tkgetOpenFile())
    if (name == "") return;

    bcnt <<- read.delim(name)
    
    if (ncol(bcnt) != 3) { 
      tkmessageBox(message="Please make sure your benthic count file is tab-delimited and has the following three fields: site ID, taxon name, taxon abundance",icon="error",type="ok")
      bcnt <<- NULL
    }
    
  }

  getfile.r <- function() {
    getfile()
    if (! is.null(bcnt)) {
      tkinsert(txt, "end", "Biological data loaded.\n")
      statevec <- c("normal", "disabled" ,"disabled")
      tkdestroy(tt)
      draw.toplev(statevec)
    }
  }
      
  select.te <- function() {
    data(coef.west.wt)
    fnames <- ls(name = globalenv(), pattern = "coef")
    tefname <<- tklist.modal("Select coefficient file:", fnames)
    if (length(tefname) > 0) {
      tkinsert(txt, "end", "Taxon-environment coefficient file selected.\n")
      tkdestroy(tt)
      draw.toplev(c("normal", "normal", "disabled"))
    }
    else {
      tkmessageBox(message = "Select coefficient file to continue",
                   icon = "error", type = "ok")
    }
  }

  compute.inf <- function() {
    tkconfigure(tt, cursor = "wait")
    coef <- get(tefname)
    data(itis.ttable)
    tkinsert(txt, "end", "Merging biological data with ITIS...\n")
    tkfocus(txt)
    tkfocus(tt)
    bcnt.tax <- get.taxonomic(bcnt, itis.ttable, gui = TRUE)
    tkinsert(txt, "end", "Assigning operational taxonomic units...\n")
    tkgrab.release(txt)
    flush.console()
    bcnt.otu <- get.otu(bcnt.tax, coef, gui = TRUE)
    tkinsert(txt, "end", "Computing inferences...\n")
    flush.console()
    ss <- makess(bcnt.otu)
    inf.out <<- mlsolve(ss, coef)
    tkinsert(txt, "end", "Inferences computed!")
    tkconfigure(tt, cursor = "arrow")
    flush.console()
    if (! is.null(inf.out)) {
      tkdestroy(tt)
      draw.toplev(c("normal", "normal", "normal"))
    }
  }

  export.res <- function() {
    name <- tclvalue(tkgetSaveFile())
    if (name != "") {
      write.table(inf.out, file = name, sep = "\t",
                  row.names = FALSE)
    }
  }

  quitinf <- function() {
    tkgrab.release(tt)
    tkdestroy(tt)
    tkdestroy(msgwin)
  }
  
  draw.toplev <- function(statevec) {
    tt <<- tktoplevel()

    tkwm.title(tt, "Biological inferences")
    button.widget1 <- tkbutton(tt, text = "Load biological data",
                               command = getfile.r)
    
    button.widget2 <- tkbutton(tt,
                               text = "Select taxon-environment relationships",
                               command = select.te, state = statevec[1])
    
    button.widget3 <- tkbutton(tt, text = "Compute inferences",
                               command = compute.inf, state = statevec[2])
    
    button.widget4 <- tkbutton(tt, text = "Export inferences",
                               command = export.res, state = statevec[3])
    button.widget5 <- tkbutton(tt, text= "Quit", command = quitinf)
    
    tkgrid(tklabel(tt, text = " "))
    tkgrid(tklabel(tt, text = "  "), button.widget1, tklabel(tt, text = "  "))
    tkgrid(tklabel(tt, text = " "))
    tkgrid(tklabel(tt, text = "  "),button.widget2, tklabel(tt, text = "  "))
    tkgrid(tklabel(tt, text = " "))
    tkgrid(tklabel(tt, text = "  "),button.widget3, tklabel(tt, text = "  "))
    tkgrid(tklabel(tt, text = " "))
    tkgrid(tklabel(tt, text = "  "),button.widget4, button.widget5,
           tklabel(tt, text = "  "))
    tkgrid(tklabel(tt, text = " "))
  }

  bcnt <<- NULL
  tefname <<- character(0)
  statevec <- c("disabled", "disabled", "disabled")
  # Set up status window
  msgwin <<- tktoplevel()
  scr <- tkscrollbar(msgwin, repeatinterval=5,
                       command=function(...)tkyview(txt,...))
  tkwm.title(msgwin, "Messages")
  txt <<- tktext(msgwin,yscrollcommand=function(...)tkset(scr,...))
  tkgrid(txt, scr)
  tkgrid.configure(scr, sticky = "ns")
  tkconfigure(txt, font = "courier")
  tkinsert(txt, "end", "*** Predicting environmental conditions from biological observations ***\n")

  draw.toplev(statevec)


}


