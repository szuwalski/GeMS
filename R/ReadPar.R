#' Pull parameters from .PAR file. For use in \code{\link{AgeStructureComp}}
#' 
#' @param wd Directory containing .PAR file from running age-structured model.
#' 
#' @return Named list of quantities from .PAR file
#'
#' @export
#' @examples
#' \dontrun{
#' WorkingDir <- "~/GeMS/inst/exdata/Cod_5_AgeStructure/Cod_AgeStructure_CTL/1/70"
#' ReadPar(WorkingDir)
#' }
#' 

ReadPar <- function(wd) {
  rawfile <- readLines(file.path(wd,"simass.par"))
  namedlines <- grep("#",rawfile)
  retlist <- list()
  temp <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[1],split=" "))))
  temp <- temp[!is.na(temp)]
  retlist[[1]] <- temp[1]
  retlist[[2]] <- temp[2]
  retlist[[3]] <- temp[3]
  namedlines <- namedlines[-1]
  
  listnames <- unlist(strsplit(rawfile[namedlines],split="# "))
  listnames <- unlist(strsplit(listnames,split=":"))
  
  for(i in 1:length(listnames)) {
    retlist[[3+i]] <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[namedlines[i]+1],split=" "))))
    if(length(retlist[[3+i]])>1) retlist[[3+i]] <- retlist[[3+i]][-1]
  }
  
  listnames <- c("Npars", "objfun", "grad", listnames)
  names(retlist) <- listnames
  
  return(retlist)
}