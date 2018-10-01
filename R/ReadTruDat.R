#' Pull parameters from TrueQuantities.DAT file. For use in \code{\link{AgeStructureComp}}
#' 
#' @param wd Directory containing TrueQuantities.DAT file from running age-structured model.
#' 
#' @return Named list of quantities from TrueQuantities.DAT file
#'
#' @export
#' @examples
#' \dontrun{
#' WorkingDir <- "~/GeMS/inst/exdata/Cod_5_AgeStructure/Cod_AgeStructure_CTL/1/70"
#' ReadTruDat(WorkingDir)
#' }
#' 

ReadTruDat <- function(wd) {
  rawfile <- readLines(file.path(wd,"TrueQuantities.DAT"))
  namedlines <- grep("#",rawfile)[-c(1:2,9)]
  retlist <- list()
  
  listnames <- gsub("#",'',rawfile[namedlines])
  listnames <- stringr::str_trim(listnames)
  listnames <- gsub(" ","_",listnames)
  
  for(i in 1:(length(namedlines)-1)) {
    cnt <- namedlines[i]+1
    retlist[[i]] <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt],split=" "))))
    while(cnt<(namedlines[i+1]-1)) {
      temp <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt+1],split=" "))))
      retlist[[i]] <- rbind(retlist[[i]],temp)
      cnt<- cnt+1
    }
    rownames(retlist[[i]]) <- c()
  }
  
  cnt <- max(namedlines)+1
  fin <- length(namedlines)
  retlist[[fin]] <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt],split=" "))))
  while(cnt<length(rawfile)) {
    temp <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt+1],split=" "))))
    retlist[[fin]] <- rbind(retlist[[fin]],temp)
    cnt <- cnt+1
  }
  
  rownames(retlist[[length(namedlines)]]) <- c()

  names(retlist) <- listnames
  
  return(retlist)
}