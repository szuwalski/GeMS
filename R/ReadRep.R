#' Pull parameters from .REP file. For use in \code{\link{AgeStructureComp}}
#' 
#' @param wd Directory containing simass.rep file from running age-structured model.
#' 
#' @return Named list of quantities from .REP file
#'
#' @export
#' @examples
#' \dontrun{
#' WorkingDir <- "~/GeMS/inst/exdata/Cod_5_AgeStructure/Cod_AgeStructure_CTL/1/70"
#' ReadPar(WorkingDir)
#' }
#' 
ReadRep <- function(wd) {
  datfile <- readLines(file.path(wd,"simass.dat"))
  yearsdat <- as.numeric(datfile[grep("survey years", datfile)[1]+1])
  rawfile <- readLines(file.path(wd,"simass.rep"))
  patterns <- c("^[[:alpha:]]","#")
  namedlines <- grep(paste(patterns,collapse="|"),rawfile)
  retlist <- list()
  
  listnames <- gsub("#","",rawfile[namedlines])
  listnames <- gsub(" ", "_", listnames)
  
  patterns <- c("pred", "bio")
  toshorten <- Reduce(`&`, lapply(patterns, grepl, listnames))
  temp <- grep("biomass",listnames)
  toshorten[temp] <- TRUE
  
  for(i in 1:(length(namedlines)-1)) {
    cnt <- namedlines[i]+1
    retlist[[i]] <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt],split=" "))))
    if(length(retlist[[i]])>1) retlist[[i]] <- retlist[[i]][-1]
    if(toshorten[i]) retlist[[i]] <- retlist[[i]][1:yearsdat]
    while(cnt<(namedlines[i+1]-1)) {
      temp <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt+1],split=" "))))[-1]
      if(toshorten[i]) temp <- temp[1:yearsdat]
      retlist[[i]] <- rbind(retlist[[i]],temp)
      cnt<- cnt+1
    }
    rownames(retlist[[i]]) <- c()
  }
  
  cnt <- max(namedlines)+1
  fin <- length(namedlines)
  retlist[[fin]] <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt],split=" "))))
  if(length(retlist[[fin]])>1) retlist[[fin]] <- retlist[[fin]][-1]
  if(toshorten[fin]) retlist[[fin]] <- retlist[[fin]][1:yearsdat]
  while(cnt<length(rawfile)) {
    temp <- suppressWarnings(as.numeric(unlist(strsplit(rawfile[cnt+1],split=" "))))
    if(toshorten[i]) temp <- temp[1:yearsdat]
    retlist[[fin]] <- rbind(retlist[[fin]],temp)
    cnt <- cnt+1
  }
  
  rownames(retlist[[length(namedlines)]]) <- c()
  
  names(retlist) <- listnames
  
  retlist$yearsDat <- yearsdat
  return(retlist)
}