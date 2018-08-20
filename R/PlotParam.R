#' Plot Life History Parameter
#' 
#' Plots life history parameter estimates from final assessment year
#' 
#' @param param Parameter of interest, named from PAR or CTL file (one of c("NatM","GrowthK","Linf","SelPars50","SelPars95","NatMn","VonKn","LinfN","sel50n","sel95n"))
#' @param OMfile Name of CTL file to base truth on
#' @param EMvec Vector of EM names
#' @param saveplot Save plot to TIFF
#' @param dolegend Logical; include legend?
#' @param ylims vector of length two containing limits for y-axis
#' @param ... any other graphical parameters, e.g., ylab,xlab
#' 
#' @return Plots plots plots; list of truth and estimates pulled from PAR and CTL files.
#'
#' @export
#'
PlotParam <- function(param,OMfile,EMvec,saveplot,plotnames=NA,dolegend=F,ylims=NA,...) {
  if (!param %in% c("NatM","GrowthK","Linf","SelPars50","SelPars95",
                   "NatMn","VonKn","LinfN","sel50n","sel95n"))
      {stop("param not listed in PAR or CTL file. It should be one of these:\n
                   NatM,GrowthK,Linf,SelPars50,SelPars95,
                   NatMn,VonKn,LinfN,sel50n,sel95n)")}
  
  if(param == "NatM" | param == "NatMn") {
    paramCTL <- "NatMn"
    paramPAR <- "NatM"
  }
  
  if(param == "GrowthK" | param == "VonKn") {
    paramCTL <- "VonKn"
    paramPAR <- "GrowthK"
  }
  
  if(param == "Linf" | param == "LinfN") {
    paramCTL <- "LinfN"
    paramPAR <- "Linf"
  }
  
  if(param == "SelPars50" | param == "sel50n") {
    paramCTL <- "sel50n"
    paramPAR <- "SelPars50"
  }
  
  if(param == "SelPars95" | param == "sel95n") {
    paramCTL <- "sel95n"
    paramPAR <- "SelPars50"
  }
  
  if(length(EMvec==1)) {
    col95 <- "grey50"
    col50 <- "grey75"
    colmed <- "grey25"
  }

  if(length(EMvec)>1) {
    if(length(EMvec)==2) {Ncols<-3}
    else{Ncols<-length(EMvec)}
    col95 <- adjustcolor(RColorBrewer::brewer.pal(Ncols,"Set1"),alpha=.25)
    col50 <- adjustcolor(RColorBrewer::brewer.pal(Ncols,"Set1"),alpha=.5)
    colmed <-adjustcolor(RColorBrewer::brewer.pal(Ncols,"Set1"),alpha=1)
  }
    
  # For own checking purposes in case ReadCTL variables change
  OM <- ReadCTLfile(OMfile)
  if(!paramCTL %in% names(OM$OM)) {stop("paramCTL not found in .par file. Check variable names from ReadCTLfile() output.")}
  SimYear <- OM$OM$SimYear
  dSimYear <- SimYear-1
  Nsim <- OM$OM$Nsim
  temp <- paste0("OM$OM$",paramCTL)
  truth <- CleanInput(eval(parse(text=temp)),SimYear)
  
  estimates <- matrix(ncol=5)
  colnames(estimates) <- c("EM","Year","Mean","Quant50","Quant95")
  temp <- matrix(nrow=Nsim,ncol=dSimYear)
  for(i in seq_along(EMvec)) {
    for(s in 1:Nsim) {
      EMdir <- file.path(getwd(),EMvec[i],s,SimYear)
      EM <- ReadPar(EMdir)
      linenums <- grep(paramPAR,names(EM))
      if(length(linenums)==0) {stop("paramPAR not found in PAR file. Check variable names from simass.par file.")}
      # BETTER WAY OF DOING THIS?
      if(length(linenums)==1) {
        if(length(grep("Sel",paramPAR))==1) {
          linenums <- c(linenums,grep(paste0("dev",strsplit(paramPAR,split="SelPars")[[1]][2]),names(EM)))
          avg_param <- EM[[linenums[1]]]
          dev_param <- EM[[linenums[2]]]
          temp[s,] <- avg_param+dev_param
          }
        if(paramPAR=="NatM") {
          linenums <- c(linenums,grep("m_dev",names(EM)))
        }
      }
      if(length(grep("SelPars",paramPAR))==0) {
        avg_param <- EM[[linenums[1]]]
        dev_param <- EM[[linenums[2]]]
        temp[s,] <- exp(avg_param+dev_param)
      }
    }
    
    tempmat <- cbind(
                rep(i,dSimYear),
                1:dSimYear,
                apply(temp,2,median,na.rm=T),
                apply(temp,2,quantile,na.rm=T,probs=.25),
                apply(temp,2,quantile,na.rm=T,probs=.75),              
                apply(temp,2,quantile,na.rm=T,probs=.05),
                apply(temp,2,quantile,na.rm=T,probs=.95)
                )
    
     if(i == 1) {
       estimates <- tempmat
     }
       
     else {
       estimates <- rbind(estimates,tempmat)
     }
  }
  
  colnames(estimates) <- c("EM","Year","Median","Quant25","Quant75","Quant05","Quant95")
  estimates <- data.frame(estimates)
  yearcoords <- c(1:dSimYear,dSimYear:1)

  if(is.na(ylims)) ylims <- range(c(estimates$Quant95,estimates$Quant05))
  plot(0,type='n',ylim=ylims,xlim=range(estimates$Year),...)
  for(em in seq_along(EMvec)) {
    tempdat <- estimates[estimates$EM==em,]
    polygon(x=yearcoords,y=c(tempdat$Quant95,rev(tempdat$Quant05)),col=col95[em],border=F)
    polygon(x=yearcoords,y=c(tempdat$Quant25,rev(tempdat$Quant75)),col=col50[em],border=F)
    lines(x=1:dSimYear,y=tempdat$Median,lwd=2,col=colmed[em])
  }
  lines(x=1:SimYear,y=truth,lwd=2,lty=2,col="black")
  if(is.na(plotnames)) plotnames <- EMvec
  if(dolegend) {legend("topleft",c("Truth",plotnames),lty=c(2,rep(1,length(EMvec))),col=c("black",colmed),lwd=2)}
  retlist <- list(truth=truth,estimates=estimates)

}
