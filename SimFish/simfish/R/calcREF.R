#' \code{calcREF} calculates reference points
#'
#' @param SRouts presumably spawner per recruit outs
#'
#' @return unclear
#' @export
calcREF<-function(SRouts)
{
  if(SRouts[[1]][1]!="Average")
  {
    testExp	<-seq(0.001,1,0.01)
    recYield	<-rep(0,length(testExp))
    projYrs	<-100
    tempN	<-matrix(nrow=projYrs,ncol=maxAge)
    tempN[1,]	<-virInit
    tempYield	<-rep(0,length(testExp))
    for(f in 1:length(testF))
    {
      tempExp<-testExp[f]
      for(j in 2:projYrs)
      {
        for (i in 2:(maxAge-1))
          tempN[j,i]		<-tempN[j-1,i-1]*(1-tempExp*vuln[i-1])*survival
        tempN[j,maxAge]	<-(tempN[j-1,(maxAge-1)])*(1-tempExp*vuln[maxAge])*survival + tempN[j-1,maxAge]*(1-tempExp*vuln[maxAge])*survival
        spbioIn			<-sum(tempN[j-1,]*mature*weight)
        if(SRouts[[1]][1]=="BevH")
          tempN[j,1]		<-calcBHrec(SRouts,spbioIn)
        if(SRouts[[1]][1]=="Rick")
          tempN[j,1]		<-calcRICKrec(SRouts,spbioIn)
      }
      tempYield[f]		<-sum(tempN[j,]*tempExp*vuln*weight)
    }
  }
  calcRICKrec(SRouts,virBio)
  if(SRouts[[1]][1]=="Average")
  {
    #STUFF
  }

}