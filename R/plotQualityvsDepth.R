#'@depend ggplot2
#'@depend reshape2
#'@depend ShortRead
#'@depend data.table
#'@name plotQualityvsDepth.R
#'@title Plot quality vs depth for an accession
#'@param acc Accension number
#'@param ref reference name
#'@param start Position start
#'@param stop1 Position end
#'@param depth Depth filter
#'@details stuff
#'@export
plotQualityVsDepth <- function(acc,ref,start  ,stop1 ,depth ) {
  
  l <- getPileUp(acc,ref,start ,stop1 ,depth )
  l$NumericalQuality <- sapply( l[,6], function(x) alphabetScore(FastqQuality(as.character(x)))/ width(FastqQuality(as.character(x))))
  l.melted <- melt(l, measure.vars = c("NumericalQuality", "PileupDepth" ))
  ggplot(l.melted, aes(ReferencePosition, value, color = variable)) + geom_line()
  
}

#plotQualityVsDepth('SRR390728','21',19152215 ,19152615 ,30)
#plotQualityVsDepth("SRR1596669",'21',19152215 ,19152615 ,30)

# testAcc <- read.table(file = '/home/user2/test/acc.test', stringsAsFactors = F)
#blah <- lapply(testAcc[,1], function(x) getPileUp(x, '21',19152215 ,19152615 ,0))
#blah <- rbindlist(blah)
#ggplot(blah, aes(ReferencePosition,PileupDepth, color = AccensionNumber)) + geom_line() + scale_y_continuous(limits = c(0, 30))



#plotQualityVsDepth('SRR390728','21',1 ,19152615 ,30)




