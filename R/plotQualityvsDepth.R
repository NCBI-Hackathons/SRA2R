

library(ggplot2)
library(reshape2)
library(ShortRead)

plotQualityVsDepth <- function(acc,ref,start  ,stop1 ,depth ) {
  
 l <- getPileUp(acc,ref,start ,stop1 ,depth )
 l$NumericalQuality <- sapply( l[,6], function(x) alphabetScore(FastqQuality(as.character(x)))/ width(FastqQuality(as.character(x))))
 l.melted <- melt(l, measure.vars = c("NumericalQuality", "PileupDepth" ))
 ggplot(l.melted, aes(ReferencePosition, value, color = variable)) + geom_line()

 }

#plotQualityVsDepth('SRR390728','21',19152215 ,19152615 ,30)




