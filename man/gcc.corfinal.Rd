\name{gcc.corfinal}
\alias{gcc.corfinal}

\title{
	get the final correlaiton and p-value of Gini method
}
\description{
	Compare two correlations produced by GCC method for one gene pair, and chose one as the final output of GCC method.
}


\usage{
	gcc.corfinal(gcccor)
}


\arguments{
  \item{gcccor}{
	a list output by cor.pair function for GCC method.
}
}



\details{
	If the p-value is "NA", the correlation with absolute maximum value is selected; otherwise, the correlation with lower p-value is chosen.
}
\value{
  \item{gcc.fcor }{the final correlation of GCC.}
  \item{gcc.fpavlue }{the final pvalue of correlation.}
%% ...
}


\author{
Chuang Ma, Xiangfeng Wang
}


\seealso{
      \code{\link{onegcc}}, \code{\link{cor.pair}}.
}
\examples{
   data(rsgcc)
   x <- rnaseq[1:50,]

   #compute correlation between 4th and 7th genes
   corpair <- cor.pair(c(4,7), GEMatrix = x, rowORcol = "row", 
                       cormethod = "GCC", pernum = 2000, 
                       sigmethod = "two.sided") 

   #get the final correlation and p-value of GCC method 
   gcc.corfinal(corpair)  

}


\keyword{correlation}

