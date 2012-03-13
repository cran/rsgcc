\name{cor.pair}
\alias{cor.pair}
\alias{gcc.cor.pair}
\alias{gcc.corpair}

\title{
	compute correlations for a gene pair
}



\description{
	This function computes the correlation with Gini correlation and other correlation methods for a pair of individuals (e.g., gene pair), and calculates the statistical significance of correlations with the permutation test method.
}


\usage{
cor.pair(idxvec, 
      GEMatrix,
      rowORcol = c("row", "col"),
      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
      pernum = 0,
      sigmethod = c("two.sided", "one.sided"))
}



\arguments{
  \item{idxvec}{
	a numer vector containing two elements for the index of genes or samples in GEMatrix.
  }
  \item{GEMatrix}{
	a data matrix containing numeric variables. Example: rows may correspond to genes and columns to samples. 
  }
  \item{rowORcol}{
	a character string ("row" or "col") indicating the gene expression profiles are extracted by rows("row") or columns ("col"). 
  }
  \item{cormethod}{
	a character string that specifies a correlation method to be used for correlation calculation. Currently, the available correlation methods include: "GCC" (Gini correlation", "PCC" (Pearson's product-moment correlation), "SCC" (Spearman's rank correlation), "KCC"(Kendall's rank correlation), and "BiWt" (Tukey's biweight correlation). 
  }
  \item{pernum}{
	the number of permutation test used for calcluating statistical significance level (i.e., p-value) of correlations. 
  }
  \item{sigmethod}{
	a character string ("two-sided" or "one-sided") specifying the method used to compute p-value for permutation test.
  }
}



\value{
	A list with the following components:
  \item{gcc.rankx }{a Gini correlation produced by using the rank information of the first gene (i.e., the first element in idxvec).}
  \item{gcc.ranky }{a Gini correlation produced by using the rank information of the second gene (i.e., the second element in idxvec).}
  \item{gcc.rankx.pvalue }{p-value of gcc.rankx.}
  \item{gcc.ranky.pvalue }{p-value of gcc.ranky.}
  \item{cor }{the correlation produced by "PCC", "SCC", "KCC" or "BiWt".}
  \item{pvalue }{the p-value of cor.}

}


\author{
Chuang Ma, Xiangfeng Wang
}


\note{
(1) To perform BiWt, the R package "biwt" should be installed in advance. 

(2) When the cormethod defined as "GCC", this function will output a list with four numeric elements: gcc.rankx, gcc.ranky, gcc.rankx.pvalue, gcc.ranky.pvalue.
}



\seealso{
   \code{\link{onegcc}}, \code{\link{cor.matrix}}, \code{\link{gcc.corfinal}}.
}




\examples{

   data(rsgcc)
   x <- rnaseq[1:50,]

   #compute correlation between 47th and 50th genes
   corpair <- cor.pair(c(47,50), GEMatrix = x, rowORcol = "row", 
                       cormethod = "GCC", pernum = 2000, 
                       sigmethod = "two.sided")  
}

\keyword{ correlation }

