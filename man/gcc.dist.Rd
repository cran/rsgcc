\name{gcc.dist}
\alias{gcc.dist}

\title{
	compute distance matrix for hierarchical clustering
}

\description{
	This function computes the distance between the rows of a data matrix with the specified distance method. 
}


\usage{
gcc.dist(x, 
      cpus = 1,
      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
      distancemethod = c("Raw", "Abs", "Sqr"))
}

\arguments{
  \item{x}{a data matrix containing numeric variables. Example: rows may correspond to genes and columns to samples.}
  \item{cpus}{the number of cpus used for computation.}
  \item{cormethod}{a character string indicating the correlation method to be used. Currently, the available correlation methods include: "GCC" (Gini correlation", "PCC" (Pearson's product-moment correlation), "SCC" (Spearman's rank correlation), "KCC"(Kendall's rank correlation), and "BiWt" (Tukey's biweight correlation).}

  \item{distancemethod}{a character string indicating the distance method to be used. Currently, three distance methods are available, include: "Raw" (1-cor)", "Abs" (1-|cor|), and "Sqr" (1-|cor|^2). }
  
}



\value{
	A list with the following components:
  \item{dist }{a data matrix containing the distances between different genes.}
  \item{pairmatrix }{a data matrix including the correlation between different genes.}
%% ...
}


\author{
Chuang Ma, Xiangfeng Wang
}


\seealso{
 \code{\link{cor.matrix}}, \code{\link{gcc.hclust}}, \code{\link{gcc.tsheatmap}}.
}



\examples{

   data(rsgcc)
   x <- rnaseq[1:10,]
   gcc.dist(x, cormethod = "GCC", distancemethod = "Raw", cpus = 1)

}

\keyword{ cluster }
