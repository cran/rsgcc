\name{adjacencymatrix}
\alias{adjacencymatrix}


\title{
adjacency matrix calculation
}
\description{
This function generates the adjacency matrix for network re-construction from gene expression data with different methods including Gini correlation (GCC), Pearson correlation (PCC), Spearman correlation (SCC), Kendall correlation (KCC) and mutual information (MI) methods. It was implemented these methods in C language and parallel mode, and thus is greatly faster than the cor.matrix function.
}
\usage{
adjacencymatrix(mat, method = c("GCC", "PCC", "SCC", "KCC", "MI"), k = 3)
}


\arguments{
  \item{mat}{a data matrix containing gene expression dataset where rows defines for genes and columns for samples.}
  \item{method}{a method used to calculate the association between a pair of genes.}
  \item{k}{the number of nearest neighbors to be considered for estimating the mutual information. Must be less than the number of columns of mat, and only work for the mutual information(MI) method.}
}


\value{
  \item{value}{a matrix recording the associations between the gene pairs.}
}


\references{
[1] Ma C and Xiang XF. Application of the Gini correlation coefficient to infer regulatory relationships in transcriptome analysis, Plant Physiology, 2012, 160(1):192-203.

[2] Sales G and Romualdi C. Parmigene-a parallel R package for mutual information estimation and gene network reconstruction. Bioinformatics, 2012, 27:1876-1877.

}


\note{
1) the mutural information estimation is based on k-nearest neighbor distance (Sales G and Romualdi C, 2012). Thus the parameter k only works for the mutual information method.

2) Two correlations can be produced by the GCC method by reciprocally using the rank and value information of one gene (or variable). Here the correlation with the maximum absolute values is selected for generating the adjacency matrix. 

3) The diagonal of the adjacency matrix is 0. 
}



\examples{

\dontrun{
 mat = matrix(rnorm(1800), nrow = 100)
 mat

 #for GCC method
 adjacencymatrix( mat, method = "GCC" )

 for MI method, k works here.
 adjacencymatrix( mat, method = "MI", k= 3)
}

}


\keyword{ correlation }

