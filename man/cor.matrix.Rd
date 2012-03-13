\name{cor.matrix}
\alias{cor.matrix}
\alias{gcc.cor.matrix}
\alias{gcc.cormatrix}

\title{
	correlation calculation for a data matrix
}

\description{
	This function can be used to compute the correlations of microarray and RNA-Seq gene expression profiles for a set of genes. The results can be further used to infer transcriptional regulation relationships between transcriptor factors and target genes, perform clustering analysis, and reconstruct biological networks.
}


\usage{cor.matrix(GEMatrix,
      cpus = 1,
      cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
      style = c("all.pairs", "pairs.between", "adjacent.pairs", "one.pair"),
      var1.id = NA,
      var2.id = NA,
      pernum = 0,
      sigmethod = c("two.sided", "one.sided"),
      output = c("matrix", "paired"))
}


%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GEMatrix}{
	a data matrix containing numeric variables. Example: rows may correspond to genes and columns to samples. 
}
  \item{cpus}{
	the number of cpus used for correlaiton calcluation. snowfall package in R needed to be installed in advance.
}
  \item{cormethod}{
	a character string that specifies a correlation method to be used for correlation calculation. Currently, the available correlation methods include: "GCC" (Gini correlation", "PCC" (Pearson's product-moment correlation), "SCC" (Spearman's rank correlation), "KCC"(Kendall's rank correlation), "BiWt" (Tukey's biweight correlation). 
}
  \item{style}{
	a character string that indicates the global or partial individuals (e.g., genes) to be used for correlation calculation.
}
  \item{var1.id}{
	a numeric vector specifying the variables to form gene pairs. For the "all.pairs" style, this parameter can be defined by the program itself. For styles of "adjacent.pairs", "pairs.between" and "adjacent.pairs", this parameter MUST be pre-defined.
}
  \item{var2.id}{
	a numeric vector indicating the variables to form gene pairs. The definiation of this parameter is similar to var1.id
}
  \item{pernum}{
	the number of permutation test used for calcluating statistical significance level (i.e., p-value) of correlations. 
}
  \item{sigmethod}{
	a character string ("two-sided" or "one-sided") that specifies the method used to compute p-value for permutation test. 
}
  \item{output}{
	a character string ("matrix" or "paired") that represents the output format of correlation data. The "matrix" outputs two matrix for correlations and p-values, respectively. The "paired" outputs only one matrix for each row with the information of paired variables, the correlation and p-value.
}

}



\details{
	Given a data matrix (e.g., microarray and RNA-Seq gene expression matrix), calculating correlation with GCC and other correlation methods for partial(or all) individuals (e.g., genes). The statistical significance (i.e., p-value) of each correlation is derived from the permutation test. Parallel computing options are also provided for speeding up correlation calculation.  	

}



\value{
	A list with the following components:
  \item{corMatrix }{correlations of gene pairs shown in matrix form. This data matrix is generated only when the output format "matrix" is specified. }
  \item{pvalueMatrix }{p-values of correlation shown in matrix form. This data matrix is generated only when the output format "matrix" is specified. }
  \item{corpvalueMatrix }{correlation and p-values listed in one form. This data matrix is generated only when the output format "paired" is indicated. }
}




\author{
Chuang Ma, Xiangfeng Wang
}


\note{
(1) var1.id and var2.id should be defined with the numeric vector format for "pairs.between", "adjacent.pairs" or "one.pair" styles. For instance, var1.id = c(1), var2.id = c(2).

(2) To perform BiWt, the R package "biwt" should be installed in advance. 

(3) To perform the parallel compution, the "snowfall" package in R should be first installed.
}


\seealso{
	\code{\link{cor.pair}}, \code{\link{onegcc}}, \code{\link{cor.test}}.
}



\examples{
   data(rsgcc)
   x <- rnaseq[1:10,]    #Just use a small subset of RNA-Seq data 

\dontrun{
   #run on one CPU.
   cor.matrix(x, cpus = 1, 
              cormethod = "GCC", style = "all.pairs", 
              pernum = 2000, sigmethod = "two.sided", 
              output = "matrix")  

   #run on two CPUs, snowfall package should be properly installed.
   cor.matrix(x, cpus = 2, 
              cormethod = "GCC", style = "all.pairs", 
              pernum = 2000, sigmethod = "two.sided", 
              output = "matrix")  

}

}

\keyword{ correlation }

