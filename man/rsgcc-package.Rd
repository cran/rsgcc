\name{rsgcc-package}
\alias{rsgcc-package}
\alias{rsgcc}
\docType{package}
\title{
	Gini-based correlation and clustering analysis of microarray and RNA-Seq gene expression data
}
\description{
	This package provides functions for calcluating Gini correlation coefficient (GCC) for normal and non-normal distributed data. Compared to Pearson correlation and Spearman correlation, GCC is more powerful to infer transcriptional regualtion relationships from microarray and RNA-Seq data, and is more robustness to outliers. This package also provides an graphical user interface (GUI) to perform clustering analysis of microarray and RNA-Seq data in a coherent step-by-step manner.
}


\details{
\tabular{ll}{
	Package: \tab rsgcc\cr
	Type: \tab Package\cr
	Version: \tab 1.0\cr
	Date: \tab 2012-03-13\cr
	License: \tab GPL(>=2)\cr
}

}

\note{

1) The implement of rsgcc requires several other R packages (e.g., biwt, ctc, stringr, fBasics, snowfall, gplots, gWidgetsRGtk2, et al.). Please properly install the packages required in the installation procession.

2) To install ctc package, please use the biocLite resource with the following commends:

      source("http://bioconductor.org/biocLite.R")

      biocLite("ctc")

3) To run the rsgcc.gui function, do remember to select the GUI toolkit "gWidgetsRGtk2".

4) Bug reports and suggestions/questions can be sent to Chuang Ma (chuangma2006@gmail.com).
	
}


\author{
	Chuang Ma, Xiangfeng Wang.
	Maintainer: Chuang Ma <chuangma2006@gmail.com>
}
\references{
	Chuang Ma, Xiangfeng Wang. (2012) Evaluation of correlation methods to infer regulatory relationships in gene expression analysis.
}


\keyword{ package }



