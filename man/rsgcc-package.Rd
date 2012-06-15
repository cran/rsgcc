\name{rsgcc-package}
\alias{rsgcc-package}
\alias{rsgcc}
\docType{package}
\title{
	Gini methodology-based correlation and clustering analysis of microarray and RNA-Seq gene expression data
}
\description{
	This package provides functions for calcluating the Gini, the Pearson, the Spearman, the Kendall and Tukey's Biweight correlations, and for correlation and clustering analysis of gene expression data profiled by the microarray and RNA-Seq technologies. Compared to other correlation methods, the GCC may perform better to detect regulatory relationships from gene expression data. In addtion, the GCC also has some other advantageous merits. Please see details in refernece (Ma and Wang, 2012). This package also provides an graphical user interface (GUI) to perform clustering analysis of microarray and RNA-Seq data in a coherent step-by-step manner.
}


\details{
\tabular{ll}{
	Package: \tab rsgcc\cr
	Type: \tab Package\cr
	Version: \tab 1.0.2\cr
	Date: \tab 2012-06-15\cr
	License: \tab GPL(>=2)\cr
}

}

\note{

1) The implement of rsgcc requires several R packages developed by other developers(e.g., biwt, cairoDevice, fBasics, snowfall, grDevices, gplots, gWidgets, gWidgetsRGtk2, stringr, ctc). Please make sure these packages have been successfully installed before installing the rsgcc package.

2) A general method to install a new package on the computer is to use the command: install.package("packagename"). Some other methods can be found at http://math.usask.ca/~longhai/software/installrpkg.html. For the installation of ctc package, please use the biocLite resource with the following commends:

      source("http://bioconductor.org/biocLite.R")

      biocLite("ctc")

3) To run the "rsgcc.gui()" function, please do remember to select the GUI toolkit "gWidgetsRGtk2".

4) Bug reports and suggestions/questions can be sent to Chuang Ma (chuangma2006@gmail.com) or Xiangfeng Wang (xwang1@cals.arizona.edu).
	
}


\author{
	Chuang Ma, Xiangfeng Wang.
	Maintainer: Chuang Ma <chuangma2006@gmail.com>
}
\references{
	Chuang Ma, Xiangfeng Wang. (2012) Application of the Gini correlation coefficient to infer regulatory relationships in transcriptome analysis.
}


\keyword{ package }



