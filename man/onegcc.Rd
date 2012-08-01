\name{onegcc}
\alias{onegcc}

\title{
	compute one Gini correlation coefficient
}
\description{
	onegcc calcluates one Gini correlation coefficient with rank information of the first variable.
}
\usage{
	onegcc(x, y)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector. }
  \item{y}{ a numeric vector with the same length of x. }
}

\details{
	This is a generic function cacluating correlation with rank information of the first variable and the actual value information of the second variable.
}
\value{
	Gini correlation coefficient (a numeric value ranged from -1.0 to 1.0).
}


\author{
	Chuang Ma, Xiangfeng Wang
}


\examples{

   data(rsgcc)
   x <- rnaseq[1:10,]     #Just use a small subset of RNA-Seq data 
   onegcc(x[1,], x[2,])   # generate one correlaiton for one gene pair
   onegcc(x[2,], x[1,])   # generate the other correlaiton for the same gene pair
}


\keyword{ correlation }
