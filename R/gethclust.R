#########################################################################
##Function: compute cluster for microarray and RNASeq gene expression
##          data with different dissimilarity methods.
##Author: Chuang Ma
##Date: 2012-02-16
#########################################################################


gcc.hclust <- function(x,
                       cpus = 1,
                       cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
                       distancemethod = c("Raw", "Abs", "Sqr"),
                       clustermethod = c("complete", "average", "median", "centroid", "mcquitty", "single", "ward") ) {
  
  if( length(cormethod) > 1 ) {
    stop("Error: only allow one correlation method")
  }
   
  if( length(distancemethod) > 1 ) {
    print(distancemethod)
    stop("Error: only allow one distance method")
  }
 
  if( length(clustermethod) > 1 ) {
    stop("Error: only allow one cluster method")
  } 
  
  if( is.null(cormethod) ) cormethod <- "GCC"
  if( is.null(distancemethod)) distancemethod <- "Raw"
  if( is.null(clustermethod)) clustermethod <- "complete"
  
    ddata <- gcc.dist(x, cpus = cpus, cormethod= cormethod, distancemethod = distancemethod )
    hcdata <- hclust(ddata$dist, method = clustermethod)
  
  return( list(hc = hcdata, dist = ddata$dist, pairmatrix = ddata$pairmatrix))
                       
}
