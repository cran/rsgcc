##########################################################################
##Function: Compute the distance (dissimilarity) between both rows.
##Notice: Here we only consider the distance between rows.
##        For distance between columns, run t(x) before devoke the DistFun
###########################################################################

gcc.dist <- function(x, 
                     cpus = 1,
                     cormethod = c("GCC", "PCC", "SCC", "KCC", "BiWt"),
                     distancemethod = c("Raw", "Abs", "Sqr")) {
  
  if( length(distancemethod) > 1 ) {
    stop("Error: only allow one distance method")
  }
  if( is.null(cormethod)) cormethod = "GCC"
  if( is.null(distancemethod)) distancemethod = "Raw"
  
  AllPairMatrix <- cor.matrix(x, cpus = cpus, cormethod= cormethod, style= "all.pairs", pernum= 0, sigmethod= "two.sided", output = "matrix")$corMatrix
 
  if( distancemethod == "Raw") {
    ad <- as.dist( 1- AllPairMatrix ) 
  }else if( distancemethod == "Abs") {
    ad <- as.dist( 1- abs(AllPairMatrix) )    
  }else if( distancemethod == "Sqr") {
    ad <- as.dist( 1-AllPairMatrix^2)
  }else {
    stop("Error: the distance method should be Raw, Abs, or Sqr")
  }
  
  return( list( dist = ad, pairmatrix = AllPairMatrix))
   
} 
