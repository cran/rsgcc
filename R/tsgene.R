##################################################################
##find tissue specific genes

if( !require(stringr) ) install.packages("stringr")
require(stringr)

uniqueTissues <- function(x) {
  
    sampleNum <- ncol(x)
    if( is.null(colnames(x)) ) {  #no annoted tissue sample, each sample belongs to different tissue by default
      tsMatrix <- diag(x = 1, nrow = sampleNum, ncol = sampleNum ) 
    }else {
      tsMatrix <- matrix(0, nrow = sampleNum, ncol = sampleNum )
      colnames(tsMatrix) <- colnames(x)
      uniTS <- c("")
      for( i in 1: sampleNum ) {
        lastdot <- sapply(gregexpr("\\.", colnames(x)[i]), tail, 1)
        if( lastdot < 0 )  { curTSName <- colnames(x)[i] 
        } else { 
          if( lastdot > 1 ) {curTSName <- str_sub( colnames(x)[i], 1, lastdot - 1) 
        }else {
          curTSName <- colnames(x)[i]
        }
      }
      
      curTSIndex <- which(uniTS == curTSName)
      if( length(curTSIndex) == 0 ) { #not record
        if( i == 1) uniTS[1] <- curTSName
        else        uniTS <- c(uniTS, curTSName)
        
        tsMatrix[length(uniTS), i] <- 1
      }else {
        tsMatrix[curTSIndex, i] <- 1
      } 
    }#end for i
    tsMatrix <- tsMatrix[1:length(uniTS),]
    rownames(tsMatrix) <- uniTS
  }#end else
    
    tsMatrix
  
}


##################################################################################
##ts gene: t measure, 1- max(nonTissue)/max(Tissue)
getsgene <- function(x, Log = FALSE, Base = 2, AddOne = FALSE, tsThreshold = 0.95, Fraction = TRUE ) {
  
  ##remove all zeros
  zeroIndex <- which( apply(x, 1, function(vec) length(which(vec==0))) == ncol(x) )
  if( length(zeroIndex) > 0) {  x <- x[-zeroIndex,]  }
   
  if(AddOne) { x <- x + 1  }
  if(Log) { x <- log(x, Base) }  
 
  #check one gene
  onets <- function(vec, tsMatrix) {
    tscorematrix <- matrix(0, nrow = dim(tsMatrix)[1], ncol = 2)
    for( i in 1: dim(tsMatrix)[1] ) {
      sampleIndex <- which(tsMatrix[i,] > 0)  #sample index for cur tissue
      meanvalue <- mean(vec[sampleIndex])
      tscorematrix[i,1] <- i
      if( meanvalue == 0 ) { tscorematrix[i,2] <- -100}
      else { tscorematrix[i,2] <- 1 - max(vec[-sampleIndex])/meanvalue }
    }#end for i
    tmax = max(tscorematrix[,2])
    tmaxidx = tscorematrix[which(tscorematrix[,2] == tmax),1][1]
    return( list( tmaxidx = tmaxidx, tmax = tmax ) )
  }
  tsMatrix <- uniqueTissues(x)  
  tt <- apply(x, 1, onets, tsMatrix = tsMatrix )
  tscorematrix <- matrix(0, nrow = dim(x)[1], ncol = 3)
  colnames(tscorematrix) <- c("GeneIndex", "tsmaxscore", "tsmaxidx")
  for( i in 1:dim(x)[1]) {
    tscorematrix[i,1] <- i
    tscorematrix[i,2] <- tt[[i]]$tmax
    tscorematrix[i,3] <- tt[[i]]$tmaxidx    
  }#end for i
  
  tscore <- tscorematrix[which(tscorematrix[,2] >= tsThreshold),]
  tsgene <- x[tscore[,1],]
  if(Fraction) {   
    tsgene <- t(apply(tsgene, 1, function(vec) vec/sum(vec)) )
  }
  
  return( list( tsgene = tsgene, tscore = tscorematrix, uniquets = tsMatrix) ) 
}
 
