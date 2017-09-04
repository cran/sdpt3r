#' Max-kCut Problem
#'
#'\code{maxkcut} creates input for sqlp to solve the Max-kCut Problem -
#' given a graph object B, determine if a cut of at least size k exists.
#'
#'@details
#' Determines if a cut of at least size k exists for a graph B. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param B A (weighted) adjacency matrix
#' @param K An integer value, the minimum number of cuts in B
#' 
#' @return 
#' Returns an object of class sqlp_input, containing the following:
#' 
#' \item{blk}{A matrix object describing the block diagonal structure of the SQLP data}
#' \item{At}{A matrix object containing constraint matrices for the primal-dual problem}
#' \item{C}{A matrix object containing the constant c matrices in the primal objective function}
#' \item{b}{A vector containing the right hand side of the equality constraints in the primal problem}
#' \item{OPTIONS}{A list object specifying the value of parbarrier}
#' 
#' @examples 
#' data(Bmaxkcut)
#'
#' out <- maxkcut(Bmaxkcut,2)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
maxkcut <- function(B,K){
  
  #Error Checking
  stopifnot(is.matrix(B), is.numeric(B), isSymmetric(B,check.attributes = FALSE), nrow(B) == ncol(B), !all(B == 0), is.numeric(K), K > 1)
  
  #Define Variables
  blk <- matrix(list(),2,2)
  At <- matrix(list(),2,1)
  C <- matrix(list(),2,1)
  
  n <- max(dim(B))
  e <- matrix(1,n,1)
  n2 <- n*(n-1)/2
  
  C[[1]] <- -(1-1/K)/2 * (Diagonal(n,B%*%e)-B)
  b <- e
  
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n
  blk[[2,1]] <- "l"
  blk[[2,2]] <- n2
  
  A <- matrix(list(),1,n)
  
  for(j in 1:n){
    A[[j]] <- matrix(0,n,n)
    A[[j]][j,j] <- 1
  }
  
  Avec <- svec(blk[1,,drop=FALSE],A,1)
  tmp <- Diagonal(n*(n+1)/2,1)
  idx <- cumsum(c(1:n))
  Atmp <- tmp[,setdiff(1:(n*(n+1)/2),idx)]
  
  At[[1,1]] <- cbind(Avec[[1,1]], Atmp/sqrt(2))
  At[[2,1]] <- cbind(Matrix(0,n2,n), -Diagonal(n2,1))
  
  b <- rbind(b, -1/(K-1) * matrix(1,n2,1))
  C[[2,1]] <- matrix(0,n2,1)
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
}