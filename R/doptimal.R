#'D-Optimal Experimental Design
#'
#'\code{doptimal} creates input for sqlp to solve the D-Optimal Experimental Design problem -
#'given an nxp matrix with p <= n, find the portion of points that maximizes det(A'A)
#'
#'@details
#' Solves the D-optimal experimental design problem. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param V a pxn matrix containing a set of n test vectors in dimension p (with p <= n) 
#' 
#' @return 
#' Returns an object of class sqlp_input, containing the following:
#' 
#' \item{blk}{A matrix object describing the block diagonal structure of the SQLP data}
#' \item{At}{A matrix object containing constraint matrices for the primal-dual problem}
#' \item{C}{A matrix object containing the constant C matrices in the primal objective function}
#' \item{b}{A vector containing the right hand side of the equality constraints in the primal problem}
#' \item{OPTIONS}{A list object specifying the value of parbarrier}
#' 
#' @examples 
#' data(DoptDesign)
#' 
#' out <- doptimal(DoptDesign)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' OPTIONS <- out$OPTIONS
#'
#' @export
doptimal <- function(V){
  
  #Error Checking
  stopifnot(is.matrix(V), nrow(V) <= ncol(V), is.numeric(V))
  
  #Define Variables
  blk <- matrix(list(), 3, 2)
  C <- matrix(list(), 3, 1)
  At <- matrix(list(), 3, 1)
  OPTIONS <- list(parbarrier = matrix(list(),3,1))
  
  n <- nrow(V)
  p <- ncol(V)
  
  b <- matrix(0, p, 1)
  
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n
  
  Ftmp <- matrix(list(),1,p)
  
  for(k in 1:p){
    Ftmp[[1,k]] <- -V[,k] %*% t(V[,k])
  }
  
  At[1] <- svec(blk[1,,drop=FALSE],Ftmp,1)
  C[[1,1]] <- Matrix(0,n,n,sparse=TRUE)
  
  blk[[2,1]] <- "l"
  blk[[2,2]] <- p
  
  At[[2,1]] <- -diag(1,p,p)
  C[[2,1]] <- matrix(0,p,1)
  
  blk[[3,1]] <- "u"
  blk[[3,2]] <- 1
  At[[3,1]] <- matrix(1, 1, p)
  C[[3,1]] <- 1
  
  OPTIONS$parbarrier[[1,1]] <- 1
  OPTIONS$parbarrier[[2,1]] <- 0
  OPTIONS$parbarrier[[3,1]] <- 0
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = OPTIONS)
  class(output) <- "sqlp_input"
  
  return(output)
}