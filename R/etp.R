#' Educational Testing Problem
#'
#'\code{etp} creates input for sqlp to solve the Educational Testing Problem -
#'given a symmetric positive definite matrix S, how much can be subtracted from the diagonal
#'elements of S such that the resulting matrix is positive semidefinite definite.
#'
#'@details
#' Solves the education testing problem. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param B A symmetric positive definite matrix
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
#' data(Betp)
#' 
#' out <- etp(Betp)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
etp <- function(B){
  
  #Error Checking
  stopifnot(is.matrix(B), is.numeric(B), isSymmetric(B,check.attributes = FALSE))
  
  #Define Variables
  n <- max(dim(B))
  
  blk <- matrix(list(),2,2)
  C <- matrix(list(),2,1)
  At <- matrix(list(),2,1)
  A <- matrix(list(),2,n)
  
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n
  blk[[2,1]] <- "l"
  blk[[2,2]] <- n
  
  b <- matrix(1,n,1)
  C[[1,1]] <- B
  C[[2,1]] <- matrix(0,n,1)
  
  for(k in 1:n){
    A[[1,k]] <- Matrix(0,n,n)
    A[[1,k]][k,k] <- 1
    A[[2,k]] <- rbind(matrix(0,k-1,1),-1,matrix(0,n-k,1))
  }
  
  At <- svec(blk,A,matrix(1,nrow(blk),1))
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
  
}