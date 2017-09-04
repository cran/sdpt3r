#' Max-Cut Problem
#'
#'\code{maxcut} creates input for sqlp to solve the Max-Cut problem -
#'given a graph B, find the maximum cut of the graph
#'
#'@details
#' Determines the maximum cut for a graph B. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param B A (weighted) adjacency matrix corresponding to a graph
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
#' data(Bmaxcut)
#'
#' out <- maxcut(Bmaxcut)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
maxcut <- function(B){
  
  #Error Checking
  stopifnot(is.matrix(B), is.numeric(B), isSymmetric(B,check.attributes = FALSE), nrow(B) == ncol(B), !all(B == 0))
  
  #Define Variables
  n <- max(dim(B))
  e <- matrix(1,n,1)
  
  C <- matrix(list(), 1, 1)
  blk <- matrix(list(),1,2)
  A <- matrix(list(), 1,n)
  
  C[[1]] <- matrix(0,n,n)
  diag(C[[1]]) <- B %*% e
  C[[1]] <- -(C[[1]] - B)/4
  b <- e
  
  blk[[1,1]] <- 's'
  blk[[1,2]] <- n
  
  for(k in 1:n){
    A[[k]] <- Matrix(0,n,n)
    A[[k]][k,k] <- 1
  }
  
  Avec <- svec(blk,M=A,isspx=matrix(0,nrow(blk),1))
  
  output <- list(blk=blk, At=Avec, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
  
  
}