#' Graph Partitioning Problem
#' 
#'\code{gpp} creates input for sqlp to solve the graph partitioning problem.
#'
#'@details
#'
#' Solves the graph partitioning problem. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param B A weighted adjacency matrix
#' @param alpha Any real value in (0,n^2)
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
#' data(Bgpp)
#' alpha <- nrow(Bgpp)
#'
#' out <- gpp(Bgpp, alpha)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#'
#' @export
gpp <- function(B, alpha){
  
  #Error Checking
  stopifnot(is.matrix(B), is.numeric(B), isSymmetric(B,check.attributes = FALSE), ncol(B)==nrow(B), !all(B == 0), is.numeric(alpha), alpha > 0, alpha < nrow(B)^2)
  
  #Define Variables
  blk <- matrix(list(),1,2)
  C <- matrix(list(),1,1)
  At <- matrix(list(),1,1)
  
  n <- max(dim(B))
  e <- matrix(1,n,1)
  C[[1]] <- -(Diagonal(n,B%*%e) - B)
  b <- rbind(alpha,e)
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n
  
  A <- matrix(list(),1,n+1)
  A[[1]] <- e %*% t(e)
  for(k in 1:n){
    A[[k+1]] <- Matrix(0,n,n)
    A[[k+1]][k,k] <- 1
  }
  
  At <- svec(blk,A,1)
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
}