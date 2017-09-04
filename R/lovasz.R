#' Lovasz Number of a Graph
#'
#'\code{lovasz} creates input for sqlp to find the Lovasz Number of a graph
#'
#'@details
#' Finds the maximum Shannon entropy of a graph, more commonly known as the Lovasz number. 
#' Mathematical and implementation details can be found in the vignette
#' 
#' @param G An adjacency matrix corresponding to a graph
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
#' data(Glovasz)
#'
#' out <- lovasz(Glovasz)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
lovasz <- function(G){
  
  #Error Checking
  stopifnot(is.matrix(G), is.numeric(G),isSymmetric(G,check.attributes = FALSE), nrow(G) == ncol(G), !all(G == 0))
  
  #Define Variables
  blk <- matrix(list(),1,2)
  C <- matrix(list(),1,1)
  At <- matrix(list(),1,1)
  
  n <- max(dim(G))
  m <- sum(G[upper.tri(G)]) + 1
  e1 <- matrix(c(1,rep(0,m-1)),m,1)
  
  C[[1]] <- matrix(-1,n,n)
  b <- e1
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n
  
  A <- matrix(list(),1,m)
  A[[1]] <- Diagonal(n,1)
  cnt <- 2
  for(i in 1:(n-1)){
    idx <- which(G[i,(i+1):n] != 0)
    idx <- idx + i
    if(length(idx) > 0){
      for(j in 1:length(idx)){
        A[[cnt]] <- Matrix(0,n,n)
        A[[cnt]][i,idx[j]] <- 1
        A[[cnt]][idx[j],i] <- 1
        cnt <- cnt + 1
      }
    }
  }
  
  At <- svec(blk,A,1)
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
  
}