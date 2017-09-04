#' Log Chebyshev Approximation
#'
#'\code{logcheby} creates input for sqlp to solve the Chebyshev Approximation Problem
#'
#'@details
#' Solves the log Chebyshev approximation problem. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param B A pxm real valued matrix with p > m
#' @param f A vector of length p
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
#' data(Blogcheby)
#' data(flogcheby)
#'
#' out <- logcheby(Blogcheby, flogcheby)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
logcheby <- function(B,f){
  
  #Error Checking
  stopifnot(is.matrix(B), is.numeric(B), is.numeric(f), nrow(B) == length(f), ncol(as.matrix(f)) == 1)
  
  for(i in 1:ncol(B)){
    stopifnot(all(B[,i]/f > 0))
  }
  
  #Define Variables
  blk <- matrix(list(),2,2)
  At <- matrix(list(),2,1)
  C <- matrix(list(),2,1)
  
  p <- nrow(B)
  m <- ncol(B)
  
  blk[[1,1]] <- "s"
  blk[[1,2]] <- 2 * matrix(1,1,p)
  blk[[2,1]] <- "l"
  blk[[2,2]] <- p
  
  E <- matrix(0,p,m)
  for(j in 1:m){
    E[,j] <- B[,j]/f
  }
  if(any(E < 1e-10)){
    stop("B/f must have all positive entries")
  }
  beta <- matrix(0,p,1)
  for(i in 1:p){
    beta[i] <- sum(E[i,])
  }
##
  temp <- matrix(0,2*p+1,1)
  temp[seq(2,2*p,2)] <- rep(1,p)
  C[[1,1]] <- Matrix(0,2*p,2*p)
  
  for(i in 1:(nrow(C[[1,1]])-1)){
    C[[1,1]][i,i+1] <- temp[i+1]
  }
  
  C[[1,1]] <- C[[1,1]] + t(C[[1,1]])
  C[[2,1]] <- matrix(0,p,1)
  b <- matrix(c(rep(0,m),1),m+1,1)
  
  A <- matrix(list(),2,m+1)
  temp <- matrix(0,2*p,1)
  for(k in 1:m){
    temp[seq(1,2*p-1,2)] <- -E[,k]
    A[[1,k]] <- Matrix(0,2*p,2*p)
    diag(A[[1,k]]) <- c(temp)
    A[[2,k]] <- as.matrix(E[,k])
  }
  temp <- matrix(0,2*p,1)
  temp[seq(2,2*p,2)] <- rep(1,p)
  A[[1,m+1]] <- Matrix(0,2*p,2*p)
  diag(A[[1,m+1]]) <- c(temp)
  A[[2,m+1]] <- rep(1,p)
  
  At <- svec(blk,A,matrix(1,nrow(blk),1))
  
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
}