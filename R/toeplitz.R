#' Toeplitz Approximation Problem
#'
#'\code{toep} creates input for sqlp to solve the Toeplitz approximation problem -
#'given a symmetric matrix F, find the nearest symmetric positive definite Toeplitz matrix.
#'
#'@details
#' For a symmetric matrix A, determines the closest Toeplitz matrix. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param A A symmetric matrix
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
#' data(Ftoep)
#' 
#' out <- toep(Ftoep)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
toep <- function(A){
  
  #Error Checking
  stopifnot(is.matrix(A), is.numeric(A), nrow(A) == ncol(A), isSymmetric(A,check.attributes = FALSE))
  
  #Define Variables
  
  n <- max(dim(A))
  gam <- sqrt(c(n,2*seq(n-1,1,-1)))
  q <- matrix(0,n,1)
  q[1] <- -sum(diag(A))
  for(k in 1:(n-1)){
    tmp <- c()
    #Get kth diagonal
    for(i in 1:(nrow(A)- k)){
      tmp <- c(tmp,A[i,i+k])
    }
    q[k+1] <- -2*sum(tmp)
  }
  beta <- norm(A,type="F")^2
  
  blk <- matrix(list(),2,2)
  C <- matrix(list(),2,1)
  At <- matrix(list(),2,1)
  
  blk[[1,1]] <- "s"
  blk[[1,2]] <- n+1
  blk[[2,1]] <- "s"
  blk[[2,2]] <- n+1
  
  b <- matrix(c(rep(0,n),-1),ncol=1)
  
  C[[1,1]] <- Matrix(0,n+1,n+1,sparse=TRUE)
  C[[2,1]] <- Diagonal(n+1,c(rep(1,n),-beta))
  
  Acell <- matrix(list(),1,n+1)
  Acell[[1]] <- -Diagonal(n+1,c(rep(1,n),0))
  tmpvec <- c(rep(-1,n),0)
  for(k in 1:(n-1)){
    tmp <- Matrix(0, n+1,n+1,sparse=TRUE)
    for(j in 1:(nrow(tmp)-k)){
      tmp[j,j+k] <- tmpvec[j+k]
    }
    Acell[[k+1]] <- tmp + t(tmp)
  }
  Acell[[n+1]] <- Matrix(0,n+1,n+1,sparse=TRUE)
  Acell[[n+1]][n+1,n+1] <- -1
  At[[1,1]] <- svec(blk[1,,drop=FALSE],Acell,1)[[1]]
  
  for(k in 1:n){
    Acell[[k]] <- Matrix(0,n+1,n+1)
    Acell[[k]][k,n+1] <- -gam[k]
    Acell[[k]][n+1,k] <- -gam[k]
    Acell[[k]][n+1,n+1] <- 2*q[k]
  }
  At[[2,1]] <- svec(blk[2,,drop=FALSE],Acell,1)[[1]]
  
  return(list(blk=blk,At=At,C=C,b=b))
  
}