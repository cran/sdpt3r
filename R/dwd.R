#' Distance Weighted Discrimination
#'
#'\code{dwd} creates input for sqlp to solve the Distance Weighted Discrimination problem -
#'Given two sets of points An and Ap, find an optimal classification rule to group the points as accurately
#'as possible for future classification.
#'
#'@details
#'
#' Solves the distance weighted discrimination problem. Mathematical and implementation
#' details can be found in the vignette
#' 
#' @param Ap An nxp point configuration matrix
#' @param An An nxp point configuration matrix
#' @param penalty A real valued scalar penalty for moving points across classification rule
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
#' data(Andwd)
#' data(Apdwd)
#' penalty <- 0.5
#' 
#' out <- dwd(Apdwd,Andwd,penalty)
#' blk <- out$blk
#' At <- out$At
#' C <- out$C
#' b <- out$b
#' 
#' @export
dwd <- function(Ap,An,penalty){
  
  #Error Checking
  stopifnot(is.matrix(Ap), is.matrix(An), is.numeric(Ap), is.numeric(An), is.numeric(penalty), nrow(Ap)==nrow(An), ncol(Ap)==ncol(An))
  
  #Data input as nxp, program needs it to be pxn
  Ap <- t(Ap)
  An <- t(An)
  
  #Define Variables
  np <- nrow(Ap)
  mp <- ncol(Ap)
  nn <- nrow(An)
  mn <- ncol(An)
  n <- np
  nv <- 1 + n + 2 + 3*(mp + mn) + mp + mn
  nc <- 1 + 2*(mp + mn)
##
  blk <- matrix(list(),2,2)
  At <- matrix(list(),2,1)
  C <- matrix(list(),2,1)
##  
  
  blk[[1,1]] <- "q"
  blk[[1,2]] <- matrix(c(n+1,2,3*rep(1,mp+mn)),nrow=1)
  blk[[2,1]] <- "l"
  blk[[2,2]] <- mp + mn
##
  A <- Matrix(0,nc,nv-mp-mn,sparse=TRUE)
  A[1:mp,2:(n+3)] <- cbind(t(Ap),matrix(0,mp,1),matrix(1,mp,1))
  A[(mp+1):(mp+mn), 2:(n+3)] <- -cbind(t(An),matrix(0,mn,1), matrix(1,mn,1))
  A[1:mp,seq(n+4,n+3+3*mp,3)] <- Diagonal(mp,-1)
  A[1:mp,seq(n+6,n+5+3*mp,3)] <- Diagonal(mp,-1)
  A[(mp+1):(mp+mn),seq(3*mp+n+4,3*mp+n+3+3*mn,3)] <- Diagonal(mn,-1)
  A[(mp+1):(mp+mn),seq(3*mp+n+6,3*mp+n+5+3*mn,3)] <- Diagonal(mn,-1)
  A[mp+mn+1,1] <- 1
  A[(mp+mn+2):(mp+mn+1+mp),seq(n+5,n+4+3*mp,3)] <- Diagonal(mp,1)
  A[(mp+mn+1+mp+1):(mp+mn+1+mp+mn),seq(3*mp+n+5,3*mp+n+4+3*mn,3)] <- Diagonal(mn,1)

  At[[1,1]] <- A
  At[[2,1]] <- rbind(as.matrix(Diagonal(mp+mn,1)),matrix(0,1+mp+mn,mp+mn))
  b <- rbind(matrix(0,mp+mn,1),matrix(1,1+mp+mn,1))
##
  ctmp <- matrix(0,nv-mp-mn,1)
  ctmp[seq(n+4,n+3+3*mp,3)] <- rep(1,mp)
  ctmp[seq(n+6,n+5+3*mp,3)] <- -rep(1,mp)
  ctmp[seq(3*mp+n+4,3*mp+n+3+3*mn,3)] <- rep(1,mn)
  ctmp[seq(3*mp+n+6,3*mp+n+5+3*mn,3)] <- -rep(1,mn)
  
  C[[1,1]] <- ctmp
  C[[2,1]] <- penalty*matrix(1,mp+mn,1)
##
  output <- list(blk=blk, At=At, b=b, C=C, OPTIONS = list())
  class(output) <- "sqlp_input"
  
  return(output)
}