% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toeplitz.R
\name{toep}
\alias{toep}
\title{Toeplitz Approximation Problem}
\usage{
toep(A)
}
\arguments{
\item{A}{A symmetric matrix}
}
\value{
\item{X}{A list containing the solution matrix to the primal problem}
\item{y}{A list containing the  solution vector to the dual problem}
\item{Z}{A list containing the  solution matrix to the dual problem}
\item{pobj}{The achieved value of the primary objective function}
\item{dobj}{The achieved value of the dual objective function}
}
\description{
\code{toep} creates input for sqlp to solve the Toeplitz approximation problem -
given a symmetric matrix F, find the nearest symmetric positive definite Toeplitz matrix.
}
\details{
For a symmetric matrix A, determines the closest Toeplitz matrix. Mathematical and implementation
details can be found in the vignette
}
\examples{
data(Ftoep)

#Not Run
#out <- toep(Ftoep)

}
