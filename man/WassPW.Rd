% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WassPW.R
\name{WassPW}
\alias{WassPW}
\title{Wasserstein pairwise warping for distribution-valued random processes}
\usage{
WassPW(Lt = NULL, Ly = NULL, Lh = NULL, Lq = NULL, optns = list())
}
\arguments{
\item{Lt}{A list of \eqn{n} vectors, where \code{Lt[[i]]} is a vector of 
length \eqn{m_i} holding the observation time points for the \eqn{i}-th 
distribution-valued process \eqn{X_i(\cdot)}. Each vector should be sorted in ascending order.}

\item{Ly}{A list of \eqn{n} lists, where \code{Ly[[i]]} is a list of 
\eqn{m_i} vectors and the \eqn{j}-th vector \code{Ly[[i]][[j]]} holds the sample of data points 
generated from (a noisy version of) the distribution \eqn{X_i(t_{ij})} with \eqn{t_{ij}} given by 
\code{Lt[[i]][j]}. 
Note that only one of the three, \code{Ly}, \code{Lh}, and \code{Lq}, needs to be input. 
If more than one of them are specified, \code{Ly} overwrites \code{Lh}, and \code{Lh} overwrites \code{Lq}.}

\item{Lh}{A list of \eqn{n} lists, where \code{Lh[[i]]} is a list of \eqn{m_i} lists 
and \code{Lh[[i]][[j]]} holds the histogram corresponding to (a noisy version of) 
the distribution \eqn{X_i(t_{ij})} with \eqn{t_{ij}} given by \code{Lt[[i]][j]}. 
Specifically, \code{Lh[[i]][[j]]} is a list of two components, 
\code{breaks} (or \code{mids}) holding the cell boundaries (or cell midpoints) and 
\code{counts} holding the counts falling in each cell.}

\item{Lq}{A list of \eqn{n} lists, where \code{Lq[[i]]} is a matrix or list 
holding the quantile functions corresponding to the \eqn{i}-th distribution-valued process. 
If \code{Lq[[i]]} is a matrix, the support of the \eqn{\sum_{i=1}^{n} m_i} quantile functions 
should be the same (i.e., \code{optns$qSup}), and \eqn{j}-th row of \code{Lq[[i]]} holds 
the quantile function corresponding to corresponding to (a noisy version of) 
the distribution \eqn{\{X_i(t_{ij})\}_{j=1}^{m_i}} with \eqn{t_{ij}} given by \code{Lt[[i]][j]}. 
If the quantile functions are evaluated on different grids, then \code{Lq[[i]]} should be a list, 
with \code{Lq[[i]][[j]]} being a list of two components \code{x} and \code{y} 
holding the support grid and the corresponding values of the quantile functions, respectively.}

\item{optns}{A list of options control parameters specified by \code{list(name=value)}, 
including all those available in \code{\link[frechet]{LocDenReg}} and \code{\link[fdapace]{WFDA}}, 
except for that the default \code{subsetProp} is 1 and default \code{kernelReg} is \code{'epan'}; 
as well as 
\describe{
\item{ngrid}{Number of grid points for evaluating the warping functions (default: 51).}
}}
}
\value{
A list of the following:
\item{h}{A matrix holding the warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
\item{hInv}{A matrix holding the inverse warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
\item{qfAligned}{A list of matrices. The \eqn{i}-th matrix holds the \eqn{i}-th aligned 
distributional trajectory on \code{workGrid} represented by quantile functions on \code{qSup}, 
where each row corresponds to one time point in \code{workGrid}.}
\item{qf}{Local Fréchet regression estimates of distributional trajectories on \code{workGrid} 
represented by quantile functions on \code{qSup}, with the same format as \code{qfAligned}.}
\item{qSup}{An increasing grid between 0 and 1 holding the support grid on which the quantile functions are evaluated.}
\item{workGrid}{A vector of length \code{optns$ngrid}. 
An equidistant grid on which (inverse) warping functions in \code{h} and \code{hInv} 
and distributional trajectories in \code{qf} and \code{qfAligned} are evaluated.}
\item{optns}{Control options used.}
\item{costs}{The mean cost associated with each trajectory.}
\item{timingWarp}{The time cost of time warping.}
\item{timingPrsm}{The time cost of local Fréchet regression.}
}
\description{
Pairwise warping for distribution-valued random processes 
\eqn{\{X_i(\cdot)\}_{i=1}^{n}} endowed with 2-Wasserstein metric, 
where each \eqn{X_i(\cdot)} is a distribution-valued random process 
estimated via local Fréchet regression.
}
\references{
\cite{Chen, Y., and Müller, H.-G. (2023). "Uniform convergence of local Fréchet regression with applications to locating extrema and time warping for metric space valued trajectories." The Annals of Statistics, 50(3), 1573--1592.}

\cite{Tang, R. and Müller, H.-G. (2008). "Pairwise curve synchronization for functional data." Biometrika 95, 875--889.}

\cite{Petersen, A., and Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
}
