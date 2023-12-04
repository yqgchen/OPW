#' @title Wasserstein pairwise warping for distribution-valued random processes
#' @description Pairwise warping for distribution-valued random processes 
#' \eqn{\{X_i(\cdot)\}_{i=1}^{n}} endowed with 2-Wasserstein metric, 
#' where each \eqn{X_i(\cdot)} is a distribution-valued random process 
#' estimated via local Fréchet regression.
#' @param Lt A list of \eqn{n} vectors, where \code{Lt[[i]]} is a vector of 
#' length \eqn{m_i} holding the observation time points for the \eqn{i}-th 
#' distribution-valued process \eqn{X_i(\cdot)}. Each vector should be sorted in ascending order.
#' @param Ly A list of \eqn{n} lists, where \code{Ly[[i]]} is a list of 
#' \eqn{m_i} vectors and the \eqn{j}-th vector \code{Ly[[i]][[j]]} holds the sample of data points 
#' generated from (a noisy version of) the distribution \eqn{X_i(t_{ij})} with \eqn{t_{ij}} given by 
#' \code{Lt[[i]][j]}. 
#' Note that only one of the three, \code{Ly}, \code{Lh}, and \code{Lq}, needs to be input. 
#' If more than one of them are specified, \code{Ly} overwrites \code{Lh}, and \code{Lh} overwrites \code{Lq}. 
#' @param Lh A list of \eqn{n} lists, where \code{Lh[[i]]} is a list of \eqn{m_i} lists 
#' and \code{Lh[[i]][[j]]} holds the histogram corresponding to (a noisy version of) 
#' the distribution \eqn{X_i(t_{ij})} with \eqn{t_{ij}} given by \code{Lt[[i]][j]}. 
#' Specifically, \code{Lh[[i]][[j]]} is a list of two components, 
#' \code{breaks} (or \code{mids}) holding the cell boundaries (or cell midpoints) and 
#' \code{counts} holding the counts falling in each cell. 
#' @param Lq A list of \eqn{n} lists, where \code{Lq[[i]]} is a matrix or list 
#' holding the quantile functions corresponding to the \eqn{i}-th distribution-valued process. 
#' If \code{Lq[[i]]} is a matrix, the support of the \eqn{\sum_{i=1}^{n} m_i} quantile functions 
#' should be the same (i.e., \code{optns$qSup}), and \eqn{j}-th row of \code{Lq[[i]]} holds 
#' the quantile function corresponding to corresponding to (a noisy version of) 
#' the distribution \eqn{\{X_i(t_{ij})\}_{j=1}^{m_i}} with \eqn{t_{ij}} given by \code{Lt[[i]][j]}. 
#' If the quantile functions are evaluated on different grids, then \code{Lq[[i]]} should be a list, 
#' with \code{Lq[[i]][[j]]} being a list of two components \code{x} and \code{y} 
#' holding the support grid and the corresponding values of the quantile functions, respectively. 
#' @param optns A list of options control parameters specified by \code{list(name=value)}, 
#' including all those available in \code{\link[frechet]{LocDenReg}} and \code{\link[fdapace]{WFDA}}, 
#' except for that the default \code{subsetProp} is 1 and default \code{kernelReg} is \code{'epan'}; 
#' as well as 
#' \describe{
#' \item{ngrid}{Number of grid points for evaluating the warping functions (default: 51).}
#' }
#' @return A list of the following:
#' \item{h}{A matrix holding the warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{hInv}{A matrix holding the inverse warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{qfAligned}{A list of matrices. The \eqn{i}-th matrix holds the \eqn{i}-th aligned 
#' distributional trajectory on \code{workGrid} represented by quantile functions on \code{qSup}, 
#' where each row corresponds to one time point in \code{workGrid}.}
#' \item{qf}{Local Fréchet regression estimates of distributional trajectories on \code{workGrid} 
#' represented by quantile functions on \code{qSup}, with the same format as \code{qfAligned}.}
#' \item{qSup}{An increasing grid between 0 and 1 holding the support grid on which the quantile functions are evaluated.}
#' \item{workGrid}{A vector of length \code{optns$ngrid}. 
#' An equidistant grid on which (inverse) warping functions in \code{h} and \code{hInv} 
#' and distributional trajectories in \code{qf} and \code{qfAligned} are evaluated.}
#' \item{optns}{Control options used.}
#' \item{costs}{The mean cost associated with each trajectory.}
#' \item{timingWarp}{The time cost of time warping.}
#' \item{timingPrsm}{The time cost of local Fréchet regression.}
#' @references
#' \cite{Chen, Y., and Müller, H.-G. (2023). "Uniform convergence of local Fréchet regression with applications to locating extrema and time warping for metric space valued trajectories." The Annals of Statistics, 50(3), 1573--1592.}
#' 
#' \cite{Tang, R. and Müller, H.-G. (2008). "Pairwise curve synchronization for functional data." Biometrika 95, 875--889.}
#' 
#' \cite{Petersen, A., and Müller, H.-G. (2019). "Fréchet regression for random objects with Euclidean predictors." The Annals of Statistics, 47(2), 691--719.}
#' 
#' @importFrom frechet LocDenReg
#' @importFrom pracma trapz
#' @importFrom stats approx optim quantile weighted.mean
#' @importFrom utils installed.packages
#' @export
#' 
WassPW <- function ( Lt=NULL, Ly=NULL, Lh=NULL, Lq=NULL, optns = list() ){
  
  # check input
  if ( !is.list (Lt) | any( !(sapply(Lt, is.numeric) & sapply(Lt, is.vector)) ) ) {
    stop("Missing input or incorrect format of Lt.")
  }
  nsubj <- length(Lt)
  ntime_per_subj <- sapply( Lt, length )
  
  if ( !is.null(Ly) ) {
    if ( !is.list(Ly) ) {
      stop( "Ly should be a list." )
    } else if ( any( !sapply(Ly, is.list) ) ) {
      stop( "Ly[[i]] should be a list for all plausible i.")
    } else if ( any( abs( sapply( Ly, length ) - ntime_per_subj ) > 0 ) ) {
      stop( "Mismatched numbers of time points per subject between Lt and Ly." )
    } else if ( length(Ly) != nsubj ) {
      stop( "Mismatched numbers of trajectories/subjects between Lt and Ly. ")
    } else if ( !all( unlist( lapply( Ly, sapply, is.vector ) ) ) | !all( unlist( lapply( Ly, sapply, is.numeric ) ) ) ) {
      stop("Incorrect format of Ly.")
    }
  } else if ( !is.null(Lh) ) {
    if ( !is.list(Lh) ) {
      stop( "Lh should be a list." )
    } else if ( any( !sapply(Lh, is.list) ) ) {
      stop( "Lh[[i]] should be a list for all plausible i.")
    } else if ( any( abs( sapply( Lh, length ) - ntime_per_subj ) > 0 ) ) {
      stop( "Mismatched numbers of time points per subject between Lt and Lh." )
    } else if ( length(Lh) != nsubj ) {
      stop( "Mismatched numbers of trajectories/subjects between Lt and Lh. ")
    } else if ( !all( unlist( lapply( Lh, sapply, is.list ) ) ) ) {
      stop("Incorrect format of Lh.")
    } else if ( any( unlist( lapply( Lh, sapply, function (h) {
      ( is.null(h$breaks) & is.null(h$mids) ) | is.null(h$counts)
    }) ) ) ) {
      stop("Incorrect format of Lh.")
    }
  } else if ( !is.null(Lq) ) {
    if ( !is.list(Lq) ) {
      stop( "Lq should be a list." )
    } else if ( any( !sapply(Lq, is.list) & !sapply(Lq, is.matrix) ) ) {
      stop( "Lq[[i]] should be a list or matrix for all plausible i.")
    } else if ( abs ( length(Lq) - nsubj ) > 0 ) {
      stop( "Mismatched numbers of trajectories/subjects between Lt and Lq. ")
    } else {
      for ( i in seq_len(nsubj) ) {
        if ( is.matrix(Lq[[i]]) ) {
          if ( abs( nrow(Lq[[i]]) - ntime_per_subj[i] ) > 0 ) {
            stop( "Mismatched numbers of time points per subject between Lt and Lq." )
          } else if ( any( apply( Lq[[i]], 1, is.unsorted ) ) ) {
            stop( "Quantile functions in Lq are not monotonic.")
          }
        } else {
          # Lq[[i]] is a list of lists
          if ( abs( length(Lq[[i]]) - ntime_per_subj[i] ) > 0 ) {
            stop( "Mismatched numbers of time points per subject between Lt and Lh." )
          } else if ( !all(sapply(Lq[[i]],is.list)) ) {
            stop("Incorrect format of Lq.")
          } else if ( any( sapply( Lq[[i]], function(q) {
            is.null(q$x) | is.null(q$y)
          }) ) ) {
            stop("Incorrect format of Lq.")
          }
        }
      }
    }
    
  } else {
    stop("At least one of Ly, Lh and Lq should be input.")
  }
  
  # set up some default options
  if ( is.null(optns$ngrid) ) optns$ngrid <- 51
  if (is.null(optns$kernelReg)) optns$kernelReg <- 'epan'
  
  tin <- unique( sort( unlist( Lt ) ) )
  trange <- range(tin)
  tin <- lapply( Lt, function(t) ( t - trange[1] ) / diff(trange) )
  M <- optns$ngrid
  workGrid = seq( 0, 1, length.out = M ) # a grid on [0,1]
  
  ## pre-smoothing
  timingPrsm <- Sys.time()
  if ( !is.null(Ly) ) {
    qf <- lapply(1:nsubj, function(i) {
      LocDenReg( xin = tin[[i]], yin = Ly[[i]], xout = workGrid, optns = optns )
    })
  } else if ( !is.null(Lh) ) {
    qf <- lapply(1:nsubj, function(i) {
      LocDenReg( xin = tin[[i]], hin = Lh[[i]], xout = workGrid, optns = optns )
    })
  } else {
    qf <- lapply(1:nsubj, function(i) {
      LocDenReg( xin = tin[[i]], qin = Lq[[i]], xout = workGrid, optns = optns )
    })
  }
  qSup <- qf[[1]]$qSup
  qf <- lapply( 1:nsubj, function(i) { 
    qf[[i]]$qout # each column corresponds to a time point
  } )
  timingPrsm <- Sys.time() - timingPrsm
  
  # time warping
  
  res <- WassPWdense( tVec = workGrid, qf = qf, qSup = qSup, optns = optns )
  
  return(list(
    h = res$h * diff(trange) + trange[1], 
    hInv = res$hInv * diff(trange) + trange[1],
    qfAligned = res$qfAligned,
    qf = qf, 
    qSup = qSup,
    workGrid = res$workGrid * diff(trange) + trange[1],
    optns = res$optns,
    costs = res$costs,
    timingWarp = res$timingWarp,
    timingPrsm = timingPrsm
  ))
}
