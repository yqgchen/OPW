#' @title Wasserstein pairwise warping for distribution-valued random processes.
#' @description Pairwise warping for distribution-valued random processes 
#' \eqn{\{X_i(\cdot)\}_{i=1}^{n}} endowed with 2-Wasserstein metric, 
#' where each \eqn{X_i(\cdot)} is a distribution-valued random process.
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
#' as well as 
#' \describe{
#' \item{ngrid}{Number of grid points for evaluating the warping functions (default: 51).}
#' }
#' @return A list of the following:
#' \item{h}{A matrix holding the warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{hInv}{A matrix holding the inversed warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{qfAligned}{A list of matrices. The \eqn{i}-th matrix holds the \eqn{i}-th aligned 
#' distributional trajectory on \code{workGrid} represented by quantile functions on \code{qSup}, 
#' where each row corresponds to one time point in \code{workGrid}.}
#' \item{qf}{Local Fréchet regression estimates of distributional trajectories on \code{workGrid} 
#' represented by quantile functions on \code{qSup}, with the same format as \code{qfAligned}.}
#' \item{qSup}{A increasing grid between 0 and 1.}
#' \item{workGrid}{A vector of length \code{optns$ngrid}. 
#' An equidistant grid on which warping functions and local Fréchet regression estimated trajectories are evaluated.}
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
  if ( is.null(optns$nknots) ) optns$nknots <- 2
  if ( is.null(optns$subsetProp) ) optns$subsetProp <- 1
  if (is.null(optns$choice)) optns$choice = 'truncated'
  if ( !(optns$choice %in% c('truncated','weighted') )) {
    stop("The estimation of warping functions can only be done by 'truncated' or 'weighted' average.")
  }
  if (is.null(optns$isPWL)) optns$isPWL <- TRUE
  if (is.null(optns$seed)) optns$seed <- 666
  if (is.null(optns$verbose)) optns$verbose <- FALSE
  if (is.null(optns$kernelReg)) optns$kernelReg <- 'epan'
  # lambda is to be set up later
  
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
    t(qf[[i]]$qout) # each column corresponds to a time point
  } )
  timingPrsm <- Sys.time() - timingPrsm
  
  numOfKcurves = min(round(optns$subsetProp * (nsubj-1)))
  gijMat <- array(dim = c(numOfKcurves,M,nsubj) ) 
  distMat <- matrix( nrow = nsubj, ncol = numOfKcurves)
  hMat <- array(dim = c(nsubj,M) )
  hInvMat <- array(dim = c(nsubj,M) )
  qfAligned <- list()
  
  ## extract the closest subset to tj of workGrid (standardized presmoothing time grid)
  get_tJ <- function(tj=workGrid){
    workGrid[round(tj * (M-1)) + 1]
  }
  
  ## extract the quantile functions for the jth country on the closest subset to tj of workGrid (standardized presmoothing time grid)
  getQtJ <- function(j, tj = workGrid){
    qf[[j]][, round(tj * (M-1)) + 1]
  }
  
  ## evaluating the warping function on tGrid according to res 
  getSol <- function(res, tGrid){
    approx(x = seq(0,1, length.out = (2 + optns$nknots)), y = c(0, sort(res), 1), xout = tGrid)$y
    #approx(x = seq(0,1, length.out = (2+ optns$nknots)), y = c(0, Rcppsort(res),1) ,n = M)$y
    #RcppPseudoApprox(X = seq(0,1, length.out = (2+ optns$nknots)), Y = c(0, Rcppsort(res),1), X_target = seq(0,1, length.out = M))
  }
  
  theCostOptim <- function(x, i, j, lambda, ti){
    tj = getSol(x, ti)
    pracma::trapz(ti, apply((getQtJ(j, tj) - getQtJ(i, ti))^2, 2, sum) + lambda * (get_tJ(tj)-ti)^2)
    #fdapace::trapzRcpp(ti, apply((getQtJ(tj, qtj) - qti)^2, 2, sum) + lambda * (tj-ti)^2)
  }
  
  getGijOptim <- function(i, j, lambda, minqaAvail ){
    s0 <- seq(0,1,length.out = (2+ optns$nknots))[2:(1+optns$nknots)]
    if( !minqaAvail ) { 
      optimRes <- optim( par = s0, fn = theCostOptim, method = 'L-BFGS-B', 
                         lower = rep(1e-6, optns$nknots), upper = rep(1 - 1e-6, optns$nknots),
                         i = i, j = j, lambda = lambda, ti = workGrid)
    } else {
      optimRes <-  minqa::bobyqa( par = s0, fn = theCostOptim,  
                                  lower = rep(1e-6, optns$nknots), upper = rep(1 - 1e-6, optns$nknots),
                                  i = i, j = j, lambda = lambda, ti = workGrid)
    }
    bestSol <- getSol(optimRes$par, seq(0,1,length.out=M))  
    return( bestSol )
  }
  
  if (is.null(optns$lambda)){
    Vy = sqrt( sum( unlist(lapply(seq_along(qf), function(i) {
      qti = getQtJ(i)
      mu = rowMeans(qti)
      pracma::trapz(workGrid, apply((qti - mu)^2, 2, sum)) * diff(trange)
      #fdapace::trapzRcpp(workGrid, apply((qti - mu)^2, 2, sum)) * diff(trange)
    }) ) )/(nsubj-1) )
    optns$lambda <- Vy*10^-4
  } 
  lambda <- optns$lambda
  
  if( !is.element('minqa', installed.packages()[,1]) && optns$isPWL){
    warning("Cannot use 'minqa::bobyqa' to find the optimal knot locations as 'minqa' is not installed. We will do an 'L-BFGS-B' search.") 
    minqaAvail = FALSE
  } else {
    minqaAvail = TRUE
  }
  
  timingWarp <- Sys.time()
  
  for(i in seq_len(nsubj)){ # For each trajectory
    if(optns$verbose){
      cat('Computing pairwise warping for trajectory #', i, 'out of', nsubj, 'trajectories.\n')
    }
    set.seed( i + optns$seed )
    candidateKcurves = sample(seq_len(nsubj)[-i], numOfKcurves)  
    
    for(j in seq_len(numOfKcurves)){ # For each of the candidate curves 
      gijMat[j, ,i] = getGijOptim(i, candidateKcurves[j], lambda, minqaAvail)
      distMat[i,j] = sqrt(pracma::trapz(workGrid, apply((getQtJ(j=candidateKcurves[j], tj = gijMat[j, ,i]) - getQtJ(i))^2, 2, sum)))
      #fdapace::trapzRcpp(workGrid, apply((getQtJ(tj = gijMat[j, ,i], qtj = qf[[candidateKcurves[j]]]) - qti)^2, 2, sum))
    }
    
    if(optns$choice == 'weighted'){
      hInvMat[i,] = apply(gijMat[, ,i] , 2, weighted.mean, 1/distMat[i,])
    } else {
      hInvMat[i,] = apply(gijMat[  (distMat[i,] <= quantile( distMat[i,], p=0.90) ),  ,i] , 2, mean)
    }
    
    hMat[i,] = approx(y = workGrid, x = hInvMat[i,], xout = workGrid)$y
    #alignedMat[[i]] = getQtJ(hMat[i,],qti)
    
    ## quantile function at h(t) with t in StartEnd$start:StartEnd$end
    qfAligned[[i]] = getQtJ(j = i, tj = hMat[i,])
  }
  
  timingWarp <- Sys.time() - timingWarp
  
  qf <- lapply(qf, t)
  qfAligned <- lapply(qfAligned, t)
  
  return(list(
    h = hMat * diff(trange) + trange[1], 
    hInv = hInvMat * diff(trange) + trange[1],
    qfAligned = qfAligned,
    qf = qf, 
    qSup = qSup,
    workGrid = workGrid * diff(trange) + trange[1],
    optns = optns,
    costs = rowMeans(distMat),
    timingWarp = timingWarp,
    timingPrsm = timingPrsm
  ))
}
