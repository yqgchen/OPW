#' @title Wasserstein pairwise warping for completely/densely observed distribution-valued random processes
#' @description Pairwise warping for distribution-valued random processes 
#' \eqn{\{X_i(\cdot)\}_{i=1}^{n}} endowed with 2-Wasserstein metric, 
#' where each \eqn{X_i(\cdot)} is a distribution-valued random process densely observed
#'  on the same grid across \eqn{i} (and hence no local Fréchet regression is involved).
#' @param tVec A vector holding the common equidistant time grid on which all \eqn{\{X_i(\cdot)\}_{i=1}^{n}} are observed.
#' @param qf A list of matrices with \code{qf[[i]]} holding the \eqn{i}-th distributional 
#' trajectory \eqn{X_i(\cdot)} on \code{tVec} represented by quantile functions on \code{qSup}, 
#' of which each row corresponds to a time point in \code{tVec}.
#' @param qSup An increasing grid between 0 and 1 holding the support grid on which the quantile functions are evaluated.
#' @param optns A list of options control parameters specified by \code{list(name=value)}, 
#' including all those available in \code{\link[fdapace]{WFDA}} except for that the default \code{subsetProp} is 1. 
#' A difference is regarding the option \code{choice}: \code{'unweighted'} (default), \code{'weighted'}, or \code{'truncated'}. 
#' If \code{choice} is \code{'unweighted'}, a simple average of pairwise warping functions are computed; otherwise see \code{\link[fdapace]{WFDA}} for details.
#' @return A list of the following:
#' \item{h}{A matrix holding the warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{hInv}{A matrix holding the inverse warping functions evaluated on \code{workGrid}, each row corresponding to a subject.}
#' \item{qfAligned}{A list of matrices. The \eqn{i}-th matrix holds the \eqn{i}-th aligned 
#' distributional trajectory on \code{workGrid} represented by quantile functions on \code{qSup}, 
#' where each row corresponds to one time point in \code{workGrid}.}
#' \item{qSup}{An increasing grid between 0 and 1 holding the support grid on which the quantile functions are evaluated.}
#' \item{workGrid}{A copy of \code{tVec}. An equidistant grid on which 
#' (inverse) warping functions in \code{h} and \code{hInv} 
#' and distributional trajectories in \code{qf} and \code{qfAligned} are evaluated.}
#' \item{optns}{Control options used.}
#' \item{costs}{The mean cost associated with each trajectory.}
#' \item{timingWarp}{The time cost of time warping.}
#' \item{g}{An \eqn{(n-1)\times T\times n} array, where \code{g[j,,i]} holds the pairwise function aligning \eqn{X_j} to \eqn{X_i} evaluated on \code{workGrid}, where \eqn{T} is the length of \code{workGrid}.}
#' @references
#' \cite{Chen, Y., and Müller, H.-G. (2023). "Uniform convergence of local Fréchet regression with applications to locating extrema and time warping for metric space valued trajectories." The Annals of Statistics, 50(3), 1573--1592.}
#' 
#' \cite{Tang, R. and Müller, H.-G. (2008). "Pairwise curve synchronization for functional data." Biometrika 95, 875--889.}
#' 
#' @importFrom pracma trapz
#' @importFrom stats approx optim quantile weighted.mean
#' @importFrom utils installed.packages
#' @export
#' 
WassPWdense <- function (tVec, qf, qSup, optns = list() ) {
  # set up some default options
  if ( is.null(optns$nknots) ) optns$nknots <- 2
  # if ( is.null(optns$subsetProp) ) optns$subsetProp <- 1
  if (is.null(optns$choice)) optns$choice = 'unweighted'
  if ( !(optns$choice %in% c('unweighted','weighted','truncated') )) {
    stop("The estimation of warping functions can only be done by 'unweighted' or 'weighted' average.")
  }
  if (is.null(optns$isPWL)) optns$isPWL <- TRUE
  if (is.null(optns$seed)) optns$seed <- 666
  if (is.null(optns$verbose)) optns$verbose <- FALSE
  # lambda is to be set up later
  
  trange <- range(tVec)
  workGrid <- ( tVec - trange[1] ) / diff(trange) # a grid on [0,1]
  M <- length(tVec)
  nsubj <- length(qf)
  qf <- lapply( qf, t ) # each column corresponds to a time point
  
  numOfKcurves = nsubj-1 # min(round(optns$subsetProp * (nsubj-1)))
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
    candidateKcurves = seq_len(nsubj)[-i]
    
    for(j in seq_len(numOfKcurves)){ # For each of the candidate curves 
      gijMat[j, ,i] = getGijOptim(i, candidateKcurves[j], lambda, minqaAvail)
      distMat[i,j] = sqrt(pracma::trapz(workGrid, apply((getQtJ(j=candidateKcurves[j], tj = gijMat[j, ,i]) - getQtJ(i))^2, 2, sum)))
      #fdapace::trapzRcpp(workGrid, apply((getQtJ(tj = gijMat[j, ,i], qtj = qf[[candidateKcurves[j]]]) - qti)^2, 2, sum))
    }
    
    if(optns$choice == 'weighted'){
      hInvMat[i,] = apply(gijMat[, ,i, drop = FALSE] , 2, weighted.mean, 1/distMat[i,])
    } else if ( optns$choice == 'truncated') {
      hInvMat[i,] = apply(gijMat[  (distMat[i,] <= quantile( distMat[i,], p=0.90) ),  ,i, drop = FALSE] , 2, mean)
    } else {
      hInvMat[i,] = apply(gijMat[ , ,i, drop = FALSE] , 2, mean)
    }
    
    hMat[i,] = approx(y = workGrid, x = hInvMat[i,], xout = workGrid)$y
    #alignedMat[[i]] = getQtJ(hMat[i,],qti)
    
    ## quantile function at h(t) with t in StartEnd$start:StartEnd$end
    qfAligned[[i]] = getQtJ(j = i, tj = hMat[i,])
  }
  
  timingWarp <- Sys.time() - timingWarp
  
  for(i in seq_len(nsubj)){
    for(j in seq_len(numOfKcurves)){
      gijMat[j, ,i] = gijMat[j, ,i] * diff(trange) + trange[1]
    }
  }
  
  return(list(
    h = hMat * diff(trange) + trange[1], 
    hInv = hInvMat * diff(trange) + trange[1],
    qfAligned = lapply( qfAligned, t ),
    qSup = qSup,
    workGrid = workGrid * diff(trange) + trange[1],
    optns = optns,
    costs = rowMeans(distMat),
    timingWarp = timingWarp,
    g = gijMat
  ))
  
}
