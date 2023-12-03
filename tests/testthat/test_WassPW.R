library(testthat)
# test for WassPW()

nsubj <- 20
ntime <- 21
qSup <- seq(0,1,0.05)
Lt <- rep( list(seq_len(ntime)), nsubj )
Ly <- lapply(seq_len(nsubj), function (i){
  lapply(seq_len(ntime), function (j) {
    rnorm(51, mean = j)
  })
})
Lh <- lapply(seq_len(nsubj), function (i){
  lapply(seq_len(ntime), function (j) {
    hist(Ly[[i]][[j]])
  })
})
Lq <- lapply(seq_len(nsubj), function (i){
  t(sapply(seq_len(ntime), function (j) {
    qbeta(qSup, j, 2)
  }))
})

expect_error( res1 <- WassPW(Lt = Lt, Ly = Ly), NA )
expect_error( res2 <- WassPW(Lt = Lt, Lh = Lh), NA )
expect_error( res3 <- WassPW(Lt = Lt, Lq = Lq), NA )
expect_lt( sum( abs( t(res1$h) - seq(1,ntime,length.out = ncol(res1$h)) ) ), 1e-8 )
expect_lt( sum( abs( t(res2$h) - seq(1,ntime,length.out = ncol(res2$h)) ) ), 1e-8 )
expect_lt( sum( abs( t(res3$h) - seq(1,ntime,length.out = ncol(res3$h)) ) ), 1e-8 )
