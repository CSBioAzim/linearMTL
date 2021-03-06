#' Test of tree-guided group lasso on synthetic data.
#'
#' Generate synthetic data as in Kim et al. 2010 and apply tree-guided group
#' lasso and lasso models.
#'
#' @export
#' @importFrom lattice levelplot
#' @importFrom grDevices gray
#' @importFrom stats rnorm
#' @importFrom RColorBrewer brewer.pal
#' @import latticeExtra
TestLinearMTL <- function() {
  set.seed(12345)
  # number of data points
  N <- 50
  # number of input dimensions
  J <- 150
  # number of tasks
  K <- 60
  # signal intensity
  signal <- 0.6

  # generate hierarchical pattern in B.truth
  setsignal <- function(start, thickness, steps, B.truth) {
    d <- ncol(B.truth)
    l <- d/steps
    offset <- start
    for (step in 0:(steps-1)) {
      for (i in 1:thickness) {
        B.truth[offset, (step * l + 1):((step + 1) * l)] <- signal
        offset <- offset + 1
      }
    }
    return(B.truth)
  }
  B.truth <- matrix(0, nrow = J, ncol = K)
  B.truth <- setsignal(start = 5, thickness = 1, steps = 3, B.truth)
  B.truth <- setsignal(start = 10, thickness = 2, steps = 6, B.truth)
  B.truth <- setsignal(start = 30, thickness = 3, steps = 12, B.truth)
  B.truth <- setsignal(start = 75, thickness = 1, steps = 60, B.truth)
  B.truth <- B.truth[, ncol(B.truth):1]
  B.truth[2, ] <- signal

  # generate training and test matrices
  Xtrain <- matrix(rnorm(N * J), nrow = N)
  Ytrain <- Xtrain %*% B.truth + matrix(rnorm(N * K), nrow = N)
  Xtest <- matrix(rnorm(200 * J), nrow = 200)
  Ytest <- Xtest %*% B.truth + matrix(rnorm(200 * K), nrow = 200)

  M <- BuildTreeHC(B.truth, 2)
  lambda.vec = 10^seq(-2, 1, length.out = 20)
  tggl <- RunGroupCrossvalidation(X = Xtrain,
                                  Y = Ytrain,
                                  groups = M$groups,
                                  weights.matrix = matrix(M$weights, nrow = 1),
                                  lambda.vec = lambda.vec,
                                  epsilon = 1e-4,
                                  mu = 1e-5, mu.adapt = 1,
                                  verbose = 0, standardize = T)
  tggl.model <- tggl$full.model

  lambda.vec = 10^seq(-2, 1, length.out = 20)
  tbt <- RunTBTCrossvalidation(X = Xtrain,
                               Y = Ytrain,
                               lambda.vec = lambda.vec,
                               standardize = T)
  tbt.model <- tbt$full.model


  # evaluate
  # compute test error and test correlation
  sbs.test.pred <- MTPredict(LMTL.model = tbt.model, X = Xtest)
  sbs.err <- MTComputeError(LMTL.model = tbt.model, Y = Ytest, X = Xtest)
  sbs.cor <- MTComputeMeanCorrelation(LMTL.model = tbt.model, Y = Ytest, X = Xtest)
  print(sprintf("SBS    - test mse: %f, test cor: %f", sbs.err, sbs.cor))

  group.test.pred <- MTPredict(LMTL.model = tggl.model, X = Xtest)
  group.err <- MTComputeError(LMTL.model = tggl.model, Y = Ytest, X = Xtest)
  group.cor <- MTComputeMeanCorrelation(LMTL.model = tggl.model, Y = Ytest, X = Xtest)
  print(sprintf("GROUP  - test mse: %f, test cor: %f", group.err, group.cor))

  opt.model <- list(B = B.truth, intercept = rep(0, ncol(B.truth)))
  opt.pred <- MTPredict(LMTL.model = opt.model, X = Xtest)
  opt.err <- MTComputeError(LMTL.model = opt.model, Y = Ytest, X = Xtest)
  opt.cor <- MTComputeMeanCorrelation(LMTL.model = opt.model, Y = Ytest, X = Xtest)
  print(sprintf("OPT    - test mse: %f, test cor: %f", opt.err, opt.cor))

  colreg <- c(RColorBrewer::brewer.pal(8, "PuBu")[8:2], 'white',
              RColorBrewer::brewer.pal(8, "OrRd")[2:8])
  scl <- seq(-1,1, length.out = length(colreg)+1)

  p1 <- lattice::levelplot(t(B.truth)[,J:1], at = scl, col.regions = colreg,
                           xlab = list(label="Tasks", ps = 12),
                           ylab = list(label="Features", ps = 12),
                           main = "")
  p2 <- lattice::levelplot(t(tbt.model$B)[,J:1], at = scl, col.regions = colreg)
  p3 <- lattice::levelplot(t(tggl.model$B)[,J:1], at = scl, col.regions = colreg)

  plots <- c("Groundtruth" = p1, "TBT-Lasso" = p2, "TGGL" = p3,
             layout = c(3,1), merge.legends = FALSE, x.same = FALSE, y.same = TRUE)
  print(plots)
}
