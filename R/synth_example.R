#' Test of tree-guided group lasso on synthetic data.
#'
#' Generate synthetic data as in Kim et al. 2010 and apply tree-guided group
#' lasso and lasso models.
#'
#' @param method Method to use. "group" (TGGL), "tbt" (task-by-task), "all" (both).
#'
#' @importFrom lattice levelplot
#' @importFrom grDevices gray
#' @export
TestLinearMTL <- function(method = "group") {
  set.seed(12345)
  # number of data points
  N <- 150
  # number of input dimensions
  J <- 200
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
  Xtrain <- matrix(sample(c(0,1,2), N * J, replace = T), nrow = N)
  Xtest <- matrix(sample(c(0,1,2), N * J, replace = T), nrow = N)
  Ytrain <- Xtrain %*% B.truth + matrix(rnorm(N * K), nrow = N)
  Ytest <- Xtest %*% B.truth + matrix(rnorm(N * K), nrow = N)

  if ((method == "group") | (method == "all")) {
    M <- BuildTreeHC(Ytrain, 0.6)
    lambda.vec = 10^seq(-3, 2, length.out = 6)
    tggl <- RunGroupCrossvalidation(X = Xtrain,
                                    Y = Ytrain,
                                    groups = M$groups,
                                    weights.matrix = matrix(M$weights, nrow = 1),
                                    lambda.vec = lambda.vec,
                                    epsilon = 1e-4,
                                    mu = 1e-5, mu.adapt = 0.99,
                                    verbose = 0)
    tggl.model <- tggl$full.model
  }
  if ((method == "tbt") | (method == "all")) {
    lambda.vec = 10^seq(-5, 2, length.out = 20)
    tbt <- RunTBTCrossvalidation(X = Xtrain,
                                 Y = Ytrain,
                                 lambda.vec = lambda.vec)
    tbt.model <- tbt$full.model
  }

  # evaluate
  # compute test error and test correlation
  if ((method == "group") | (method == "all")) {
    group.test.pred <- MTPredict(LMTL.model = tggl.model, X = Xtest)
    group.err <- MTComputeError(LMTL.model = tggl.model, Y = Ytest, X = Xtest)
    group.cor <- MTComputeMeanCorrelation(LMTL.model = tggl.model, Y = Ytest, X = Xtest)
    print(sprintf("GROUP  - test mse: %f, test cor: %f", group.err, group.cor))
  }
  if ((method == "tbt") | (method == "all")) {
    sbs.test.pred <- MTPredict(LMTL.model = tbt.model, X = Xtest)
    sbs.err <- MTComputeError(LMTL.model = tbt.model, Y = Ytest, X = Xtest)
    sbs.cor <- MTComputeMeanCorrelation(LMTL.model = tbt.model, Y = Ytest, X = Xtest)
    print(sprintf("SBS    - test mse: %f, test cor: %f", sbs.err, sbs.cor))
  }
  opt.model <- list(B = B.truth, intercept = rep(0, ncol(B.truth)))
  opt.pred <- MTPredict(LMTL.model = opt.model, X = Xtest)
  opt.err <- MTComputeError(LMTL.model = opt.model, Y = Ytest, X = Xtest)
  opt.cor <- MTComputeMeanCorrelation(LMTL.model = opt.model, Y = Ytest, X = Xtest)
  print(sprintf("OPT    - test mse: %f, test cor: %f", opt.err, opt.cor))

  # plot coefficients
  scl <- seq(-0.3,0.5, by = 0.01)
  colreg <- grDevices::gray(1:(length(scl)+1)/(length(scl)+1))
  if ((method == "group") | (method == "all")) {
    print(lattice::levelplot(tggl.model$B, at = scl, col.regions = colreg))
  }
  if ((method == "tbt") | (method == "all")) {
    print(lattice::levelplot(tbt.model$B, at = scl, col.regions = colreg))
  }
  print(lattice::levelplot(B.truth, at = scl, col.regions = colreg))

  # plot error
  if (method == "all") {
    plot(colSums((sbs.test.pred - Ytest)^2) / N, col = "red", type = "p")
    lines(colSums((group.test.pred - Ytest)^2) / N, type = "p", col = "blue")
    lines(colSums((Xtest %*% B.truth - Ytest)^2) / N, type = "p", col = "black")
  }
}
