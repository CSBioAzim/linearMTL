TGGLMixtureModelSelection <- function(X = NULL, task.specific.features = list(), Y,
                                      M.vec, lambda.vec, restarts = 1,
                                      groups, weights, gam = 0,
                                      num.threads = 1, verbose = -1, ...) {
  ##################
  # error checking #
  ##################

  # check if any input data was supplied
  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # check dimensions of shared features
  J1 <- 0
  if (!is.null(X)) {
    if (nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows!")
    }
    J1 <- ncol(X)
  }

  # check dimensions of task specific features
  J2 <- 0
  if (length(task.specific.features) > 0) {
    for (task.matrix in task.specific.features) {
      if (nrow(task.matrix) != nrow(Y)) {
        stop("Task specific feature matrices and Y must have the same number of rows!")
      }
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  # check dimensions of groups and weights
  if (ncol(Y) != ncol(groups)) {
    stop("Y and groups must have the same number of columns!")
  }
  if (length(weights) != nrow(groups)) {
    stop("Length of weights has to equal number of rows of groups!")
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  # build parameter grid
  parameter.grid <- cbind(rep(M.vec, each = length(lambda.vec)),
                          rep(lambda.vec, times = length(M.vec)))

  RunParameter <- function (ind) {
    M <- parameter.grid[ind, 1]
    lambda <- parameter.grid[ind, 2]

    top.model <- NULL
    top.lik <- -Inf
    for (r in 1:restarts) {
      tggl.mix <- TGGLMixture(X = X,
                              task.specific.features = task.specific.features,
                              Y = Y, M = M, groups = groups, weights = weights,
                              lambda = lambda, gam = gam)
      if (verbose > 0) {
        print(sprintf('Param [M = %d, lambda = %.e]: LL: %.3f', M, lambda, tggl.mix$likelihood))
      }
      if (tggl.mix$likelihood > top.lik) {
        top.model <- tggl.mix
        top.lik <- tggl.mix$likelihood
      }
    }

    tggl.mix <- top.model

    # compute BIC
    tggl.param <- 0
    for (m in 1:M) {
      model.coef <- tggl.mix$models[[m]]
      model.coef <- model.coef[-1, ]
      for (v in 1:nrow(groups)) {
        group <- as.logical(groups[v, ])
        group.size <- sum(group)
        subcoef.norm <- sqrt(rowSums(model.coef[, group, drop = FALSE]^2))
        subcoef.norm <- subcoef.norm[subcoef.norm != 0]
        tggl.param <- tggl.param + sum(1 + (group.size - 1)/(1 + weights[v] * tggl.mix$prior[m]^gam * lambda/subcoef.norm))
      }
    }
    # tggl.coef <- do.call(cbind, tggl.mix$models)
    # tggl.param <- sum(tggl.coef[-1, ] != 0)

    # parameters are:
    # - M*K variance parameters
    # - M*K intercept parameters
    # - M-1 prior probabilities
    # - parameters of TGGL (estimated as nonzero coefficients)

    num.of.parameters <- 2*M*K + (M-1) + tggl.param
    BIC = -2*tggl.mix$likelihood + log(N)*num.of.parameters
    if (verbose > -1) {
      print(sprintf('Param [M = %d, lambda = %.e]: BIC: %.3f, Num of Param. : %.2f, L: %3f.',
                    M, lambda, BIC, num.of.parameters, tggl.mix$likelihood))
    }
    return(list(model = top.model, BIC = BIC))
  }

  if (verbose > -2) {
    print(sprintf("Running TGGL-Mixture for %d parameter settings ... ", nrow(parameter.grid)))
  }
  sel.start.time <- Sys.time()

  doMC::registerDoMC(num.threads)
  param.inds <- 1:nrow(parameter.grid)
  bic.results <- foreach(l = param.inds) %dopar% RunParameter(l)
  opt.model.idx <- which.min(sapply(bic.results, FUN = function(x){x$BIC}))

  sel.end.time <- Sys.time()
  sel.time <- as.numeric(sel.end.time - sel.start.time, units = "mins")
  if (verbose > -2) {
    print(sprintf("Minutes to run Model Selection %0.1f. Best Parameter Setting: M = %d, lambda = %.e",
                  sel.time, parameter.grid[opt.model.idx, 1], parameter.grid[opt.model.idx, 2]))
  }
  return(bic.results[[opt.model.idx]]$model)
}
