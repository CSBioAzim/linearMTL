#' Fit a latent variable linar multi-task model (Zhang et al. 2008).
#'
#' A Bayesian linear multi-task model where the regression matrix is assumed to
#' be composed of latent factors.
#'
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param S H by K loading matrix.
#' @param (Optional) max.iter Maximum number of iterations.
#' @param (Optional) epsilon Desired accuracy. If error change drops below
#'   epsilon, the algorithm terminates.
#' @param XTX (Optional) Precomputed matrices t(X)*X as for example produced by
#'   PrepareMatrices.
#' @param XTY (Optional) Precomputed matrices t(X)*Y as for example produced by
#'   PrepareMatrices
#' @param verbose (Optional) Integer in {0,1,2}. verbose = 0: No output. verbose
#'   = 1: Print summary at the end of the optimization. verbose = 2: Print
#'   progress during optimization.
#'
#' @return List containing
#'   \item{Gamma}{Estimated mixing matrix.}
#'   \item{sigma2}{Estimated sigma^2.}
#'   \item{Psi}{Estimated Psi.}
#'   \item{S}{Loading matrix used.}
#'   \item{B}{MAP estimate of the regression coefficients.}
#'   \item{early.termination}{Boolean indicating whether the algorithm exceeded
#'   max.iter iterations.}
#'
#' @seealso \code{\link{RunGroupCrossvalidation}}
#' @export
LatentVariableRegression <- function (task.specific.features = list(), Y, S,
                                  max.iter = 10000, epsilon = 1e-5,
                                  XTX = NULL, XTY = NULL,
                                  verbose = 1) {
  if (!is.null(XTX) | !is.null(XTY)) {
    if (xor(!is.null(XTX), !is.null(XTY))) {
      stop("When using immutables, both XTX and XTY need to be specified!")
    }
    use.cached.immutables <- TRUE
  } else {
    use.cached.immutables <- FALSE
  }

  # precompute matrices if necessary
  if (!use.cached.immutables) {
    mats <- PrepareMatrices(Y = Y, task.specific.features = task.specific.features,
                            standardize = FALSE)
    XTX <- mats$XTX
    XTY <- mats$XTY
  }

  H <- nrow(S)
  K <- ncol(Y)
  N <- nrow(Y)
  J <- ncol(task.specific.features[[1]])

  Gamma <- matrix(0, nrow = J, ncol = H)
  GammaS <- Gamma %*% S
  Psi <- diag(J)
  sigma2 <- 100

  SST.inverse <- solve(S %*% t(S))

  delta <- Inf
  iter <- 0
  gradient <- 0
  start.time <- Sys.time()
  while (delta > epsilon && iter < max.iter) {
    Gamma.old <- Gamma
    Psi.old <- Psi
    sigma2.old <- sigma2

    Psi.inverse <- solve(Psi)

    V <- list()
    M <- matrix(0, nrow = J, ncol = K)
    XTXV.traces <- 0
    # E-step
    for (k in 1:K) {
      V[[k]] <- solve(Psi.inverse + 1/sigma2 * XTX[[k]])
      M[, k] <- V[[k]] %*% (Psi.inverse %*% GammaS[, k] + 1/sigma2 * XTY[, k])
      XTXV.traces <- XTXV.traces + sum(diag(XTX[[k]] %*% V[[k]]))
    }
    # M-step
    Gamma <- M %*% t(S) %*% SST.inverse
    GammaS <- Gamma %*% S
    V.sum <- Reduce("+", V)
    Psi <- 1/K * V.sum + 1/K * (M - GammaS) %*% t(M - GammaS)
    sigma2 <- MTComputeError(LMTL.model = list(B = M, intercept = rep(0, ncol(M))),
                             Y = Y, task.specific.features = task.specific.features,
                             normalize = TRUE) + XTXV.traces / (N*K)

    diff.gamma <- sqrt(sum((Gamma.old - Gamma)^2))/sqrt(sum(Gamma.old^2))
    diff.psi <- sqrt(sum((Psi.old - Psi)^2))/sqrt(sum(Psi.old^2))
    diff.sigma2 <- abs(sigma2.old - sigma2)/sigma2.old
    delta <- max(diff.gamma, diff.psi, diff.sigma2)
    iter <- iter + 1
    if (verbose == 2) {
      print(sprintf("Iter: %d, Max Relative Parameter Change: %f.", iter, delta))
    }
  }
  end.time <- Sys.time()
  if (verbose >= 1) {
    print(sprintf("Total: Iter: %d. Time: %f", iter, as.numeric(end.time - start.time, units = "secs")))
  }

  B <- matrix(0, nrow = J, ncol = K)
  Psi.inverse <- solve(Psi)
  for (k in 1:K) {
    B[, k] <- solve(Psi.inverse + 1/sigma2 * XTX[[k]]) %*% (Psi.inverse %*% GammaS[, k] + 1/sigma2 * XTY[, k])
  }

  return(list(Gamma = Gamma,
              sigma2 = sigma2,
              Psi = Psi,
              S = S,
              B = B,
              early.termination = iter > max.iter))
}
