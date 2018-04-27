#' RunTGGLMixSelection is a convenience function which runs TGGLMix repeatedly
#' for all possible combinations of parameters in M.vec and lambda.vec.
#'
#' Fit a tree-guided group lasso mixture model. Restart with different random
#' initializations and keep the model with the lowest objective value. Optional:
#' Evaluate best model on test set. By default all data is used for training. If
#' validation.ids is not NULL, exclude corresponding indices from training and
#' use them for validating parameters instead.
#'
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param M.vec Vector with numbers of clusters.
#' @param validation.ids (Optional) Indices of data points to be used for
#'   validation. Needs to be supplied if more than one parameter pair (M,
#'   lambda) is given. Default is to use all data for training.
#' @param num.starts (Optional) Number of starts. Default is 1 (no restarts).
#' @param num.threads (Optional) Number of threads to be used. Default is 1.
#' @param groups Binary V by K matrix determining group membership: Task k in
#'   group v iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param lambda.vec Vector with regularization parameters.
#' @param verbose (Optional) Integer in {0,1,2,3}. verbose = 0: No output.
#'   verbose = 1: Print final summary. verbose = 2: Print summary for each
#'   parameter. verbose = 3: Print summary for each restart.
#' @param gam (Optional) Regularization parameter for component m will be lambda
#'   times the prior for component m to the power of gam.
#' @param homoscedastic (Optional) Force variance to be the same for all tasks
#'   in a component. Default is FALSE.
#' @param EM.max.iter (Optional) Maximum number of iterations for EM algorithm.
#' @param EM.epsilon (Optional) Desired accuracy. Algorithm will terminate if
#'   change in penalized negative log-likelihood drops below EM.epsilon.
#' @param EM.verbose (Optional) Integer in {0,1,2}. verbose = 0: No output.
#'   verbose = 1: Print summary at the end of the optimization. verbose = 2:
#'   Print progress during optimization.
#' @param sample.data (Optional) Sample data according to posterior probability
#'   or not.
#' @param TGGL.mu (Optional) Mu parameter for TGGL.
#' @param TGGL.epsilon (Optional) Epsilon parameter for TGGL.
#' @param shrink.mu (Optional) Multiply mu by min(lambda, 1).
#'
#' @return List containing
#' \item{results}{List of TGGLMix models for each parameter setting.}
#' \item{top.model}{Model with highest predictive likelihood.}
#'
#' @seealso \code{\link{TGGLMix}}
#' @export
RunTGGLMixSelection <- function(X = NULL, task.specific.features = list(), Y, M.vec,
                       validation.ids = NULL, num.starts = 1, num.threads = NULL,
                       groups, weights, lambda.vec, verbose = 0,
                       gam = 1, homoscedastic = FALSE,
                       EM.max.iter = 200, EM.epsilon = 1e-5,
                       EM.verbose = 0, sample.data = FALSE,
                       TGGL.mu = 1e-5, TGGL.epsilon = 1e-5, shrink.mu = TRUE) {
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

  # separate into training and validation set
  if (!is.null(validation.ids)) {
    train.ids <- setdiff(1:N, validation.ids)
    if (J1 > 0) {
      X.test <- X[validation.ids, ]
      X <- X[train.ids, ]
    } else {
      X.test <- NULL
    }
    if (J2 > 0) {
      tsf.test <- lapply(task.specific.features, FUN = function(x){x[validation.ids, ]})
      task.specific.features <- lapply(task.specific.features, FUN = function(x){x[train.ids, ]})
    } else {
      tsf.test <- list()
    }
    Y.test <- Y[validation.ids, ]
    Y <- Y[train.ids, ]
  } else {
    if (nrow(parameter.grid) > 1) {
      stop("Please provide validation.ids when passing multiple model parameters.")
    }
  }

  RunParameter <- function(ind) {
    # run model for entry of parameter grid
    M <- parameter.grid[ind, 1]
    lambda <- parameter.grid[ind, 2]

    if (shrink.mu) {
      mu <- TGGL.mu * min(lambda, 1)
    } else {
      mu <- TGGL.mu
    }

    train.start.time <- Sys.time()
    # run TGGLMix multiple times and keep best local optimum.
    top.model <- NULL
    for (st in 1:num.starts) {
      tggl.mix <- TGGLMix(X = X, task.specific.features = task.specific.features,
                          Y = Y, M = M, groups = groups, weights = weights,
                          lambda = lambda,
                          homoscedastic = homoscedastic, gam = gam,
                          EM.max.iter = EM.max.iter, EM.epsilon = EM.epsilon,
                          EM.verbose = EM.verbose, sample.data = sample.data,
                          TGGL.mu = mu, TGGL.epsilon = TGGL.epsilon)
      if (verbose > 2) {
        print(sprintf('[M = %d, lambda = %.e] Start %d - PenNegLL: %.3f, LL: %.3f.',
                      M, lambda, st, tggl.mix$obj, tggl.mix$loglik))
      }
      if (is.null(top.model)) {
        top.model <- tggl.mix
      } else if (top.model$obj > tggl.mix$obj) {
        top.model <- tggl.mix
      }
    }
    train.end.time <- Sys.time()
    train.time <- as.numeric(train.end.time - train.start.time, units = "mins")

    # determine validation likelihood (if applicable)
    validation.loglik <- NA
    if (!is.null(validation.ids)) {
      test.stats <- ComputeLogLikelihood(top.model$models,
                                         top.model$prior,
                                         top.model$sigmas,
                                         X = X.test,
                                         task.specific.features = tsf.test,
                                         Y = Y.test)
      validation.loglik <- test.stats$loglik
    }
    if (verbose > 1) {
      if (is.null(validation.ids)) {
        print(sprintf('Trained TGGLMix [M = %d, lambda = %.e], %0.1f min. PenNegLL: %.3f, LL: %.3f.',
                      M, lambda, train.time, top.model$obj, top.model$loglik))
      } else {
        print(sprintf('Trained TGGLMix [M = %d, lambda = %.e], %0.1f min. PenNegLL: %.3f, Train-LL: %.3f, Validation-LL: %.3f.',
                      M, lambda, train.time, top.model$obj, top.model$loglik, test.stats$loglik))
      }
    }
    return(list(models = top.model$models,
                posterior = top.model$posterior,
                prior = top.model$prior,
                sigmas = top.model$sigmas,
                obj = top.model$obj,
                train.loglik = top.model$loglik,
                validation.loglik = validation.loglik,
                groups = groups,
                weights = weights,
                lambda = lambda))
  }

  if (verbose > 0) {
    print(sprintf("Running TGGL-Mixture for %d parameter settings ... ", nrow(parameter.grid)))
  }
  sel.start.time <- Sys.time()
  doMC::registerDoMC(num.threads)
  results <- foreach(l = 1:nrow(parameter.grid)) %dopar% RunParameter(l)
  sel.end.time <- Sys.time()
  sel.time <- as.numeric(sel.end.time - sel.start.time, units = "mins")

  if (is.null(validation.ids)) {
    opt.model.idx <- 1
  } else {
    opt.model.idx <- which.max(sapply(results, FUN = function(x){x$validation.loglik}))
  }

  if (verbose > 0) {
    print(sprintf("Minutes to run Model Selection %0.1f. Best Parameter Setting: M = %d, lambda = %.e.",
                  sel.time, parameter.grid[opt.model.idx, 1], parameter.grid[opt.model.idx, 2]))
  }

  return(list(results = results, top.model = results[[opt.model.idx]]))
}
