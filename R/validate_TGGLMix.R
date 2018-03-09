#' Convenience function for running a tree-guided group lasso mixture model
#' model (TGGLMix) multiple times.
#'
#' Fit a tree-guided group lasso mixture model. Restart with different random
#' initializations and keep the model with the lowest objective value. Optional:
#' Evaluate best model on test set. By default all data is used for training. If
#' test.ids is not NULL, exclude corresponding indices from training and use
#' them for testing instead.
#'
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param M Number of Clusters.
#' @param test.ids (Optional) Indices of data points to be used for testing.
#'   Default is to use all data for training.
#' @param num.starts (Optional) Number of starts. Default is 1 (no restarts).
#' @param num.threads (Optional) Number of threads to be used. Default is 1.
#' @param groups Binary V by K matrix determining group membership: Task k in
#'   group v iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param lambda Regularization parameter.
#' @param verbose (Optional) Integer in {0,1,2}. verbose = 0: No output.
#'   verbose = 1: Print summary at the end of the optimization. verbose = 2:
#'   Print progress during optimization.
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
#'
#'
#' @return List containing
#' \item{models}{List of TGGL models for each component.}
#' \item{posterior}{N by M Matrix containing posterior probabilities.}
#' \item{prior}{Vector with prior probabilities for each component.}
#' \item{sigmas}{M by K Matrix with standard deviations for each component.}
#' \item{obj}{Penalized negative likelihood (final objective value).}
#' \item{train.loglik}{Likelihood for training data.}
#' \item{test.loglik}{Likelihood for test data.}
#'
#' @seealso \code{\link{TGGLMix}}
#' @export
RunTGGLMix <- function(X = NULL, task.specific.features = list(), Y, M,
                       test.ids = NULL, num.starts = 1, num.threads = NULL,
                       groups, weights, lambda, verbose = 0,
                       gam = 1, homoscedastic = FALSE,
                       EM.max.iter = 200, EM.epsilon = 1e-5,
                       EM.verbose = 0, sample.data = FALSE,
                       TGGL.mu = 1e-5, TGGL.epsilon = 1e-5) {
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

  # separate into training and test set
  if (!is.null(test.ids)) {
    train.ids <- setdiff(1:N, test.ids)
    if (J1 > 0) {
      X.test <- X[test.ids, ]
      X <- X[train.ids, ]
    } else {
      X.test <- NULL
    }
    if (J2 > 0) {
      tsf.test <- lapply(task.specific.features, FUN = function(x){x[test.ids, ]})
      task.specific.features <- lapply(task.specific.features, FUN = function(x){x[train.ids, ]})
    } else {
      tsf.test <- list()
    }
    Y.test <- Y[test.ids, ]
    Y <- Y[train.ids, ]
  }

  RunModel <- function(start.id) {
    # Run TGGLMixture
    tggl.mix <- TGGLMix(X = X, task.specific.features = task.specific.features,
                        Y = Y, M = M, groups = groups, weights = weights,
                        lambda = lambda,
                        homoscedastic = homoscedastic, gam = gam,
                        EM.max.iter = EM.max.iter, EM.epsilon = EM.epsilon,
                        EM.verbose = EM.verbose, sample.data = sample.data,
                        TGGL.mu = TGGL.mu, TGGL.epsilon = TGGL.epsilon)
    if (verbose > 1) {
      print(sprintf('Start %d - PenNegLL: %.3f, LL: %.3f.',
                    start.id, tggl.mix$obj, tggl.mix$loglik))
    }
    return(tggl.mix)
  }

  if (verbose > 1) {
    print(sprintf("Running TGGLMix [M = %d, lambda = %.e]. Number of starts: %d ... ",
                  M, lambda, num.starts))
  }
  doMC::registerDoMC(num.threads)
  train.start.time <- Sys.time()
  results <- foreach(l = 1:num.starts) %dopar% RunModel(l)
  train.end.time <- Sys.time()
  top.model.idx <- which.min(sapply(results, FUN = function(x){x$obj}))
  top.model <- results[[top.model.idx]]
  train.time <- as.numeric(train.end.time - train.start.time, units = "mins")
  test.loglik <- NULL
  if (!is.null(test.ids)) {
    test.stats <- ComputeLogLikelihood(top.model$models, top.model$prior, top.model$sigmas,
                                        X = X.test, task.specific.features = tsf.test, Y = Y.test)
  }
  if (verbose > 0) {
    if (is.null(test.ids)) {
      print(sprintf('Trained TGGLMix [M = %d, lambda = %.e], %0.1f min. PenNegLL: %.3f, LL: %.3f.',
                    M, lambda, train.time, top.model$obj, top.model$loglik))
    } else {
      print(sprintf('Trained TGGLMix [M = %d, lambda = %.e], %0.1f min. PenNegLL: %.3f, Train-LL: %.3f, Test-LL: %.3f.',
                    M, lambda, train.time, top.model$obj, top.model$loglik, test.stats$loglik))
    }
  }
  return(list(models = top.model$models,
              posterior = top.model$posterior,
              prior = top.model$prior,
              sigmas = top.model$sigmas,
              obj = top.model$obj,
              train.loglik = top.model$loglik,
              test.loglik = test.stats$loglik))
}
