#' Cross-validation wrapper for task-by-task lasso model.
#'
#' Perform k-fold cross-validation for each task using \code{\link[glmnet]{cv.glmnet}}.
#'
#' @param X Column centered N by J input matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param Y Column centered N by K output matrix for every task.
#' @param lambda.vec Vector of regularization parameters.
#' @param num.folds Number of folds.
#' @param num.threads Number of threads to use.
#' @param ... Additional parameters passed to
#'   \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return List containing
#'    \item{lambda}{Best lambda for each task.}
#'    \item{beta}{Final regression coefficients for the model fitted on the full data set.}
#'    \item{error}{Cross-validation errors for each task.}
#'
#' @export
RunTBTCrossvalidation <- function (X = NULL, task.specific.features = list(), Y,
                                   lambda.vec, num.folds = 10, num.threads = 1, ...) {

  # initialization and error checking
  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # check for shared features
  J1 <- 0
  if (!is.null(X)) {
    if (nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows!")
    }
    J1 <- ncol(X)
  }

  # check for task specific features
  J2 <- 0
  if (length(task.specific.features > 0)) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  K <- ncol(Y)
  N <- nrow(Y)
  J <- J1 + J2

  RunParameter <- function (ind) {
    # Run crossvalidation for given task.
    #
    # Args:
    #   ind: Task index.
    #
    # Returns:
    #   List with best lambda, regression coefficients for the model
    #   trained on the full data set and the cv error.

    print(sprintf('Crossvalidation for task: %d', ind))

    error <- 0
    early.termination <- TRUE
    task.start.time <- Sys.time()

    if (J2 == 0) {
      # no task specific features
      mat <- X
    } else {
      # with task specific features
      if (is.null(X)) {
        mat <- task.specific.features[[ind]]
      } else {
        mat <- cbind(X, task.specific.features[[ind]])
      }
    }
    cv.results <- glmnet::cv.glmnet(x = mat, y = Y[, ind],
                                    lambda = lambda.vec, intercept = FALSE, ...)

    task.end.time <- Sys.time()
    print(sprintf('Minutes to run cv for task %d: %0.1f',
                  ind, as.numeric(task.end.time-task.start.time, units = "mins")))


    min.err <- which.min(cv.results$cvm)
    error <- min(cv.results$cvm)
    beta <- cv.results$glmnet.fit$beta[, min.err]
    lambda <- cv.results$lambda[min.err]
    return(list(lambda = lambda, beta = beta, error = error))
  }

  tbt.beta <- matrix(0, nrow = J, ncol = K)
  # store best lambda for each task
  lambda <- rep(0, ncol(X))
  # store MSE for each task
  err <- rep(0, ncol(X))

  doMC::registerDoMC(num.threads)
  cv.results <- foreach::foreach(task = 1:K) %dopar% RunParameter(task)
  for (task in 1:K) {
    tbt.beta[, task] <- cv.results[[task]]$beta
    lambda[task] <- cv.results[[task]]$lambda
    err[task] <- cv.results[[task]]$error
  }

  return(list(lambda = lambda, beta = tbt.beta, error = err))
}
