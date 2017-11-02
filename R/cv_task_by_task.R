#' Cross-validation wrapper for task-by-task lasso model.
#'
#' Perform k-fold cross-validation for each task using \code{\link[glmnet]{cv.glmnet}}.
#'
#' @param X N by J input matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to
#'   each task. Each entry contains an N by J2 matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param Y N by K output matrix for every task.
#' @param lambda.vec Vector of regularization parameters.
#' @param num.folds Number of folds.
#' @param num.threads Number of threads to use.
#' @param ... Additional parameters passed to
#'   \code{\link[glmnet]{cv.glmnet}}.
#'
#' @return List containing
#'     \item{cv.results}{data.frame with cross-validation errors for
#'     different parameters.}
#'     \item{full.model}{Full model trained on the whole data set.}
#'
#' @importFrom foreach foreach %dopar%
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
  if (length(task.specific.features) > 0) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  K <- ncol(Y)
  N <- nrow(Y)
  J <- J1 + J2

  RunTask <- function (ind) {
    # Run crossvalidation for given task.
    #
    # Args:
    #   ind: Task index.
    #
    # Returns:
    #   List with lambdas and cv errors and regression coefficients for the
    #   best model trained on the full data set.

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

    cv.results <- glmnet::cv.glmnet(x = mat, y = Y[, ind], lambda = lambda.vec, ...)

    task.end.time <- Sys.time()
    print(sprintf('Minutes to run cv for task %d: %0.1f',
                  ind, as.numeric(task.end.time-task.start.time, units = "mins")))

    cvm <- cv.results$cvm[order(cv.results$lambda)]
    top.coef <- coef(cv.results, s = "lambda.min")
    return(list(cvm = cv.results$cvm, lambda.min = cv.results$lambda.min,
                B = top.coef[-1], intercept = top.coef[1]))
  }

  # run all tasks
  doMC::registerDoMC(num.threads)
  all.cv.results <- foreach(task = 1:K) %dopar% RunTask(task)

  # extract full model
  B <- matrix(0, nrow = J, ncol = K)
  intercept <- rep(0, K)
  lambda <- rep(0, K)
  cv.results <- matrix(0, nrow = length(lambda.vec), ncol = K + 1)
  cv.results[, 1] <- sort(lambda.vec)
  for (task in 1:K) {
    B[, task] <- all.cv.results[[task]]$B
    intercept[task] <- all.cv.results[[task]]$intercept
    lambda[task] <- all.cv.results[[task]]$lambda.min
    cv.results[, task + 1] <- all.cv.results[[task]]$cvm
  }
  colnames(cv.results) <- c("lambda", paste("Task", 1:K, sep = ""))
  full.model <- list(lambda = lambda, B = B, intercept = intercept)

  return(list(cv.results = cv.results, full.model = full.model))
}
