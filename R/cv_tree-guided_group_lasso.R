#' Cross-validation wrapper for \code{\link{TreeGuidedGroupLasso}}.
#'
#' Perform k-fold cross-validation to find optimal parameters settings.
#'
#' @param X Column centered N by J input matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param Y Column centered N by K output matrix for every task.
#' @param groups V by K matrix determining group membership: Task k in group v
#'   iff groups[v,k] == 1 weights.matrix: Numerical matrix where each row is a
#'   numerical vector mapping groups to their weight.
#' @param weights.matrix Numerical matrix where each row is a numerical vector
#'   mapping groups to their weight.
#' @param lambda.vec Vector of regularization parameters.
#' @param num.folds Number of folds.
#' @param num.threads Number of threads to use.
#' @param ... Additional parameters passed to
#'   \code{\link{TreeGuidedGroupLasso}}.
#'
#' @return List containing
#'     \item{cv.results}{data.frame with cross-validation errors for
#'     different parameters.}
#'     \item{full.model}{Full model trained on the whole data set.}
#'
#' @importFrom foreach foreach %dopar%
#' @seealso \code{\link{TreeGuidedGroupLasso}}
#' @export
RunGroupCrossvalidation <- function (X = NULL, task.specific.features = list(), Y,
                                     groups, weights.matrix, lambda.vec, num.folds = 10,
                                     num.threads = 1, ...) {
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

  if (ncol(Y) != ncol(groups)) {
    stop("Y and groups must have the same number of columns!")
  }
  if (ncol(weights.matrix) != nrow(groups)) {
    stop("Length of weights has to equal number of rows of groups!")
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  # divide data into folds
  cv.folds <- split(sample(N), 1:num.folds)

  # precompute matrices to save computation time
  print("Computing immutable matrices ... ")
  precomp.start.time <- Sys.time()

  XTX.global <- list()
  XTY.global <- list()
  doMC::registerDoMC(num.threads)

  PrepareMatricesForFold <- function(idx) {
    ids <- cv.folds[[idx]]
    # exclude indices for this fold
    fold.ex.task.specific.features <- lapply(task.specific.features, function (x) {x[-ids, ]} )
    # prepare matrices
    pm.res <- PrepareMatrices(Y = Y[-ids, ], X = X[-ids, ],
                              task.specific.features = fold.ex.task.specific.features,
                              standardize = TRUE)
    return(pm.res)
  }

  mats <- foreach(i = seq_along(cv.folds)) %dopar% PrepareMatricesForFold(i)
  for (i in seq_along(mats)) {
    XTX.global[[i]] <- mats[[i]]$XTX
    XTY.global[[i]] <- mats[[i]]$XTY
  }

  # build parameter grid
  parameter.grid <- cbind(lambda.vec, matrix(rep(weights.matrix, each = length(lambda.vec)),
                                             ncol = ncol(weights.matrix)))

  RunParameter <- function (ind) {
    # Run crossvalidation on given row of parameter grid.
    #
    # Args:
    #   ind: row of parameter grid.
    #
    # Returns:
    #   List containing hold-out error for given parameter set and
    #   whether or not the optimization procedure ended prematurely.

    print(sprintf('Crossvalidation for parameter set: %d', ind))

    lambda <- parameter.grid[ind, 1]
    weights <- parameter.grid[ind, 2:ncol(parameter.grid)]

    error <- 0
    early.termination <- TRUE
    param.start.time <- Sys.time()
    for (i in 1:length (cv.folds)) {
      fold <- cv.folds[[i]]
      # restrict to subset of the data
      if (J2 > 0) {
        fold.ex.task.specific.features <- lapply(task.specific.features, function (x) {x[-fold, ]} )
        fold.task.specific.features <- lapply(task.specific.features, function (x) {x[fold, ]} )
      } else {
        fold.ex.task.specific.features <- list()
        fold.task.specific.features <- list()
      }

      # train model
      fold.result <- TreeGuidedGroupLasso(X = X[-fold, ], task.specific.features = fold.ex.task.specific.features,
                                          Y = Y[-fold,,drop = FALSE], groups = groups, weights = weights, lambda = lambda,
                                          XTX = XTX.global[[i]], XTY = XTY.global[[i]], ...)
      early.termination <- early.termination & fold.result$early.termination
      error <- error + MTComputeError(LMTL.model = fold.result, Y = Y[fold, , drop = FALSE],
                                      X = X[fold, ], task.specific.features = fold.task.specific.features)
    }

    # compute hold-out error
    error <- error / num.folds

    param.end.time <- Sys.time()
    print(sprintf('Minutes to run parameter set %d: %0.1f', ind, as.numeric(param.end.time-param.start.time, units = "mins")))

    return(list(lambda = lambda, weights = weights, error = error, early.termination = early.termination))
  }

  print(sprintf("Running group crossvalidation on %d parameter settings ... ", nrow(parameter.grid)))
  cv.start.time <- Sys.time()

  # run cv on all parameter settings
  doMC::registerDoMC(num.threads)
  param.inds <- 1:nrow(parameter.grid)
  cv.results <- foreach(l = param.inds) %dopar% RunParameter(l)

  # turn result into data.frame
  weight.names <- 1:length(cv.results[[1]]$weights)
  cv.results <- lapply(cv.results, unlist)
  cv.results <- data.frame(t(sapply(cv.results, c)))
  colnames(cv.results) <- c("lambda", weight.names, "cv.error", "early.termination")

  cv.end.time <- Sys.time()
  print(sprintf("Minutes to run crossvalidation : %0.1f", as.numeric(cv.end.time - cv.start.time, units = "mins")))

  print("Training model on full data set ... ")
  train.start.time <- Sys.time()

  # identify best parameters
  min.idx <- which.min(cv.results$cv.error)
  lambda <- parameter.grid[min.idx, 1]
  weights <- parameter.grid[min.idx, 2:ncol(parameter.grid)]

  # retrain model
  full.model <- TreeGuidedGroupLasso(X = X, task.specific.features = task.specific.features,
                                     Y = Y, groups = groups, weights = weights, lambda = lambda, ...)
  train.end.time <- Sys.time()
  print(sprintf("Minutes to run train full model : %0.1f", as.numeric(train.end.time - train.start.time, units = "mins")))

  return(list(cv.results = cv.results, full.model = full.model))
}
