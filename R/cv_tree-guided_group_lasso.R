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
#' @param verbose (Optional) Integer in {-2, -1, 0,1,2}. verbose = -2: No
#'   output. verbose = -1: Display total elapsed time. verbose = 0: Display
#'   elapsed time for every parameter. verbose = 1: Print summary at the end of
#'   the optimization. verbose = 2: Print progress during optimization.
#' @param row.weights (Optional) Use weighted MSE.
#' @param standardize (Optional) Default is TRUE. Standardize data (using R
#'   function scale()). Coefficients will be returned on original scale.
#' @param fit.intercept (Optional) Default is TRUE. Include intercept.
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
                                     num.threads = 1, verbose = -1, row.weights = NULL,
                                     standardize = TRUE, fit.intercept = TRUE, ...) {
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

  # set input weights
  if (is.null(row.weights)) {
    row.weights <- rep(1, N)
  }
  # divide data into folds
  cv.folds <- split(sample(N), 1:num.folds)

  # build parameter grid
  parameter.grid <- cbind(lambda.vec, matrix(rep(weights.matrix, each = length(lambda.vec)),
                                             ncol = ncol(weights.matrix)))

  if (verbose > -2) {
    print(sprintf("Running group crossvalidation for %d parameter settings ... ", nrow(parameter.grid)))
  }
  cv.start.time <- Sys.time()

  if (!is.null(row.weights)) {
    # apply row weights
    row.weights.sqrt <- sqrt(row.weights)
    if (J1 > 0) {
      X <- X * row.weights.sqrt
    }
    if (J2 > 0) {
      task.specific.features <- lapply(task.specific.features,
                                       FUN = function(A) {A * row.weights.sqrt})
    }
    Y <- Y * row.weights.sqrt
  } else {
    row.weights.sqrt <- rep(1, N)
  }

  if (standardize) {
    return.inputs <- TRUE
  } else {
    return.inputs <- FALSE
  }

  PrepareMatricesForFold <- function(idx) {
    ids <- cv.folds[[idx]]
    # exclude indices for this fold
    fold.ex.task.specific.features <- lapply(task.specific.features, function (x){x[-ids, ]})
    # prepare matrices
    pm.res <- PrepareMatrices(Y = Y[-ids, , drop = FALSE], X = X[-ids, ],
                              task.specific.features = fold.ex.task.specific.features,
                              standardize = standardize,
                              row.weights = NULL,
                              return.inputs = return.inputs)
    return(pm.res)
  }

  # precompute matrices to save computation time
  if (verbose > -1) {
    print("Computing immutable matrices ... ")
  }
  precomp.start.time <- Sys.time()
  doMC::registerDoMC(num.threads)
  cached.mats <- foreach(i = seq_along(cv.folds)) %dopar% PrepareMatricesForFold(i)

  MSE.Lipschitz.list <- list()
  for (i in seq_along(cached.mats)) {
    if (J2 > 0) {
      # too expensive:
      #L1 <- max(unlist(lapply(XTX, FUN = function(M) {max(eigen(M)$values)})))
      # instead just consider first matrix
      MSE.Lipschitz.list[[i]] <- max(eigen(cached.mats[[i]]$XTX[[1]])$values)
    } else {
      MSE.Lipschitz.list[[i]] <- max(eigen(cached.mats[[i]]$XTX)$values)
    }
  }

  precomp.end.time <- Sys.time()
  precomp.time <- as.numeric(precomp.end.time - precomp.start.time, units = "mins")

  RunParameter <- function (ind) {
    # Run crossvalidation on given row of parameter grid.
    #
    # Args:
    #   ind: row of parameter grid.
    #
    # Returns:
    #   List containing hold-out error for given parameter set and
    #   whether or not the optimization procedure ended prematurely.

    if (verbose > -1) {
      print(sprintf('Crossvalidation for parameter set: %d', ind))
    }

    lambda <- parameter.grid[ind, 1]
    weights <- parameter.grid[ind, 2:ncol(parameter.grid)]

    error <- rep(0, K)
    early.termination <- TRUE
    param.start.time <- Sys.time()
    for (i in 1:length (cv.folds)) {
      fold <- cv.folds[[i]]
      # restrict to subset of the data
      if (standardize) {
        # use scaled matrices
        fold.ex.X <- cached.mats[[i]]$X
        fold.ex.task.specific.features <- cached.mats[[i]]$task.specific.features
        fold.ex.Y <- cached.mats[[i]]$Y
      } else {
        fold.ex.X <- X[-fold, ]
        fold.ex.task.specific.features <- lapply(task.specific.features, function(x){x[-fold, ]})
        fold.ex.Y <- Y[-fold, , drop = FALSE]
      }

      # train model
      fold.result <- TreeGuidedGroupLasso(X = fold.ex.X,
                                          task.specific.features = fold.ex.task.specific.features,
                                          Y = fold.ex.Y,
                                          groups = groups, weights = weights, lambda = lambda,
                                          cached.mats = cached.mats[[i]],
                                          MSE.Lipschitz = MSE.Lipschitz.list[[i]],
                                          verbose = max(verbose, 0),
                                          standardize = standardize,
                                          fit.intercept = fit.intercept,
                                          row.weights = row.weights, ...)
      early.termination <- early.termination & fold.result$early.termination

      fold.task.specific.features <- lapply(task.specific.features, function(x){x[fold, ]})
      error <- error + MTComputeError(LMTL.model = fold.result,
                                      Y = Y[fold, , drop = FALSE], X = X[fold, ],
                                      task.specific.features = fold.task.specific.features,
                                      aggregate.tasks = FALSE)
    }

    # compute hold-out error
    error <- error / num.folds

    param.end.time <- Sys.time()
    if (verbose > -1) {
      print(sprintf('Minutes to run parameter set %d: %0.1f',
                    ind, as.numeric(param.end.time-param.start.time, units = "mins")))
    }

    return(list(lambda = lambda, weights = weights,
                error = error, early.termination = early.termination))
  }

  # run cv on all parameter settings
  doMC::registerDoMC(num.threads)
  param.inds <- 1:nrow(parameter.grid)
  cv.results <- foreach(l = param.inds) %dopar% RunParameter(l)

  # turn result into data.frame
  weight.names <- 1:length(cv.results[[1]]$weights)
  cv.results <- lapply(cv.results, unlist)
  cv.results <- data.frame(t(sapply(cv.results, c)))
  colnames(cv.results) <- c("lambda", weight.names, paste("Task", 1:K, sep = ""), "early.termination")

  cv.end.time <- Sys.time()
  cv.time <- as.numeric(cv.end.time - cv.start.time, units = "mins")

  if (verbose > -1) {
    print("Training model on full data set ... ")
  }
  train.start.time <- Sys.time()

  # identify best parameters
  min.idx <- which.min(rowMeans(cv.results[, paste("Task", 1:K, sep = "")]))
  lambda <- parameter.grid[min.idx, 1]
  weights <- parameter.grid[min.idx, 2:ncol(parameter.grid)]

  # retrain model
  full.model <- TreeGuidedGroupLasso(X = X, task.specific.features = task.specific.features,
                                     Y = Y, groups = groups, weights = weights, lambda = lambda,
                                     verbose = max(verbose, 0),
                                     standardize = standardize,
                                     fit.intercept = fit.intercept,
                                     row.weights = row.weights, ...)
  train.end.time <- Sys.time()
  train.time <- as.numeric(train.end.time - train.start.time, units = "mins")
  if (verbose > -2) {
    print(sprintf("Minutes to run CV (precomp. + CV + final model): %0.1f + %0.1f + %0.1f = %0.1f.",
                  precomp.time, cv.time, train.time, precomp.time + cv.time + train.time))
  }

  return(list(cv.results = cv.results, full.model = full.model))
}
