################################################################################
# Utility functions for linear multi-task models.
################################################################################

#' Compute predictions for linear multi-task model.
#'
#' Compute linear predictions from predictor variables and regression
#' coefficients.
#'
#' @param LMTL.model Linear multi-task learning model (list containing B and
#'   intercept).
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#'
#' @return N by J predictions matrix.
#'
#' @export
MTPredict <- function(LMTL.model, X = NULL, task.specific.features = list()) {

  if (length(intersect(names(LMTL.model), c("B", "intercept"))) < 2) {
    stop("No valid model supplied.")
  }
  K <- ncol(LMTL.model$B)

  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # check for shared features
  J1 <- 0
  if (!is.null(X)) {
    J1 <- ncol(X)
  }

  # check for task specific features
  J2 <- 0
  if (length(task.specific.features) > 0) {
    if (length(task.specific.features) != K) {
      stop("B must have one column per element in task.specific.features!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  J <- J1 + J2
  if (nrow(LMTL.model$B) != J) {
    stop("Incorrect number of dimensions: B must have one row per feature!")
  }

  # compute predictions for shared feature matrix
  if (J1 > 0) {
    prediction <- X %*% LMTL.model$B[1:J1, ]
  } else {
    # we must have task specific features
    C <- task.specific.features[[1]]
    prediction <- matrix(0, nrow = nrow(C), ncol = K)
  }

  # compute predictions for task specific feature matrices
  if (J2 > 0) {
    # rows J1 + 1 and below of B correspond to task specific features
    tsf.idx <- (J1 + 1):J
    for (k in 1:K) {
      C <- task.specific.features[[k]]
      prediction[, k] <- prediction[, k] + C %*% LMTL.model$B[tsf.idx, k]
    }
  }

  # add intercepts
  prediction <- sweep(prediction, 2, FUN = "+", LMTL.model$intercept)
  return(prediction)
}

#' Compute (mean) squared error for linear multi-task model.
#'
#' @param LMTL.model Linear multi-task learning model (list containing B and
#'   intercept).
#' @param Y N by K output matrix for every task.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 matrix for one particular task
#'   (where columns are features). List has to be ordered according to the
#'   columns of Y.
#' @param pred Predicted output matrix. If NULL, compute predictions using input
#'   features.
#' @param normalize Compute mean (TRUE) or sum (FALSE).
#' @param aggregate.tasks Aggregate results over all tasks (TRUE) or return task
#'   specific errors (FALSE).
#'
#' @return The (mean) squared error between predictions for each task and Y.
#'
#' @export
MTComputeError <- function (LMTL.model, Y, X = NULL,
                            task.specific.features = list(),
                            pred = NULL, normalize = TRUE,
                            aggregate.tasks = TRUE) {

  if (is.null(pred)) {
    pred <- MTPredict(LMTL.model = LMTL.model, X = X,
                      task.specific.features = task.specific.features)
  } else {
    if (!isTRUE(all.equal(dim(pred), dim(Y)))) {
      stop("Dimensions of pred and Y have to coincide!")
    }
  }

  D <- (Y - pred)^2
  if (normalize) {
    # compute mean
    if (aggregate.tasks) {
      return(mean(D))
    } else {
      return(colMeans(D))
    }
  } else {
    # compute sum
    if (aggregate.tasks) {
      return(sum(D))
    } else {
      return(colSums(D))
    }
  }
}

#' Compute (mean) correlation for linear multi-task model.
#'
#' @param LMTL.model Linear multi-task learning model (list containing B and
#'   intercept).
#' @param Y N by K output matrix for every task.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param pred Predicted output matrix. If NULL, compute predictions using input
#'   features.
#' @param method Correlation method to use.
#'
#' @return The mean correlation between predictions for each task and Y.
#'
#' @importFrom stats cor
#' @export
MTComputeMeanCorrelation <- function (LMTL.model, Y, X = NULL,
                                      task.specific.features = list(),
                                      pred = NULL, method = "spearman") {

  if (is.null(pred)) {
    pred <- MTPredict(LMTL.model = LMTL.model, X = X,
                      task.specific.features = task.specific.features)
  } else {
    if (!isTRUE(all.equal(dim(pred), dim(Y)))) {
      stop("Dimensions of pred and Y have to coincide!")
    }
  }
  return(mean(diag(cor(pred, Y, method = method)), na.rm = TRUE))
}




#' Precompute matrices XTX and XTY for \code{\link{TreeGuidedGroupLasso}}.
#'
#' Computes XTX and XTY. If task specific features are supplied, the result will
#' be a list containing the respective matrices for each task.
#'
#' @param Y N by K output matrix for every task.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 matrix for one particular task
#'   (where columns are features). List has to be ordered according to the
#'   columns of Y.
#' @param standardize Standardize data (default is TRUE).
#' @param row.weights Observation weights.
#' @param return.inputs Return transformed (scaled and weighted, if applicable)
#'   version of X, Y and task.specific.features.
#'
#' @return List containing Y, X, task.specific.features, XTX, XTY, XT1 and YT1.
#'   Also return means and standard deviations for Y and all inputs.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}}, \code{\link{RunGroupCrossvalidation}}.
#' @importFrom stats sd
#' @export
PrepareMatrices <- function(Y, X = NULL,
                            task.specific.features = list(),
                            standardize = TRUE,
                            row.weights = NULL,
                            return.inputs = FALSE) {

  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # initialize variables
  K <- ncol(Y)
  N <- nrow(Y)
  J1 <- 0
  if (!is.null(X)) {
    J1 <- ncol(X)
  }
  J2 <- 0
  if(length(task.specific.features) > 0) {
    J2 <- ncol(task.specific.features[[1]])
  }
  J <- J1 + J2

  if (!is.null(row.weights)) {
    #####################
    # apply row.weights #
    #####################
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

  ##########################
  # compute sample moments #
  ##########################
  X.sds <- NULL
  X.mus <- NULL
  tsf.sds <- list()
  tsf.mus <- list()
  Y.sds <- NULL
  Y.mus <- NULL

  if (J1 > 0) {
    X.sds <- apply(X, 2, sd)
    X.mus <- apply(X, 2, mean)
  }
  if (J2 > 0) {
    tsf.sds <- lapply(task.specific.features, FUN = function(A){apply(A, 2, sd)})
    tsf.mus <- lapply(task.specific.features, FUN = function(A){apply(A, 2, mean)})
  }
  Y.mus <- apply(Y, 2, mean)
  Y.sds <- apply(Y, 2, sd)

  sample.moments <- list(X.mus = X.mus, X.sds = X.sds,
                         tsf.mus = tsf.mus, tsf.sds = tsf.sds,
                         Y.mus = Y.mus, Y.sds = Y.sds)

  ###########################
  # standardize if required #
  ###########################
  if (standardize) {
    if (J1 > 0) {
      X <- scale(X)
    }
    if (J2 > 0) {
      task.specific.features <- lapply(task.specific.features, scale)
    }
    Y <- scale(Y)
  }

  ##############################
  # precompute matrix products #
  ##############################
  YT1 <- t(Y) %*% row.weights.sqrt
  if (J2 == 0) {
    # no task specific features
    XTX <- t(X) %*% X
    XT1 <- t(X) %*% row.weights.sqrt
    XTY <- t(X) %*% Y
  } else {
    # task specific features
    XTX <- list()
    XT1 <- list()
    XTY <- matrix(0, J, K)
    for (k in 1:K) {
      mat <- cbind(X, task.specific.features[[k]])
      XTX[[k]] <- t(mat) %*% mat
      XT1[[k]] <- t(mat) %*% row.weights.sqrt
      XTY[, k] <- t(mat) %*% Y[, k]
    }
  }
  if (return.inputs) {
    return(list(XTX = XTX, XTY = XTY, XT1 = XT1, YT1 = YT1,
                task.specific.features = task.specific.features,
                X = X, Y = Y,
                sample.moments = sample.moments))
  } else {
    return(list(XTX = XTX, XTY = XTY, XT1 = XT1, YT1 = YT1,
                sample.moments = sample.moments))
  }
}
