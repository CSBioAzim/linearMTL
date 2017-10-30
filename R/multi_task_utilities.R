################################################################################
# Utility functions for linear multi-task models.
################################################################################

#' Compute predictions for linear multi-task model.
#'
#' Compute linear predictions from predictor variables and regression
#' coefficients.
#'
#' @param B J by K matrix of regression coefficients, J = J1 + J2.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#'
#' @return N by J predictions matrix.
#'
#' @export
MTPredict <- function(B, X = NULL, task.specific.features = list()) {

  K <- ncol(B)

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
  if (nrow(B) != J) {
    stop("Incorrect number of dimensions: B must have one row per feature!")
  }

  # compute predictions for shared feature matrix
  if (J1 > 0) {
    prediction <- X %*% B[1:J1, ]
  } else {
    # we must have task specific features
    C <- task.specific.features[[1]]
    prediction <- matrix(0, nrow = nrow(C), ncol = ncol(B))
  }

  # compute predictions for task specific feature matrices
  if (J2 > 0) {
    # rows J1 + 1 and below of B correspond to task specific features
    tsf.idx <- (J1 + 1):nrow(B)
    for (k in 1:K) {
      C <- task.specific.features[[k]]
      prediction[, k] <- prediction[, k] + C %*% B[tsf.idx, k]
    }
  }
  return(prediction)
}

#' Compute (mean) squared error for linear multi-task model.
#'
#' @param Y Column centered N by K output matrix for every task.
#' @param B J by K matrix of regression coefficients, J = J1 + J2.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param pred Predicted output matrix. If NULL, compute predictions using input features.
#' @param normalize Compute mean (TRUE) or sum (FALSE).
#'
#' @return The (mean) squared error between predictions for each task and Y.
#'
#' @export
MTComputeError <- function (Y, B, X = NULL, task.specific.features = list(),
                            pred = NULL, normalize = TRUE) {

  if (is.null(pred)) {
    pred <- MTPredict(B = B, X = X, task.specific.features = task.specific.features)
  } else {
    if (!isTRUE(all.equal(dim(pred), dim(Y)))) {
      stop("Dimensions of pred and Y have to coincide!")
    }
  }

  if (normalize) {
    return(mean((Y - pred)^2))
  } else {
    return(sum((Y - pred)^2))
  }
}

#' Compute (mean) correlation for linear multi-task model.
#'
#' @param Y Column centered N by K output matrix for every task.
#' @param B J by K matrix of regression coefficients, J = J1 + J2.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param pred Predicted output matrix. If NULL, compute predictions using input features.
#' @param method Correlation method to use.
#'
#' @return The mean correlation between predictions for each task and Y.
#'
#' @importFrom stats cor
#' @export
MTComputeMeanCorrelation <- function (Y, B, X = NULL, task.specific.features = list(),
                                      pred = NULL, method = "spearman") {

  if (is.null(pred)) {
    pred <- MTPredict(B = B, X = X, task.specific.features = task.specific.features)
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
#'   each task. Each entry contains an N by J2 matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param standardize Standardize data (default is TRUE).
#'
#' @return List with XTX and XTY if no task specific features are supplied, or
#'   list of lists otherwise.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}}, \code{\link{RunGroupCrossvalidation}}.
#' @export
PrepareMatrices <- function(Y, X = NULL, task.specific.features = list(),
                            standardize = TRUE) {

  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # initialize variables
  K <- ncol(Y)
  N <- nrow(Y)

  J1 <- 0
  if (!is.null(X)) {
    J1 <- ncol(X)
    if (standardize) {
      # center and scale
      X <- scale(X)
    }
  }
  J2 <- 0
  if(length(task.specific.features) > 0) {
    J2 <- ncol(task.specific.features[[1]])
    if (standardize) {
      # center and scale
      task.specific.features <- lapply(task.specific.features, scale)
    }
  }
  J <- J1 + J2

  if (standardize) {
    # center response
    Y <- scale(Y, scale = FALSE)
  }

  if (J2 == 0) {
    # no task specific features
    XTX <- t(X) %*% X
    XTY <- t(X) %*% Y
  } else {
    # task specific features
    XTX <- list ()
    XTY <- matrix(0, J, K)
    for (k in 1:K) {
      if (!is.null(X)) {
        mat <- cbind(X, task.specific.features[[k]])
      } else {
        mat <- task.specific.features[[k]]
      }
      XTX[[k]] <- t(mat) %*% mat
      XTY[, k] <- t(mat) %*% Y[, k]
    }
  }
  return(list(XTX = XTX, XTY = XTY))
}
