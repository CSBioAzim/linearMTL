################################################################################
# Utility functions for linear multi-task models.
################################################################################

#' Compute predictions for linear multi-task model.
#'
#' Compute linear predictions from predictor variables and regression
#' coefficients.
#'
#' @param beta J by K matrix of regression coefficients, J = J1 + J2.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#'
#' @return N by J predictions matrix.
#'
#' @export
MTPredict <- function(beta, X = NULL, task.specific.features = list()) {

  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  K <- ncol(beta)
  # compute predictions for shared feature matrix
  if (!is.null(X)) {
    J1 <- ncol(X)
    # first J1 rows of beta correspond to features in X
    prediction <- X %*% beta[1:J1, ]
  } else {
    J1 <- 0
    C <- task.specific.features[[1]]
    prediction <- matrix(0, nrow = nrow(C), ncol = ncol(beta))
  }

  # compute predictions for task specific feature matrices
  if (length(task.specific.features) > 0) {
    # rows J1 + 1 and below of beta correspond to task specific features
    tsf.idx <- (J1 + 1):nrow(beta)
    for (k in 1:K) {
      C <- task.specific.features[[k]]
      prediction[, k] <- prediction[, k] + C %*% beta[tsf.idx, k]
    }
  }
  return(prediction)
}

#' Compute (mean) squared error for linear multi-task model.
#'
#' @param Y Column centered N by K output matrix for every task.
#' @param beta J by K matrix of regression coefficients, J = J1 + J2.
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
MTComputeError <- function (Y, beta, X = NULL, task.specific.features = list(),
                            pred = NULL, normalize = TRUE) {

  if (is.null(pred)) {
    pred <- MTPredict(beta = beta, X = X, task.specific.features = task.specific.features)
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
#' @param beta J by K matrix of regression coefficients, J = J1 + J2.
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
MTComputeMeanCorrelation <- function (Y, beta, X = NULL, task.specific.features = list(),
                                      pred = NULL, method = "spearman") {

  if (is.null(pred)) {
    pred <- MTPredict(beta = beta, X = X, task.specific.features = task.specific.features)
  }
  return(mean(diag(cor(pred, Y, method = method))))
}




#' Precompute matrices XTX and XTY for \code{\link{TreeGuidedGroupLasso}}.
#'
#' Computes XTX and XTY. If task specific features are supplied, the result will
#' be a list containing the respective matrices for each task.
#'
#' @param Y Column centered N by K output matrix for every task.
#' @param X Column-centered N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param idx Vector of indices (data points) to exclude.
#'
#' @return List with XTX and XTY if no task specific features are supplied, or
#'   list of lists otherwise.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}}, \code{\link{RunGroupCrossvalidation}}.
#' @export
PrepareMatrices <- function(Y, X = NULL, task.specific.features = list(), idx = NULL) {

  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # initialize variables
  K <- ncol(Y)
  N <- nrow(Y)
  J1 <- if (is.null(X)) 0 else nrow(X)
  J2 <- if(length(task.specific.features) == 0) 0 else ncol(task.specific.features[[1]])
  J <- J1 + J2

  if (is.null(idx)) {
    idx <- -seq_len(N)
  }

  if (J2 == 0) {
    # no task specific features
    XTX <- t(X[-idx, ]) %*% X[-idx, ]
    XTY <- t(X[-idx, ]) %*% Y[-idx, , drop = FALSE]
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
      XTX[[k]] <- t(mat[-idx, ]) %*% mat[-idx, ]
      XTY[, k] <- t(mat[-idx, ]) %*% Y[-idx, k]
    }
  }
  return (list(XTX = XTX, XTY = XTY))
}
