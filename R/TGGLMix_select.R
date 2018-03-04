#' TGGLMixSelection is a convenience function for which calls RunTGGLMix for all
#' possible combinations of parameters in M.vec and lambda.vec.
#'
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param M.vec Vector with number of clusters to be tested.
#' @param lambda.vec Vector with lambda values to be tested.
#' @param test.ids Mandatory: Indices of data points to be used for testing.
#' @param groups Binary V by K matrix determining group membership: Task k in
#'   group v iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param verbose (Optional) Integer in {0,1,2}. verbose = 0: No output.
#'   verbose = 1: Print summary at the end of the optimization. verbose = 2:
#'   Print progress during optimization.
#' @param ... Additional parameters passed to
#'   \code{\link{RunTGGLMix}}.
#'
#' @return List containing
#' \item{results}{List of TGGLMix models for each parameter setting.}
#' \item{top.model}{Model with highest predictive likelihood.}
#' @export
TGGLMixSelection <- function(X = NULL, task.specific.features = list(), Y, M.vec, lambda.vec,
                             test.ids, groups, weights, verbose = 0, ...) {

  # build parameter grid
  parameter.grid <- cbind(rep(M.vec, each = length(lambda.vec)),
                          rep(lambda.vec, times = length(M.vec)))

  if (verbose > 0) {
    print(sprintf("Running TGGL-Mixture for %d parameter settings ... ", nrow(parameter.grid)))
  }
  sel.start.time <- Sys.time()
  results <- list()
  for (ind in 1:nrow(parameter.grid)) {
    M <- parameter.grid[ind, 1]
    lambda <- parameter.grid[ind, 2]
    results[[ind]] <- RunTGGLMix(X = X, task.specific.features = task.specific.features, Y = Y, M = M,
                               test.ids = test.ids, groups = groups, weights = weights,
                               lambda = lambda, verbose = verbose, ...)
  }
  opt.model.idx <- which.max(sapply(results, FUN = function(x){x$test.loglik}))
  sel.end.time <- Sys.time()
  sel.time <- as.numeric(sel.end.time - sel.start.time, units = "mins")
  if (verbose > 0) {
    print(sprintf("Minutes to run Model Selection %0.1f. Best Parameter Setting: M = %d, lambda = %.e.",
                  sel.time, parameter.grid[opt.model.idx, 1], parameter.grid[opt.model.idx, 2]))
  }
  return(list(results = results, top.model = results[[opt.model.idx]]))
}
