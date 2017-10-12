#' Fit a tree-guided group lasso model (Kim et al. 2010).
#'
#' Fit a tree-guided group lasso model via smoothed proximal gradient descent
#' (Kim et al. 2012). May be trained on shared or task specific feature
#' matrices.
#'
#' @param X Column-centered N by J1 matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param Y Column centered N by K output matrix for every task.
#' @param groups V by K matrix determining group membership: Task k in group v
#'   iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param lambda Regularization parameter.
#' @param max.iter Maximum number of iterations.
#' @param epsilon Desired accuracy. If error change drops below epsilon, the
#'   algorithm terminates. Also used to compute mu if mu == NULL.
#' @param mu Determines accuracy of smooth approximation to the group penalty
#'   term.
#' @param mu.adapt Multiply mu with a factor of mu.adapt every iteration.
#' @param XTX Precomputed matrices t(X)*X as for example produced by
#'   PrepareMatrices.
#' @param XTY Precomputed matrices t(X)*Y as for example produced by
#'   PrepareMatrices
#' @param init.beta J by K matrix with initializations for the regression
#'   coefficients.
#' @param verbose Integer in {0,1,2}. verbose = 0: No output. verbose = 1: Print
#'   summary at the end of the optimization. verbose = 2: Print progress during
#'   optimization.
#'
#' @return List containing
#'   \item{lambda}{Regularization parameter used.}
#'   \item{weights}{Node weights used.}
#'   \item{beta}{Final estimate of the regression coefficients.}
#'   \item{obj}{Final objective value.}
#'   \item{early.termination}{Boolean indicating whether the algorithm exceeded
#'   max.iter iterations.}
#'
#' @seealso \code{\link{RunGroupCrossvalidation}}
#' @export
TreeGuidedGroupLasso <- function (X = NULL, task.specific.features = list(), Y,
                                  groups, weights, lambda, max.iter = 10000, epsilon = 1e-5,
                                  mu = NULL, mu.adapt = 1, XTX = NULL, XTY = NULL, init.beta = NULL, verbose = 1) {

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

  if (ncol(Y) != ncol(groups)) {
    stop("X and groups must have the same number of columns!")
  }
  if (length(weights) != nrow(groups)) {
    stop("Length of weights has to equal number of rows of groups!")
  }

  if (!is.null(XTX) | !is.null(XTY)) {
    if (xor(!is.null(XTX), !is.null(XTY))) {
      stop("When using immutables, both XTX and XTY need to be specified!")
    }
    use.cached.immutables <- TRUE
  } else {
    use.cached.immutables <- FALSE
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  # isolate singleton groups
  singleton.groups <- rowSums(groups) == 1
  singleton.tasks <- apply(groups[singleton.groups,], 1, FUN=function(x){which(x != 0)})
  singleton.weights <- rep(0, K)
  singleton.weights[singleton.tasks] <- weights[singleton.groups]

  # retain inner nodes
  groups <- groups[!singleton.groups, ]
  weights <- weights[!singleton.groups]
  V <- nrow(groups) # number of inner nodes

  # precompute matrices if necessary
  if (!use.cached.immutables) {
    mats <- PrepareMatrices(Y = Y, X = X, task.specific.features = task.specific.features)
    XTX <- mats$XTX
    XTY <- mats$XTY
  }

  # build C
  group.sizes <- rowSums(groups)
  group.sums <- sum(groups)
  group.cum.sizes <- cumsum(group.sizes)
  group.ranges <- cbind(c(1, group.cum.sizes[1:V-1] + 1), group.cum.sizes, group.sizes)

  J.idx <- rep(0, group.sums)
  W.entries <- rep(0, group.sums)
  for (v in 1:V) {
    J.idx[group.ranges[v,1]:group.ranges[v,2]] <- which(groups[v,] != 0)
    W.entries[group.ranges[v,1]:group.ranges[v,2]] <- weights[v]
  }
  C <- matrix(0, nrow = group.sums, ncol = K)
  C[cbind(1:group.sums, J.idx)] <- W.entries
  C <- C * lambda

  # set mu
  D <- J * V / 2
  if (is.null(mu)) {
    mu <- epsilon / (2*D)
  }

  # compute Lipschitz constant
  C.norm.squared <- t(tcrossprod(rep(1, K), weights)) * groups;
  C.norm.squared <- max(colSums(C.norm.squared^2));
  if (length(task.specific.features) > 0) {
    # expensive
    #L1 <- max(unlist(lapply(XTX, FUN = function(M) {max(eigen(M)$values)})))
    L1 <- max(eigen(XTX[[1]])$values)
  } else {
    L1 <- max(eigen(XTX)$values)
  }

  # Lipschitz constant
  L <- L1 + C.norm.squared / mu

  # initialize variables
  if (is.null(init.beta)) {
    B <- matrix(0, nrow = J, ncol = K)
  } else {
    B <- init.beta
  }
  W <- B
  iter <- 0
  theta <- 1
  delta <- epsilon + 1
  obj.old <- ComputeObjective(Y = Y, beta = B, X = X,
                              task.specific.features = task.specific.features,
                              C = C, group.ranges = group.ranges)
  dh <- matrix(0, nrow = J, ncol = K)
  start.time <- Sys.time();
  while (iter < max.iter & delta > epsilon) {
    # optimize
    A <- Shrink(C %*% t(W) / mu, group.ranges)
    df <- t(A) %*% C
    if (length(task.specific.features) > 0) {
      # task specific features
      for (k in 1:K) {
        dh[, k] <- XTX[[k]] %*% W[, k] - XTY[[k]]
      }
      dh <- dh + df
    } else {
      # no task specific features
      dh <- XTX %*% W - XTY + df
    }
    V <- W - 1/L * dh
    B.new  <- sign(V) * pmax(0, sweep(abs(V), 2, lambda/L * singleton.weights))
    theta.new <- 2 / (iter + 3)
    W <- B.new + (1 - theta) / theta * theta.new * (B.new - B)

    # compute new objective
    obj <- ComputeObjective(Y = Y, beta = B.new, X = X,
                            task.specific.features = task.specific.features,
                            C = C, group.ranges = group.ranges)
    delta <- abs(obj - obj.old) / abs(obj.old)
    obj.old <- obj
    if (verbose == 2) {
      print(sprintf("Iter: %d, Obj: %f, delta: %f", iter, obj,delta))
    }

    # adapt learning rate
    mu <- mu * mu.adapt
    L <- L1 + C.norm.squared / mu

    B <- B.new
    theta <- theta.new
    iter <- iter + 1
  }
  end.time <- Sys.time()
  if (verbose >= 1) {
    print(sprintf("Total: Iter: %d, Obj: %f, Time: %f", iter, obj, as.numeric(end.time - start.time, units = "secs")))
  }

  early.termination <- iter > max.iter
  if (early.termination) {
    print(sprintf("Warning: Reached maximum number of iterations (%d).", max.iter))
  }
  return(list(lambda = lambda, weights = weights, beta = B, obj = obj, early.termination = early.termination))
}

Shrink <- function(A, group.ranges) {
  # Applies shrinkage operator to matrix A.
  V <- nrow(group.ranges)
  for (v in 1:V) {
    idx <- group.ranges[v, 1]:group.ranges[v, 2]
    gnorm <- sqrt(colSums(A[idx,]^2))
    A[idx, ] <- sweep(A[idx, ], 2, pmax(gnorm , 1), FUN = "/")
  }
  return(A)
}

CalculateInnerGroupPenalty <- function(A, group.ranges) {
  # Computes the penalty term for non-singleton groups.
  V <- nrow(group.ranges)
  s <- 0
  for (v in 1:V) {
    idx <- group.ranges[v, 1]:group.ranges[v, 2]
    gnorm <- sqrt(colSums(A[idx,]^2))
    s <- s + sum(gnorm)
  }
  return(s)
}

ComputeObjective <- function(Y, beta, X = NULL, task.specific.features = list(), C, group.ranges) {
  # Compute the optimization objective
  obj <- 1/2 * MTComputeError(Y = Y, beta = beta, X = X,
                              task.specific.features = task.specific.features,
                              normalize = FALSE)
  obj <- obj + CalculateInnerGroupPenalty(C %*% t(beta), group.ranges)
  return(obj)
}


#' Construct tree for \code{\link{TreeGuidedGroupLasso}} using hierarchical
#' clustering.
#'
#' Compute tree and node weights from the results of applying hierarchical
#' clustering to the (response) matrix Y.
#'
#' @param Y N by K matrix for every task.
#' @param threshold Threshold used for cutting nodes with high dissimilarity.
#' @param linkage Linkage type used for hierarchical clustering.
#'
#' @return List containing \item{groups}{V by K matrix determining group
#'   membership: Task k in group v iff groups[v,k] == 1.} \item{weights}{V
#'   dimensional weight vector.}
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}},
#'   \code{\link{ComputeGroupWeights}}.
#' @export
BuildTreeHC <- function(Y, threshold = 0.9, linkage = "complete") {

  K <- ncol(Y)
  hc <- hclust(as.dist(1 - cor(Y)), method = linkage)

  # (equivalent to matlab hc representation)
  clusters <- hc$merge
  clusters[clusters > 0] <- clusters[clusters > 0] + K
  clusters <- abs(clusters)

  s <- hc$height / max(hc$height)
  V <- K + sum(s < threshold)

  I.idx <- 1:K
  J.idx <- 1:K
  for (i in 1:(V-K)) {
    leaves <- FindLeaves(clusters, i, K, leaves = c())
    I.idx <- c(I.idx, rep(i + K, length(leaves)))
    J.idx <- c(J.idx, leaves)
  }
  groups <- matrix(0, nrow = V, ncol = K)
  groups[cbind(I.idx, J.idx)] <- rep(1, length(I.idx))

  # compute node weights
  weights <- ComputeGroupWeights(groups, c(rep(1, K), s[s < threshold]))
  return(list(groups = groups, weights = weights))
}

FindLeaves <- function(clusters, i, K, leaves) {
  # Used by BuildTreeHC to find tasks in cluster i.
  for (c in clusters[i, ]) {
    if (c > K) {
      leaves <- FindLeaves(clusters, c - K, K, leaves)
    } else {
      leaves <- c(leaves, c)
    }
  }
  return(leaves)
}




#' Compute node weights using scheme from Kim et al [2010].
#'
#' @param groups V by K matrix determining group membership: Task k in group v iff groups[v,k] == 1.
#' @param s V dimensional vector with weights. Entries for singleton groups will be ignored.
#'
#' @return Vector containing node weights.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}},
#'   \code{\link{BuildTreeHC}}.
#' @export
ComputeGroupWeights <- function(groups, s) {

  V <- nrow(groups)
  K <- ncol(groups)

  if (sum(rowSums(groups) == 1) != K) {
    stop("groups does not contain K singletons.")
  }

  # make sure that singleton groups are at position 1:K
  group.size.ordering <- order(rowSums(groups))
  groups <- groups[group.size.ordering, ]
  s <- s[group.size.ordering]

  weights <- rep(1, V)
  for (i in (K+1):V) {
    weights[i] <- weights[i] * (1 - s[i])
    leaves <- which(groups[i, ] == 1)
    R <- rowSums(groups[, leaves])
    # assuming that groups encodes a tree, two groups are either disjoint
    # or one is contained in the other. The first condition excludes
    # group i and its ancestors; the second condition includes children of group i.
    children <- which((R < length(leaves)) & (R >= 1))
    weights[children] <- weights[children] * s[i]
  }

  # return weights in original order of groups matrix
  return(weights[order(group.size.ordering)])
}
