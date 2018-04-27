#' Fit a tree-guided group lasso model (Kim et al. 2010).
#'
#' Fit a tree-guided group lasso model via smoothed proximal gradient descent
#' (Kim et al. 2012). May be trained on shared or task specific feature
#' matrices.
#'
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param groups Binary V by K matrix determining group membership: Task k in
#'   group v iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param lambda Regularization parameter.
#' @param max.iter (Optional) Maximum number of iterations.
#' @param epsilon (Optional) Desired accuracy. If error change drops below
#'   epsilon, the algorithm terminates.
#' @param mu (Optional) Determines accuracy of smooth approximation to the group
#'   penalty term. If NULL, mu will be determined based on desired accuracy
#'   epsilon. However, this may lead to numerical issues and it is thus
#'   recommended to tune mu by hand.
#' @param mu.adapt (Optional) Multiply mu with a factor of mu.adapt every
#'   iteration. Default is no adaptation (mu.adapt = 1).
#' @param cached.mats Precomputed matrices as produced by PrepareMatrices.
#' @param MSE.Lipschitz (Optional) Lipschitz constant for MSE. If missing, use
#'   maximum Eigenvalue of XTX.
#' @param init.B (Optional) (J+1) by K matrix with initializations for the
#'   regression coefficients and intercept.
#' @param verbose (Optional) Integer in {0,1,2}. verbose = 0: No output. verbose
#'   = 1: Print summary at the end of the optimization. verbose = 2: Print
#'   progress during optimization.
#' @param standardize (Optional) Default is TRUE. Standardize data (using R
#'   function scale()). Coefficients will be returned on original scale.
#' @param fit.intercept (Optional) Default is TRUE. Include intercept.
#' @param row.weights (Optional) Use weighted MSE. When cached.mats is supplied,
#'   it is assumed that the rows of X and Y were already weighted appropriately.
#'
#' @return List containing
#' \item{lambda}{Regularization parameter used.}
#' \item{weights}{Node weights used.}
#' \item{B}{Final estimate of the regression coefficients.}
#' \item{intercept}{Final estimate of the intercept terms.}
#' \item{obj}{Final objective value.}
#' \item{early.termination}{Boolean indicating whether the algorithm exceeded
#' max.iter iterations.}
#'
#' @seealso \code{\link{RunGroupCrossvalidation}}
#' @export
TreeGuidedGroupLasso <- function (X = NULL, task.specific.features = list(), Y,
                                  groups, weights, lambda,
                                  max.iter = 10000, epsilon = 1e-5,
                                  mu = NULL, mu.adapt = 1,
                                  cached.mats = NULL,
                                  MSE.Lipschitz = NULL,
                                  init.B = NULL,
                                  verbose = 1,
                                  standardize = TRUE,
                                  fit.intercept = TRUE,
                                  row.weights = NULL) {

  # This implementation minimizes the objective
  #   1/(2N)||Y - XB - intercept|| + lambda*Omega(B)
  # where Omega denotes the tree-guided group penalty.
  # Variable names are mainly based on the original paper.

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

  # check for precomputed matrices
  if (!is.null(cached.mats)) {
    use.cached.immutables <- TRUE
  } else {
    use.cached.immutables <- FALSE
  }

  if (standardize) {
    # we will recover the intercept at the end
    fit.intercept <- FALSE
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  ####################################################
  # transform inputs and outputs if necessary;       #
  # precompute matrix products                       #
  ####################################################

  if (!use.cached.immutables) {
    cached.mats <- PrepareMatrices(Y = Y, X = X,
                                   task.specific.features = task.specific.features,
                                   standardize = standardize, row.weights = row.weights,
                                   return.inputs = TRUE)
    X <- cached.mats$X
    task.specific.features <- cached.mats$task.specific.features
    Y <- cached.mats$Y
  }
  XTX <- cached.mats$XTX
  XTY <- cached.mats$XTY
  XT1 <- cached.mats$XT1
  YT1 <- cached.mats$YT1

  if (is.null(row.weights)) {
    row.weights <- rep(1, K)
  }
  row.weights.sqrt <- sqrt(row.weights)
  row.weights.sum <- sum(row.weights)

  ###################
  # initializations #
  ###################

  # determine groups consisting of only one task
  singleton.groups <- rowSums(groups) == 1
  singleton.tasks <- apply(groups[singleton.groups,], 1, FUN=function(x){which(x != 0)})
  singleton.weights <- rep(0, K)
  singleton.weights[singleton.tasks] <- weights[singleton.groups]

  # only keep inner nodes in groups matrix and weights
  groups <- groups[!singleton.groups, , drop = FALSE]
  weights <- weights[!singleton.groups]
  V <- nrow(groups) # number of inner nodes

  # build matrix C
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

  # compute squared norm of C
  C.norm.squared <- t(tcrossprod(rep(1, K), weights)) * groups;
  C.norm.squared <- lambda^2 * max(colSums(C.norm.squared^2));

  # compute mu if necessary
  if (is.null(mu)) {
    D <- J * V / 2
    mu <- epsilon / (2*D)
  }

  # compute Lipschitz constant
  if (is.null(MSE.Lipschitz)) {
    if (J2 > 0) {
      # too expensive:
      # L1 <- max(unlist(lapply(XTX, FUN = function(M){max(eigen(M)$values)})))
      # instead just consider first matrix:

      # divide by N because of 1/(2N) factor in the objective
      L1 <- max(eigen(XTX[[1]])$values) / N
    } else {
      L1 <- max(eigen(XTX)$values) / N
    }
  } else {
    L1 <- MSE.Lipschitz / N
  }
  # Lipschitz constant of the smooth part of the objective function
  L <- L1 + C.norm.squared / mu

  # initialize variables
  if (is.null(init.B)) {
    B <- matrix(0, nrow = J, ncol = K)
  } else {
    B <- init.B
  }

  intercept <- rep(0, K)

  W <- B
  iter <- 0
  theta <- 1
  delta <- epsilon + 1

  # compute initial objective function
  obj.old <- ComputeObjective(Y = Y, B = B, intercept = intercept, X = X,
                              task.specific.features = task.specific.features,
                              C = C, group.ranges = group.ranges, lambda = lambda,
                              singleton.weights = singleton.weights)
  dh <- matrix(0, nrow = J, ncol = K)
  start.time <- Sys.time();

  early.termination <- FALSE
  ################################################
  # optimization using proximal gradient descent #
  ################################################
  while (iter < max.iter & delta > epsilon) {
    # apply shrinkage operator
    A <- RcppShrink(C %*% t(W) / mu, group.ranges)

    # gradient of the smooth approximation to the nonsmooth penalty
    df <- t(A) %*% C

    # compute gradient of the smooth part of the objective
    if (length(task.specific.features) > 0) {
      # task specific features
      for (k in 1:K) {
        dh[, k] <- 1/N * (XTX[[k]] %*% W[, k] - XTY[, k] + XT1[[k]] * intercept[k])
      }
    } else {
      # no task specific features
      dh <- 1/N * (XTX %*% W - XTY + XT1 %*% t(intercept))
    }
    dh <- dh + df

    # Q corresponds to the matrix V in the original paper
    Q <- W - 1/L * dh

    # apply soft-thresholding
    B.new  <- sign(Q) * pmax(0, sweep(abs(Q), 2, lambda/L * singleton.weights))

    theta.new <- 2 / (iter + 3)
    W <- B.new + (1 - theta) / theta * theta.new * (B.new - B)

    if (fit.intercept) {
      intercept.new <- ComputeIntercept(XT1, YT1, B.new, row.weights.sum)
    } else {
      intercept.new <- rep(0, K)
    }

    # compute delta
    delta <- sqrt(sum((B.new - B)^2) + sum((intercept.new - intercept)^2))
    delta <- delta / max(sqrt(sum(B.new^2) + sum(intercept.new^2)), 1)

    if (is.nan(delta)) {
      print("Premature stopping due to Inf parameter values. Consider retuning mu.")
      early.termination <- TRUE
      break()
    }

    if (((iter %% 100) == 0) & (iter > 0)) {
      # compute new objective
      obj <- ComputeObjective(Y = Y, B = B.new, intercept = intercept.new, X = X,
                              task.specific.features = task.specific.features,
                              C = C, group.ranges = group.ranges, lambda = lambda,
                              singleton.weights = singleton.weights)
      delta <- min(delta, abs(obj - obj.old) / abs(obj.old))
      obj.old <- obj
      if (verbose == 2) {
        print(sprintf("Iter: %d, Obj: %f, delta: %f", iter, obj, delta))
      }
    }

    # adapt mu
    mu <- mu * mu.adapt
    L <- L1 + C.norm.squared / mu

    B <- B.new
    intercept <- intercept.new
    theta <- theta.new
    iter <- iter + 1
  }
  end.time <- Sys.time()

  # compute objective
  obj <- ComputeObjective(Y = Y, B = B, intercept = intercept, X = X,
                          task.specific.features = task.specific.features,
                          C = C, group.ranges = group.ranges, lambda = lambda,
                          singleton.weights = singleton.weights)
  if (verbose >= 1) {
    print(sprintf("Total: Iter: %d, Obj: %f, Time: %f", iter, obj,
                  as.numeric(end.time - start.time, units = "secs")))
  }

  early.termination <- (iter > max.iter) | early.termination
  if (early.termination & (verbose > 1)) {
    print(sprintf("Warning: Reached maximum number of iterations (%d).", max.iter))
  }

  if (standardize) {
    # recover coefficients on original scale and compute intercept
    B.transformed <- ComputeInverseArtesi(B, cached.mats$sample.moments)
    B <- B.transformed[2:nrow(B.transformed), ]
    intercept <- B.transformed[1, ]
  }

  return(list(lambda = lambda, weights = weights,
              B = B, intercept = intercept,
              obj = obj, early.termination = early.termination))
}

# Use Rcpp implementation instead: 'RcppShrink'.
Shrink <- function(A, group.ranges) {
  # Applies shrinkage operator to matrix A.
  #
  # A is matrix of size (cumsum(inner group sizes)) by J.
  # Shrink all vectors A[group ids, j].
  V <- nrow(group.ranges)
  for (v in 1:V) {
    idx <- group.ranges[v, 1]:group.ranges[v, 2]
    gnorm <- sqrt(colSums(A[idx, ,drop = FALSE]^2))
    A[idx, ] <- sweep(A[idx, ,drop = FALSE], 2,
                      pmax(gnorm , 1), FUN = "/")
  }
  return(A)
}

CalculateInnerGroupPenalty <- function(A, group.ranges) {
  # Computes the penalty term for non-singleton groups.
  V <- nrow(group.ranges)
  s <- 0
  for (v in 1:V) {
    idx <- group.ranges[v, 1]:group.ranges[v, 2]
    gnorm <- sqrt(colSums(A[idx, ,drop = FALSE]^2))
    s <- s + sum(gnorm)
  }
  return(s)
}

ComputeObjective <- function(Y, B, intercept = intercept,
                             X = NULL, task.specific.features = list(),
                             C, group.ranges, lambda, singleton.weights) {
  # Computes the optimization objective.
  N <- nrow(Y)
  K <- ncol(Y)
  # compute MSE
  pred <- MTPredict(LMTL.model = list(B = B, intercept = intercept), X = X,
                    task.specific.features = task.specific.features)
  obj <- 1/(2*N) * sum((Y - pred)^2)
  # compute penalty for inner nodes
  obj <- obj + CalculateInnerGroupPenalty(C %*% t(B), group.ranges)
  # compute penalty for leaf nodes (singletons)
  obj <- obj + sum(sweep(abs(B), 2, lambda * singleton.weights, FUN = "*"))

  return(obj)
}

ComputeIntercept <- function(XT1, YT1, B.new, row.weights.sum) {
  K <- ncol(B.new)
  intercept <- rep(0, K)
  if (is.list(XT1)) {
    for (k in 1:K) {
      intercept[k] <- YT1[k] - t(XT1[[k]]) %*% B.new[, k]
    }
  } else {
    intercept <- t(YT1) - t(XT1) %*% B.new
  }
  return(as.vector(intercept/row.weights.sum))
}

ComputeInverseArtesi <- function(B, sample.moments) {
  # Compute inverse artesi transformation to obtain coefficients for
  # variables on the original scale.
  B.transformed <- matrix(0, nrow = nrow(B) + 1, ncol = ncol(B))
  K <- ncol(B)

  if (is.null(sample.moments$Y.sds)) {
    Y.sds <- rep(1, length(sample.moments$Y.sds))
  }

  for (k in 1:K) {
    input.sds <- sample.moments$X.sds
    input.mus <- sample.moments$X.mus
    if (length(sample.moments$tsf.sds) > 0) {
      input.sds <- c(input.sds, sample.moments$tsf.sds[[k]])
      input.mus <- c(input.mus, sample.moments$tsf.mus[[k]])
    }
    B.transformed[2:nrow(B.transformed), k] <- sample.moments$Y.sds[k] * B[, k] / input.sds
    B.transformed[1, k] <-  - sample.moments$Y.sds[k] * sum(B[, k] * input.mus / input.sds) + sample.moments$Y.mus[k]
  }
  return(B.transformed)
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
#' @importFrom stats as.dist cor hclust
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

#' Construct tree for \code{\link{TreeGuidedGroupLasso}} from prior grouping.
#'
#' Generalize hierarchical clustering distances to compute node weights from the
#' (response) matrix Y for a given grouping. Weights will be normalized using
#' ComputeGroupWeights.
#'
#' @param Y N by K matrix for every task.
#' @param groups V by K matrix determining group membership: Task k in group v
#'   iff groups[v,k] == 1.
#' @param linkage "complete" or "average". Linkage type used (distance measure
#'   between clusters).
#' @param distance Function for computing pairwise distances between columns of
#'   Y. Default is 1-cor(a,b).
#' @return Vector containing node weights.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}},
#'   \code{\link{ComputeGroupWeights}}.
#' @export
ComputeWeightsFromTree <- function(Y, groups,
                                   linkage = "complete",
                                   distance = function(a,b){1-cor(a,b)}) {

  V <- nrow(groups)
  K <- ncol(groups)

  if (sum(rowSums(groups) == 1) != K) {
    stop("groups does not contain K singletons!")
  }

  # make sure that singleton groups are at position 1:K
  group.size.ordering <- order(rowSums(groups))
  groups <- groups[group.size.ordering, ]

  s <- rep(0, V)
  for (i in (K+1):V) {
    # determine children of node i
    node <- groups[i, ]
    node.size <- sum(node)
    descendants <- groups[1:i, ]
    descendants <- descendants[rowSums(descendants[, as.logical(node)]) > 0, ]
    children <- c()
    l <- nrow(descendants)
    # descendants is ordered according to group size
    while (sum(children) < node.size) {
      l <- l - 1
      if (l < 1) {
        stop("Invalid tree. Does tree contain all singletons?")
      }
      c.temp <- rbind(descendants[l, ], children)
      if (!any(colSums(c.temp) > 1)) {
        children <- c.temp
      }
    }
    s[i] <- ComputeDistance(children, Y, linkage, distance)
  }
  if (max(s) > 1) {
    s <- s / max(s)
  }

  weights <- ComputeGroupWeights(groups, s)
  # return weights in original order of groups matrix
  return(weights[order(group.size.ordering)])
}

ComputeDistance <- function(children, Y, linkage, distance) {
  # Used by ComputeWeightsFromTree to determine distance between clusters in
  # children.
  pairwise.distances <- c()
  for (c1 in 1:(nrow(children)-1)) {
    for (c2 in (c1+1):nrow(children)) {
      group1 <- which(children[c1, ] == 1)
      group2 <- which(children[c2, ] == 1)
      for (t1 in group1) {
        for (t2 in group2) {
          pairwise.distances <- c(pairwise.distances, distance(Y[, t1], Y[, t2]))
        }
      }
    }
  }
  if (linkage == "average") {
    return(mean(pairwise.distances))
  } else if (linkage == "complete") {
    return(max(pairwise.distances))
  } else {
    stop("Unknown linkage.")
  }
}

#' Compute node weights using scheme from Kim et al [2010].
#'
#' @param groups V by K matrix determining group membership: Task k in group v
#'   iff groups[v,k] == 1.
#' @param s V dimensional vector with weights. Entries for singleton groups will
#'   be ignored.
#'
#' @return Vector containing node weights.
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}}, \code{\link{BuildTreeHC}}.
#' @export
ComputeGroupWeights <- function(groups, s) {

  V <- nrow(groups)
  K <- ncol(groups)

  if (length(s) != V) {
    stop("s does not contain one entry for every group!")
  }

  if (sum(rowSums(groups) == 1) != K) {
    stop("groups does not contain K singletons!")
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
