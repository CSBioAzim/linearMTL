#' Fit a tree-guided group lasso mixture model model (TGGLMix).
#'
#' Fit a tree-guided group lasso mixture model using a generalized EM
#' algorithm. May be trained on shared or task specific feature matrices.
#'
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#' @param M Number of Clusters.
#' @param groups Binary V by K matrix determining group membership: Task k in
#'   group v iff groups[v,k] == 1.
#' @param weights V dimensional vector with group weights.
#' @param lambda Regularization parameter.
#' @param gam (Optional) Regularization parameter for component m will be lambda
#'   times the prior for component m to the power of gam.
#' @param homoscedastic (Optional) Force variance to be the same for all tasks
#'   in a component. Default is FALSE.
#' @param EM.max.iter (Optional) Maximum number of iterations for EM algorithm.
#' @param EM.epsilon (Optional) Desired accuracy. Algorithm will terminate if
#'   change in penalized negative log-likelihood drops below EM.epsilon.
#' @param EM.verbose (Optional) Integer in {0,1,2}. verbose = 0: No output.
#'   verbose = 1: Print summary at the end of the optimization. verbose = 2:
#'   Print progress during optimization.
#' @param sample.data (Optional) Sample data according to posterior probability
#'   or not.
#' @param TGGL.mu (Optional) Mu parameter for TGGL.
#' @param TGGL.epsilon (Optional) Epsilon parameter for TGGL.
#'
#' @return List containing
#' \item{models}{List of TGGL models for each component.}
#' \item{posterior}{N by M Matrix containing posterior probabilities.}
#' \item{prior}{Vector with prior probabilities for each component.}
#' \item{sigmas}{M by K Matrix with standard deviations for each component.}
#' \item{obj}{Penalized negative log-likelihood (final objective value).}
#' \item{loglik}{Likelihood for training data.}
#'
#' @seealso \code{\link{TreeGuidedGroupLasso}}
#' @export
#' @importFrom stats runif rmultinom
TGGLMix <- function(X = NULL, task.specific.features = list(), Y, M,
                    groups, weights, lambda,
                    gam = 1, homoscedastic = FALSE,
                    EM.max.iter = 200, EM.epsilon = 1e-5,
                    EM.verbose = 0, sample.data = FALSE,
                    TGGL.mu = 1e-5, TGGL.epsilon = 1e-5) {

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

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  # penalized negative log likelihood
  obj <- Inf

  model.list <- list()
  for (m in 1:M) {
    # initialize coefficients to zero
    model.list[[m]] <- list(B = matrix(0, nrow = J, ncol = K),
                            intercept = rep(0, K))
    # predictions
    model.list[[m]]$pred <- matrix(0, nrow = N, ncol = K)
    # penalty term
    model.list[[m]]$pen <- 0
  }

  # inverse variance matrix for every component / task
  rho <- matrix(2, nrow = M, ncol = K)

  # mixture weights (priors)
  prior <- rep(1/M, M)

  # posterior probabilities for every data point / component
  tau <- matrix(0, nrow = N, ncol = M)
  tau[cbind(1:N, sample(1:M, N, replace = TRUE))] <- 0.9
  tau[tau == 0] <- 0.1
  tau <- tau / rowSums(tau)

  # delta will store the difference
  # between the objective values of two
  # successive iterations
  delta <- Inf

  # if prior probability drops below prior.min
  # abort optimization
  prior.min <- 1/N

  iter <- 1

  #####################
  # optimize using EM #
  #####################
  while (iter < EM.max.iter) {

    # save old penalized likelihood
    obj.old <- obj
    post.old <- tau
    obj <- 0

    # compute new mixing proportions
    tau.sums <- colSums(tau)
    if (sum(1/N * tau.sums < prior.min) > 0 | sum(prior < prior.min) > 0) {
      if (EM.verbose >= 1) {
        print("Component with prior 0. Abort optimization.")
      }
      # TODO add option for component removal
      likelihood <- -Inf
      break()
    }
    prior <- ComputePiMix(model.list, tau.sums, lambda, gam, prior, N)

    if (sample.data & (M > 1)) {
      # assign data point to one of the M components
      # according to its posterior distribution
      tau <- t(apply(tau, 1, FUN=function(x){rmultinom(1,1,x)}))
    }
    for (m in 1:M) {
      sqrt.tau <- sqrt(tau[, m])
      Nm <- sum(tau[, m])

      if (Nm > 0) {
        # optimize rho
        weighted.Y <- Y * sqrt.tau
        weighted.Y.norm.squared <- colSums(weighted.Y^2)
        weighted.pred <- model.list[[m]]$pred * sqrt.tau
        # inner product of weighted response and prediction
        inner.products <- colSums(weighted.Y * weighted.pred)

        if (homoscedastic) {
          scriptC <- sum(inner.products)
          scriptY <- sum(weighted.Y.norm.squared)
          rho[m, ] <- scriptC + sqrt(scriptC^2 + 4*K*Nm*scriptY)
          rho[m, ] <- rho[m, ] / (2*scriptY)
        } else {
          rho[m, ] <- inner.products + sqrt(inner.products^2 + 4*Nm*weighted.Y.norm.squared)
          rho[m, ] <- rho[m, ] / (2*weighted.Y.norm.squared)
        }

      }

      Y.rho <- sweep(Y, 2, rho[m, ], "*")

      if (Nm > 0) {
        # optimize phi
        model.list[[m]] <- TreeGuidedGroupLasso(X = X,
                                                task.specific.features = task.specific.features,
                                                Y = Y.rho,
                                                groups = groups, weights = weights,
                                                lambda = prior[m]^gam * lambda,
                                                mu = prior[m]^gam * TGGL.mu,
                                                epsilon = TGGL.epsilon,
                                                init.B = model.list[[m]]$B,
                                                verbose = 0, standardize = FALSE,
                                                row.weights = tau[, m])
        # compute new (unweighted) predictions
        pred <- MTPredict(model.list[[m]], X = X,
                          task.specific.features = task.specific.features)
        model.list[[m]]$pred <- pred
      }

      squared.resid <- (Y.rho - model.list[[m]]$pred)^2

      # if (Nm > 0) {
      #   old.sigma <- 1/rho[m, 1]
      #   if (homoscedastic) {
      #     sigma.corrected <- sqrt( sum(rowSums(squared.resid) * tau[, m]) / (K * Nm * rho[m, 1]^2))
      #   } else {
      #     sigma.corrected <- sqrt( colSums(sweep(squared.resid, 2, rho[m, ]^2, "/") * tau[, m]) / Nm)
      #   }
      #   rho[m, ] <- 1 /sigma.corrected
      #   print(sprintf("BiasCorrect: Old: %.2f -> New: %.2f", old.sigma, sigma.corrected))
      #   Y.rho <- sweep(Y, 2, rho[m, ], "*")
      #   squared.resid <- (Y.rho - model.list[[m]]$pred)^2
      # }

      # compute TGGL penalty by subtracting MSE from the returned objective value
      pen <- model.list[[m]]$obj - 1/(2*N) * sum(rowSums(squared.resid) * tau[, m])
      model.list[[m]]$pen <- pen

      # update objective by adding TGGL penalty term for component m
      obj <- obj + pen

      # compute log of unscaled tau
      tau[, m] <- rowSums(-1/2 * squared.resid) + sum(log(rho[m, ])) - K/2*log(2*pi)
      tau[, m] <- tau[, m] + log(prior[m])
    }

    # calculate new posterior probabilities
    if (sum(is.infinite(tau)) > 0) {
      if (EM.verbose > 1) {
        print("Warning (some probabilities are Inf)")
      }
      # reset posterior for data points which
      # have Inf entries for all components
      inf.rows <- which(rowSums(is.infinite(tau)) == M)
      tau[inf.rows, ] <- log(1/M)
    }
    row.maxima <- tau[cbind(1:nrow(tau), max.col(tau, "first"))]
    logsums <- row.maxima + log(rowSums(exp(tau - row.maxima)))
    tau <- exp(tau - logsums)

    # compute objective
    loglik <- sum(logsums)
    obj <- obj - 1/N * loglik

    delta <- obj.old - obj
    s <- paste(sprintf("%.2f", prior), collapse = " ")

    if (delta > EM.epsilon) {
      rejects <- 0
      if (EM.verbose > 1) {
        print(sprintf("Iter %d. Obj: %.5f. Mixing Proportions: %s", iter, obj, s))
      }
    } else {
      if (sample.data & (delta < 0)) {
        # for hard assignments, the likelihood
        # is not guaranteed to decrease.
        # use simulated annealing
        alpha <- 1 + 3*iter
        # with increasing number of iterations,
        # alpha tends to infinity and
        # accepting a worse objective becomes
        # more and more unlikely
        if (runif(1) > exp(-alpha * (-delta))) {
          # reject step, reset to
          # previous estimates
          rejects <- rejects + 1
          tau <- post.old
          obj <- obj.old
          if (rejects > 20) {
            break()
          }
        }
        if (EM.verbose > 1) {
          print(sprintf("Iter %d. Obj: %.5f. Mixing Proportions: %s", iter, obj, s))
        }
      } else {
        break()
      }
    }
    iter <- iter + 1
  }

  # compute original parameterization
  sig <- 1/rho
  for (m in 1:M) {
    model.list[[m]]$intercept <- model.list[[m]]$intercept * sig[m, ]
    model.list[[m]]$B <- sweep(model.list[[m]]$B, 2, sig[m, ], "*")
  }
  if (EM.verbose > 0) {
    print(sprintf("Total iterations: %d. Obj: %.5f. LL: %.5f. Mixing Proportions: %s",
                  iter, obj, loglik, s))
  }
  return(list(models = model.list,
              posterior = post.old,
              prior = prior,
              sigmas = sig,
              obj = obj,
              loglik = loglik))
}


ComputePiMix <- function(model.list, tau.sums, lambda, gam, prior, N) {
  # Computes a new estimate for the prior probabilities
  # using a line search
  d <- 0.1
  tau.means <- 1/N * tau.sums
  if (gam == 0) {
    # closed form optimal solution
    return(tau.means)
  } else {
    M <- length(model.list)
    ComputePiMixObj <- function(new.prior) {
      # Compute objective value for the
      # optimization of prior probabilities
      obj <- -1/N * sum(tau.sums * log(new.prior))
      for (m in 1:M) {
        # compute penalty for new.prior
        obj <- obj + (new.prior[m]/prior[m])^gam * model.list[[m]]$pen
      }
      return(obj)
    }
    obj.old <- ComputePiMixObj(prior)
    obj <- ComputePiMixObj(tau.means)
    new.prior <- tau.means
    l <- 0
    # Find largest stepsize d^l in (0, 1] such that
    # moving in the direction of (tau.means - prior)
    # does not increase the objective
    while (obj > obj.old) {
      l <- l + 1
      new.prior <- prior + d^l * (tau.means - prior)
      obj <- ComputePiMixObj(new.prior)
    }
    return(new.prior)
  }
}

#' Compute the log-likelihood for a given mixture model.
#'
#' @param models List of TGGL models.
#' @param prior Vector of prior probabilities for each component.
#' @param sigmas M by K Matrix of standard deviations.
#' @param X N by J1 matrix of features common to all tasks.
#' @param task.specific.features List of features which are specific to each
#'   task. Each entry contains an N by J2 matrix for one particular task (where
#'   columns are features). List has to be ordered according to the columns of
#'   Y.
#' @param Y N by K output matrix for every task.
#'
#' @return Log-likelihood.
#' @export
ComputeLogLikelihood <- function(models, prior, sigmas,
                                 X = NULL, task.specific.features = list(), Y) {
  loglik <- 0
  M <- length(models)
  N <- nrow(Y)
  K <- ncol(Y)

  # posterior probabilities
  tau <- matrix(0, nrow = N, ncol = M)
  for (m in 1:M) {
    pred <- MTPredict(LMTL.model = models[[m]], X = X,
                      task.specific.features = task.specific.features)
    scaled.squared.resid <- sweep((Y - pred), 2, sqrt(2)*sigmas[m, ], "/")^2
    # compute log of unscaled tau
    tau[, m] <- -rowSums(scaled.squared.resid) - sum(log(sigmas[m, ])) - K/2*log(2*pi)
    tau[, m] <- tau[, m] + log(prior[m])
  }
  row.maxima <- tau[cbind(1:nrow(tau), max.col(tau, "first"))]
  logsums <- row.maxima + log(rowSums(exp(tau - row.maxima)))
  tau <- exp(tau - logsums)
  return(list(loglik = sum(logsums), posterior = tau))
}
