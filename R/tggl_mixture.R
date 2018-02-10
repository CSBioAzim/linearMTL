TGGLMixture <- function(X = NULL, task.specific.features = list(), Y, M,
                        groups, weights, lambda, gam = 0,
                        max.iter = 200, eps = 1e-5, verbose = -1, sample.data = TRUE) {

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
  pen.negloglik <- Inf

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

  iter <- 1

  #####################
  # optimize using EM #
  #####################
  while (iter < max.iter) {

    # save old penalized likelihood
    old.lik <- pen.negloglik

    pen.negloglik <- 0
    post.old <- tau
    if (sample.data) {
      # assign data point to one of the M components
      # according to its posterior distribution
      tau <- t(apply(tau, 1, FUN=function(x){rmultinom(1,1,x)}))
    }
    # compute new mixing proportions
    prior <- ComputePiMix(model.list, tau, lambda, gam, prior)
    for (m in 1:M) {
      sqrt.tau <- sqrt(tau[, m])
      Nm <- sum(tau[, m])

      # optimize rho
      weighted.Y <- Y * sqrt.tau
      weighted.Y.norm.sqrt <- colSums(weighted.Y^2)
      weighted.pred <- model.list[[m]]$pred * sqrt.tau
      # inner product of weighted response and prediction
      inner.products <- colSums(weighted.Y * weighted.pred)
      rho[m, ] <- inner.products + sqrt(inner.products^2 + 4*Nm*weighted.Y.norm.sqrt)
      rho[m, ] <- rho[m, ] / (2*weighted.Y.norm.sqrt)

      # optimize phi
      Yrho <- sweep(Y, 2, rho[m, ], "*")
      model.list[[m]] <- TreeGuidedGroupLasso(X = X,
                                              task.specific.features = task.specific.features,
                                              Y = Yrho,
                                              groups = groups, weights = weights,
                                              lambda = prior[m]^gam * lambda,
                                              max.iter = 10000, epsilon = 1e-5,
                                              mu = 1e-5, mu.adapt = 0.99,
                                              init.B = model.list[[m]]$B,
                                              verbose = 0, standardize = FALSE,
                                              row.weights = tau[, m])
      # compute (unweighted) predictions
      pred <- MTPredict(model.list[[m]], X = X,
                        task.specific.features = task.specific.features)
      model.list[[m]]$pred <- pred
      squared.resid <- (Yrho - pred)^2

      # compute TGGL penalty by subtracting MSE from the returned objective value
      pen <- model.list[[m]]$obj - 1/(2*N) * sum(squared.resid*tau[, m])
      model.list[[m]]$pen <- pen

      # update objective by adding TGGL penalty term for component m
      pen.negloglik <- pen.negloglik + pen

      # compute log of unscaled tau
      tau[, m] <- rowSums(-1/2 * squared.resid) + sum(log(rho[m, ])) - K/2*log(2*pi)
      tau[, m] <- tau[, m] + log(prior[m])
    }

    # calculate new posterior probabilities
    if (sum(is.infinite(tau)) > 0) {
      if (verbose > 1) {
        print("Warning (some probabilities are Inf)")
      }
      # reset posterior for data points which
      # have Inf entries for all components
      inf.rows <- which(rowSums(is.infinite(R)) == M)
      tau[inf.rows, ] <- log(1/M)
    }
    row.maxima <- tau[cbind(1:nrow(tau), max.col(tau, "first"))]
    logsums <- row.maxima + log(rowSums(exp(tau - row.maxima)))
    tau <- exp(tau - logsums)

    # compute objective
    likelihood <- sum(logsums)
    pen.negloglik <- pen.negloglik - 1/N * likelihood

    delta <- old.lik - pen.negloglik
    s <- paste(sprintf("%.2f", prior), collapse = " ")
    if (delta > eps) {
      if (verbose >= 1) {
        print(sprintf("Iter %d. PenNegLog: %.5f. Mixing Proportions: %s", iter, pen.negloglik, s))
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
          tau <- post.old
          pen.negloglik <- old.lik
        }
        if (verbose >= 1) {
          print(sprintf("Iter %d. PenNegLog: %.5f. Mixing Proportions: %s", iter, pen.negloglik, s))
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
    model.list[[m]] <- rbind(model.list[[m]]$intercept * sig[m, ],
                             sweep(model.list[[m]]$B, 2, sig[m, ], "*"))

  }
  return(list(models = model.list,
              posterior = post.old,
              prior = prior,
              likelihood = likelihood))
}


ComputePiMix <- function(model.list, tau, lambda, gam, prior) {
  # Computes a new estimate for the prior probabilities
  # using a line search
  d <- 0.1
  M <- length(model.list)
  tau.means <- colMeans(tau)
  if (gam == 0) {
    # closed form optimal solution
    return(tau.means)
  } else {
    N <- nrow(tau)
    ComputePiMixObj <- function(new.prior) {
      # Compute objective value for the
      # optimization of prior probabilities
      obj <- -1/N * sum(sweep(tau, 2, log(new.prior), "*"))
      for (m in 1:M) {
        # compute penalty for new.prior
        obj <- obj + (new.prior[m]/prior[m])^gam * model.list[[m]]$pen
      }
      return(obj)
    }
    old.obj <- ComputePiMixObj(prior)
    obj <- ComputePiMixObj(tau.means)
    new.prior <- tau.means
    l <- 0
    # Find largest stepsize d^l in (0, 1] such that
    # moving in the direction of (tau.means - prior)
    # does not increase the objective
    while (obj > old.obj) {
      l <- l + 1
      new.prior <- prior + d^l * (tau.means - prior)
      obj <- ComputePiMixObj(new.prior)
    }
    return(new.prior)
  }
}
