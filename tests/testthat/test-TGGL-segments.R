context("Segments of tggl optimization")

test_that("Shrink correctly computes shrinkage operator.", {
  group.ranges <- rbind(c(1, 2, 2),
                        c(3, 5, 3),
                        c(6, 6, 1))
  A <- rbind(c(3,0),
             c(4,-3),
             c(4,0),
             c(0,0),
             c(3,3),
             c(14,-5))

  expected.output <- rbind(c(0.6,0),
                           c(0.8,-1),
                           c(0.8,0),
                           c(0,0),
                           c(0.6,1),
                           c(1,-1))
  expect_equal(Shrink(A = A, group.ranges = group.ranges),
               expected.output)
})

test_that("CalculateInnerGroupPenalty correctly computes penalty.", {
  group.ranges <- rbind(c(1, 2, 2),
                        c(3, 5, 3),
                        c(6, 6, 1))
  A <- rbind(c(3,0),
             c(4,-3),
             c(4,0),
             c(0,0),
             c(3,3),
             c(14,-5))

  euc.norm <- function(x) sqrt(sum(x^2))
  expected.output <- euc.norm(c(3,4)) + euc.norm(c(0,-3)) +
    euc.norm(c(4,0,3)) + euc.norm(c(0,0,3)) +
    euc.norm(14) + euc.norm(-5)
  expect_equal(CalculateInnerGroupPenalty(A = A, group.ranges = group.ranges),
               expected.output)
})

test_that("ComputeInverseArtesi computes correct transformation.", {
  # train a lasso model with glmnet and use
  # inverse Artesi transformation to recover the coefficients
  # for variables on the orignal scale

  library(glmnet)
  X <- matrix(runif(100, 0, 1), ncol = 2)
  y <- 1 - 2*X[, 1] + X[, 2]

  # glmnet uses biased sd estimator
  sd.alt <- function(a) {
    return(sqrt(mean((a-mean(a))^2)))
  }

  X.sds <- apply(X, 2, sd.alt)
  X.mus <- apply(X, 2, mean)
  Y.sds <- sd.alt(y)
  Y.mus <- mean(y)

  Xst <- scale(X, center = X.mus, scale = X.sds)
  yst <- scale(y, center = Y.mus, scale = Y.sds)

  glm.orig.scale <- glmnet(X, y, lambda = 0.5 * sd.alt(y))
  glm.std <- glmnet(Xst, yst, lambda = 0.5, standardize = FALSE)

  coef.orig.scale <- coefficients(glm.orig.scale)
  coef.std <- coefficients(glm.std)

  B.no_std <- ComputeInverseArtesi(B = as.matrix(coef.std[-1]),
                                   Y.sds = Y.sds, Y.mus = Y.mus,
                                   X.sds = X.sds, X.mus = X.mus,
                                   tsf.sds = list(), tsf.mus = list())
  expect_equal(as.vector(B.no_std), as.vector(coef.orig.scale))
})
