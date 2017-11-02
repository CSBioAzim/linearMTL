context("Compute error with MTComputeError")

test_that("MTComputeError can compute the MSE.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0, 5, -10, 4), c(2,-5, 0, 3, 4, 3))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_equal(MTComputeError(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                              Y = Y, X = X,
                              task.specific.features = tsf),
               125/4)
})

test_that("MTComputeError can compute the SE.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0, 5, -10, 4), c(2,-5, 0, 3, 4, 3))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_equal(MTComputeError(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                              Y = Y, X = X,
                              task.specific.features = tsf, normalize = FALSE),
               125)
})

test_that("MTComputeError can compute the MSE for given predictions.", {
  pred <- cbind(c(10, -93), c(-14, -53))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_equal(MTComputeError(Y = Y, pred = pred),
               125/4)
})

test_that("MTComputeError recognizes faulty dimensions.", {
  pred <- as.matrix(c(10, -93))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_error(MTComputeError(Y = Y, pred = pred),
               "Dimensions of pred and Y have to coincide!")
})
