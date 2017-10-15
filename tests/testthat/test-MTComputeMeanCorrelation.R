context("Compute error with MTComputeMeanCorrelation")

test_that("MTComputeMeanCorrelation can compute the mean correlation.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0, 5, -10, 4), c(2,-5, 0, 3, 4, 3))
  pred <- cbind(c(10, -93), c(-14, -53))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_equal(MTComputeMeanCorrelation(Y = Y, B = B, X = X,
                              task.specific.features = tsf), 1)
})

test_that("MTComputeMeanCorrelation can compute the mean correlation for pred.", {
  pred <- cbind(c(10, -93), c(-14, -53))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_equal(MTComputeMeanCorrelation(Y = Y, pred = pred), 1)
})

test_that("MTComputeMeanCorrelation recognizes faulty dimensions.", {
  pred <- as.matrix(c(10, -93))
  Y <- cbind(c(10, -83), c(-9, -53))
  expect_error(MTComputeMeanCorrelation(Y = Y, pred = pred),
               "Dimensions of pred and Y have to coincide!")
})
