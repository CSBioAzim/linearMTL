context("Matrix precomputation with PrepareMatrices")

test_that("PrepareMatrices returns correct result for shared features only.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  Y <- cbind(c(10, -83), c(-9, -53))

  mats <- PrepareMatrices(X = X, Y = Y,
                          standardize = FALSE)
  expect_equal(mats$XTX, t(X) %*% X)
  expect_equal(mats$XTY, t(X) %*% Y)
})

test_that("PrepareMatrices returns correct result for task specific features only.", {
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  Y <- cbind(c(10, -83), c(-9, -53))

  mats <- PrepareMatrices(Y = Y, task.specific.features = tsf,
                          standardize = FALSE)
  expect_equal(mats$XTX[[1]], t(tsf[[1]]) %*% tsf[[1]])
  expect_equal(mats$XTX[[2]], t(tsf[[2]]) %*% tsf[[2]])
  expect_equal(mats$XTY[, 1, drop = FALSE], t(tsf[[1]]) %*% Y[, 1])
  expect_equal(mats$XTY[, 2, drop = FALSE], t(tsf[[2]]) %*% Y[, 2])
})

test_that("PrepareMatrices returns correct result for mixed features.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  Y <- cbind(c(10, -83), c(-9, -53))

  mats <- PrepareMatrices(Y = Y, X = X, task.specific.features = tsf,
                          standardize = FALSE)
  expect_equal(mats$XTX[[1]], t(cbind(X, tsf[[1]])) %*% cbind(X, tsf[[1]]))
  expect_equal(mats$XTX[[2]], t(cbind(X, tsf[[2]])) %*% cbind(X, tsf[[2]]))
  expect_equal(mats$XTY[, 1, drop = FALSE], t(cbind(X, tsf[[1]])) %*% Y[, 1])
  expect_equal(mats$XTY[, 2, drop = FALSE], t(cbind(X, tsf[[2]])) %*% Y[, 2])
})

test_that("PrepareMatrices returns correct result for mixed features with standardization.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  Y <- cbind(c(10, -83), c(-9, -53))

  mats <- PrepareMatrices(Y = Y, X = X, task.specific.features = tsf)
  expect_equal(mats$XTX[[1]], t(scale(cbind(X, tsf[[1]]))) %*% scale(cbind(X, tsf[[1]])))
  expect_equal(mats$XTX[[2]], t(scale(cbind(X, tsf[[2]]))) %*% scale(cbind(X, tsf[[2]])))
  expect_equal(mats$XTY[, 1, drop = FALSE], t(scale(cbind(X, tsf[[1]]))) %*% scale(Y[, 1]))
  expect_equal(mats$XTY[, 2, drop = FALSE], t(scale(cbind(X, tsf[[2]]))) %*% scale(Y[, 2]))
})

test_that("PrepareMatrices returns correct result for mixed features, single task.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)))
  Y <- as.matrix(c(10, -83))

  mats <- PrepareMatrices(Y = Y, X = X, task.specific.features = tsf,
                          standardize = FALSE)
  expect_equal(mats$XTX[[1]], t(cbind(X, tsf[[1]])) %*% cbind(X, tsf[[1]]))
  expect_equal(mats$XTY, t(cbind(X, tsf[[1]])) %*% Y)
})
