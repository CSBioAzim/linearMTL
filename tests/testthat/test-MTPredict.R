context("Prediction with MTPredict")

test_that("MTPredict recognizes missing data.", {
  B <- as.matrix(c(1,-4, 0))
  expect_error(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B)))),
               "No input data supplied.")
})

test_that("MTPredict recognizes mismatch in number of rows of B and number of features.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0), c(2,-5, 0))
  expect_error(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         X = X, task.specific.features = tsf),
               "Incorrect number of dimensions: B must have one row per feature!")
})

test_that("MTPredict recognizes mismatch in number of columns of B and length of tsf list.", {
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0))
  expect_error(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         task.specific.features = tsf),
               "B must have one column per element in task.specific.features!")
})

test_that("MTPredict works with only shared features for one task", {
  X <- rbind(c(5, 4, 6), c(3, 0, -10))
  B <- as.matrix(c(1,-4, 0))
  intercept <- 5
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = intercept),
                         X = X),
               as.matrix(c(-6, 8)))
})

test_that("MTPredict works with only shared features and multiple tasks", {
  X <- rbind(c(5, 4, 6), c(3, 0, -10))
  B <- cbind(c(1,-4, 0), c(2,-5, 0))
  intercept <- c(-10, 10)
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = intercept),
                         X = X),
               cbind(c(-21, -7), c(0, 16)))
})

test_that("MTPredict works with only specific features for one task", {
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)))
  B <- as.matrix(c(1,-4, 0))
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         task.specific.features = tsf),
               as.matrix(c(-11, 3)))
})

test_that("MTPredict works with only specific features for multiple tasks", {
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0), c(2,-5, 0))
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         task.specific.features = tsf),
               cbind(c(-11, 3), c(37, -23)))
})

test_that("MTPredict works with both shared and specific features one task.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)))
  B <- cbind(c(1,-4, 0, 5, -10, 4))
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         X = X, task.specific.features = tsf),
               cbind(c(10, -93)))
})

test_that("MTPredict works with both shared and specific features for multiple tasks.", {
  X <- rbind(c(1, 0, -14), c(12, 20, -5))
  tsf <- list(rbind(c(5, 4, 6), c(3, 0, -10)), rbind(c(1, -7, 3), c(1, 5, 0)))
  B <- cbind(c(1,-4, 0, 5, -10, 4), c(2,-5, 0, 3, 4, 3))
  expect_equal(MTPredict(LMTL.model = list(B = B, intercept = rep(0, ncol(B))),
                         X = X, task.specific.features = tsf),
               cbind(c(10, -93), c(-14, -53)))
})
