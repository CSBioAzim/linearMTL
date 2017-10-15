context("Automatic tree and weight construction using HC")

test_that("BuildTreeHC correctly computes tree", {
  Y <- cbind(c(1,0,0), c(2,0,0), c(0,1,1))

  hc.result <- BuildTreeHC(Y, threshold = 1)
  expect_equal(hc.result$groups, rbind(c(1,0,0),
                                       c(0,1,0),
                                       c(0,0,1),
                                       c(1,1,0)))
})

test_that("BuildTreeHC correctly computes weights", {
  Y <- cbind(c(1,0,0), c(2,0,0), c(0,1,1))

  hc.result <- BuildTreeHC(Y, threshold = 1)
  expect_equal(hc.result$weights, c(0,0,1,1))
})
