context("Node weight calculation")

test_that("ComputeGroupWeights correctly computes weights.", {
  groups <- rbind(c(1,1,1),
                  c(1,1,0),
                  c(1,0,0),
                  c(0,1,0),
                  c(0,0,1))
  # only first two values are relevant. remaining ones should be ignored
  s <- c(0.2, 0.7, -3, 1e7, Inf)
  expect_equal(ComputeGroupWeights(groups = groups, s = s),
               c(0.8, 0.2*0.3, 0.2*0.7, 0.2*0.7, 0.2))
})
