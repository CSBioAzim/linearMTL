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
