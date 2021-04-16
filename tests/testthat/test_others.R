context("Test utilities function")

A = matrix(runif(9), 3, 3)
eig = eigen(A)
N = nrow(A)
ev = eig$vectors[, which.max(abs(Re(eig$values)))] %>% Re %>% abs 
pv = power_method(A, 1e-7) %>% abs
test_that("Power method return correct value", {
  expect_equal(ev, pv)
})

