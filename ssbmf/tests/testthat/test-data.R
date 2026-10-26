mk <- function(p, n, prefix = "g", sp = "s") {
  matrix(rnorm(p * n), p, n, dimnames = list(paste0(prefix, 1:p), paste0(sp, 1:n)))
}

test_that("constant feature is rejected", {
  X <- mk(4, 10); X[2, ] <- 3
  expect_error(ssbmf_data(X, rexp(10), rep(1, 10)), "constant")
})

test_that("missing values are rejected with the feature named", {
  X <- mk(4, 10); X[3, 2] <- NA
  expect_error(ssbmf_data(X, rexp(10), rep(1, 10)), "g3")
})

test_that("block sample names must match", {
  A <- mk(4, 10); B <- mk(6, 10, "c", "t")
  expect_error(ssbmf_data(list(expression = A, methylation = B), rexp(10), rep(1, 10)), "sample")
})

test_that("blocks are aligned by sample name", {
  A <- mk(4, 10); B <- mk(6, 10, "c")
  dd <- ssbmf_data(list(a = A, b = B[, 10:1]), rexp(10), rep(1, 10))
  expect_identical(rownames(dd$Y$b), rownames(dd$Y$a))
  expect_equal(dd$Y$b, t(B))
})

test_that("matrix input becomes one block; orientation is honored", {
  X <- mk(4, 10)
  dd <- ssbmf_data(X, rexp(10), rep(1, 10))
  expect_identical(names(dd$Y), "block1"); expect_identical(dim(dd$Y$block1), c(10L, 4L))
  dd2 <- ssbmf_data(t(X), rexp(10), rep(1, 10), orientation = "samples_by_features")
  expect_identical(dd$Y, dd2$Y)
})

test_that("survival inputs are validated", {
  X <- mk(4, 10)
  expect_error(ssbmf_data(X, rexp(9), rep(1, 10)), "one value per sample")
  expect_error(ssbmf_data(X, -rexp(10), rep(1, 10)), "positive")
  expect_error(ssbmf_data(X, rexp(10), rep(2, 10)), "0")
  expect_error(ssbmf_data(X, rexp(10), rep(0, 10)), "no events")
})
