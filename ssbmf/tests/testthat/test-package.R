test_that("package loads", {
  expect_true(is.character(utils::packageVersion("ssbmf") |> as.character()))
})
