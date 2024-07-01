library(testthat)
library(voivod)

test_that("VoI_log_updatable same as VoI_log when coherent", {
  updatable <- VoI_log_updatable(.1, .02, .01, punotc(.01, .02, .1))
  normal <- VoI_log(.1, .02, .01)
  expect_equal(updatable, normal)
})

test_that("When incoherent, updatable is HIGHER than normal", {
  updatable <- VoI_log_updatable(.1, .02, .01, .08)
  normal <- VoI_log(.1, .02, .01)
  expect_true(updatable > normal)

  updatable2 <- VoI_log_updatable(.3, .2, .5, .02)
  normal2 <- VoI_log(.3, .2, .5)
  expect_true(updatable2 > normal2)
})
