test_that("the fast purity accuracy matches caret's mean balanced accuracy", {
  classes = c("EZB", "MCD", "BN2", "Other")
  truth = c("EZB", "EZB", "MCD", "MCD", "Other", "Other", "Other", "EZB")
  pred  = c("EZB", "Other", "MCD", "EZB", "Other", "MCD", "Other", "EZB")   # no BN2 samples: left out, as caret gives NA
  cm = caret::confusionMatrix(factor(pred, classes), factor(truth, classes))
  expect_equal(GAMBLR.predict:::.purity_accuracy(truth, pred, classes, "balanced_accuracy"),
               mean(cm$byClass[, "Balanced Accuracy"], na.rm = TRUE))
  expect_equal(GAMBLR.predict:::.purity_accuracy(truth, pred, classes, "overall_accuracy"), 5 / 8)
})

test_that("grid snapping picks the nearest value and the higher one on ties", {
  grid = seq(3, 0, -0.05)
  expect_equal(GAMBLR.predict:::.nearest_grid(0.62, grid), 0.6)
  expect_equal(GAMBLR.predict:::.nearest_grid(0.625, grid), 0.65)
})
