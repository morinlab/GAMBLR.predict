test_that("folds are stratified by class and keep groups together", {
  md = data.frame(sample_id = paste0("s", 1:60),
                  patient_id = paste0("p", rep(1:30, each = 2)),
                  lymphgen = rep(c("EZB", "MCD", "Other"), each = 20))
  f = DLBCLone_cv_folds(md, k = 5)
  expect_equal(f$sample_id, md$sample_id)
  expect_setequal(unique(f$fold), 1:5)
  # each class spread evenly: 20 samples -> 4 per fold
  tab = table(md$lymphgen, f$fold)
  expect_true(all(tab == 4))
  # same seed, same folds
  expect_equal(DLBCLone_cv_folds(md, k = 5), f)
  # grouped: a patient's two samples share a fold
  g = DLBCLone_cv_folds(md, k = 5, group_column = "patient_id")
  expect_true(all(tapply(g$fold, md$patient_id, function(x) length(unique(x))) == 1))
  expect_error(DLBCLone_cv_folds(md, group_column = "nope"), "no nope column")
})

test_that("cv summary scores held-out calls and leaves out unlisted labels", {
  cv = list(predictions = data.frame(sample_id = paste0("s", 1:6), fold = 1,
                                     truth = c("EZB", "EZB", "MCD", "Other", "Other", "EZB/MCD"),
                                     DLBCLone_wo = c("EZB", "Other", "MCD", "EZB", "Other", "EZB"),
                                     DLBCLone_wog = c("EZB", "Other", "MCD", "Other", "Other", "EZB")),
            settings = list(truth_classes = c("EZB", "MCD", "Other"), other_class = "Other"))
  s = DLBCLone_cv_summary(cv)
  expect_equal(s$overall$samples, c(5, 5))
  expect_equal(s$overall$agreement, c(3 / 5, 4 / 5))
  expect_equal(s$overall$other_given_class, c(1, 0))
  ezb = s$per_class[s$per_class$call == "DLBCLone_wo" & s$per_class$class == "EZB", ]
  expect_equal(c(ezb$recall, ezb$precision), c(0.5, 0.5))
  expect_false("Other" %in% s$per_class$class)
})

test_that("cv results from separate folds combine, and duplicates are refused", {
  folds = data.frame(sample_id = c("a", "b"), fold = c(1, 2))
  one = function(f) list(predictions = data.frame(sample_id = folds$sample_id[f], fold = f),
                         params = data.frame(fold = f, k = 5), features_used = setNames(list("x"), f),
                         folds = folds, settings = list())
  cv = DLBCLone_cv_combine(list(one(1), one(2)))
  expect_equal(nrow(cv$predictions), 2)
  expect_error(DLBCLone_cv_combine(list(one(1), one(1))), "more than once")
})
