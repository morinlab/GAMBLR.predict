# DLBCLone evidence weights and gate, on a small synthetic training set:
# class A has a specific feature (a_def, only in A) and a shared one (a_com, in A and Other);
# class B has two weaker features (b1, b2) and shares the BCL6-like feature (shared) with Other.
make_training = function(){
  n = 60
  ids = paste0("s", seq_len(3 * n))
  y = rep(c("A", "B", "Other"), each = n)
  # k of the n samples of a class carry the feature
  k = function(m) c(rep(1, m), rep(0, n - m))
  x = data.frame(
    a_def  = c(k(48), k(0),  k(1)),    # specific to A: 48/49 carriers in A
    a_com  = c(k(42), k(3),  k(24)),   # associated with A, 42/69 = 61% of carriers in A
    b1     = c(k(0),  k(36), k(12)),
    b2     = c(k(0),  k(36), k(12)),
    shared = c(k(3),  k(42), k(30)),   # associated with B, 42/75 = 56% in B
    noise  = rep(c(1, 0, 0, 0, 0), length.out = 3 * n),
    row.names = ids
  ) * 2
  list(x = x, labels = data.frame(sample_id = ids, lymphgen = y))
}

test_that("class evidence keeps associated features and weights specific ones 2", {
  tr = make_training()
  ev = DLBCLone_class_evidence(tr$x, tr$labels, truth_column = "lymphgen")
  expect_setequal(unique(ev$class), c("A", "B"))
  expect_false("noise" %in% ev$feature)
  expect_false("Other" %in% ev$class)
  expect_equal(ev$weight[ev$class == "A" & ev$feature == "a_def"], 2)
  expect_equal(ev$weight[ev$class == "A" & ev$feature == "a_com"], 1)
  # associated with B but common in Other: weight 1, not class-specific
  expect_equal(ev$weight[ev$class == "B" & ev$feature == "shared"], 1)
  # a named label vector and a sample_id column in the features give the same result
  x2 = cbind(sample_id = rownames(tr$x), tr$x)
  ev2 = DLBCLone_class_evidence(x2, setNames(tr$labels$lymphgen, tr$labels$sample_id))
  expect_equal(ev2, ev)
})

test_that("the gate keeps calls with enough weighted evidence and adds columns", {
  ev = data.frame(class = c("A", "A", "B", "B", "B"), feature = c("a_def", "a_com", "b1", "b2", "shared"),
                  weight = c(2, 1, 1, 1, 1))
  x = data.frame(a_def = c(2, 0, 0, 0, 0, 0), a_com = c(0, 2, 0, 0, 0, 2), b1 = c(0, 0, 2, 2, 0, 0),
                 b2 = c(0, 0, 2, 0, 0, 0), shared = c(0, 0, 0, 0, 2, 0),
                 row.names = paste0("t", 1:6))
  pred = data.frame(sample_id = paste0("t", 1:6), DLBCLone_wo = c("A", "A", "B", "B", "B", "Other"))
  out = suppressMessages(DLBCLone_evidence_gate(pred, x, ev))
  # t1: one specific feature (2) kept; t2: one weight-1 feature dropped; t3: two features kept;
  # t4, t5: one weight-1 feature dropped; t6: Other untouched
  expect_equal(out$DLBCLone_wog, c("A", "Other", "B", "Other", "Other", "Other"))
  expect_equal(out$DLBCLone_evidence, c(2, 1, 2, 1, 1, NA))
  expect_equal(out$DLBCLone_wo, pred$DLBCLone_wo)
  # per-class minimum with a default
  out2 = suppressMessages(DLBCLone_evidence_gate(pred, x, ev, min_evidence = c(default = 2, B = 1)))
  expect_equal(out2$DLBCLone_wog, c("A", "Other", "B", "B", "B", "Other"))
  expect_equal(out2$DLBCLone_evidence_min, c(2, 2, 1, 1, 1, NA))
  # calls of classes without weights (e.g. composites) are kept
  pred$DLBCLone_wo[6] = "A/B"
  out3 = suppressMessages(DLBCLone_evidence_gate(pred, x, ev))
  expect_equal(out3$DLBCLone_wog[6], "A/B")
  expect_true(is.na(out3$DLBCLone_evidence[6]))
})

test_that("the gate errors on missing samples and counts missing features as absent", {
  ev = data.frame(class = "A", feature = c("a_def", "gone"), weight = c(2, 1))
  x = data.frame(a_def = c(0, 2), row.names = c("t1", "t2"))
  expect_error(DLBCLone_evidence_gate(data.frame(sample_id = "t3", DLBCLone_wo = "A"), x, ev), "no row in features")
  expect_message(out <- DLBCLone_evidence_gate(data.frame(sample_id = c("t1", "t2"), DLBCLone_wo = "A"), x, ev),
                 "not in features")
  expect_equal(out$DLBCLone_wog, c("Other", "A"))
})
