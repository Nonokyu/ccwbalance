test_that("plot functions return ggplot objects", {
  w <- small_weighted()
  feats <- c("X", "L", "L_sq")
  r <- rrb(w, feats, method = c("unweighted", "weighted"))
  expect_s3_class(plot_love_dynamic(r), "ggplot")
  expect_s3_class(plot_balance_heatmap(r[r$method == "weighted", ]), "ggplot")
  expect_error(plot_balance_heatmap(r), "single method")
  expect_s3_class(plot_weight_support(weight_support(w)), "ggplot")
  expect_s3_class(plot_ecdf_pre_post(w, "L", "g3"), "ggplot")
  expect_s3_class(plot_mahalanobis(mahalanobis_rr(w, feats)), "ggplot")
})

test_that("plot_dashboard composes all panels", {
  w <- small_weighted()
  db <- plot_dashboard(w, c("X", "L", "L_sq"))
  plots <- attr(db, "plots")
  expect_named(plots, c("love", "heatmap", "mahalanobis", "support", "ecdf"))
  expect_true(all(vapply(plots, ggplot2::is.ggplot, logical(1))))
})
