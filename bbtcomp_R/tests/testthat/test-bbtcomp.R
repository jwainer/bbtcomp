


test_that("some examples", {
  options(mc.cores = parallel::detectCores(logical = FALSE))
  m1 <- bbtcomp(ss, lrope = T)
  expect_true(is_bbt_model(m1))
  tp <- table_pwin(m1)
  expect_equal(dim(tp),c(10,6))
  expect_s3_class(tp, "bbt_pwin_table")
  expect_true(all(c("larger", "smaller") %in% names(tp)))
  printed <- utils::capture.output(print(tp))
  expect_true(any(grepl(">", printed, fixed = TRUE)))
  expect_false("larger" %in% strsplit(printed[1], "\\s+")[[1]])
  expect_s3_class(plot_pwin(m1),"ggplot")
  skip_if_not_installed("gridExtra")
  expect_s3_class(plot_ppc(m1),"bayesplot_grid" )
  expect_equal(dim(table_ppc(m1)),c(4,2))

  expect_true(is_bbt_model(bbtcomp(ss, lrope = T, paired = F)))
  expect_true(is_bbt_model(bbtcomp(ss, lrope_value = 0.2, deal_with_ties = "f")))
  expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d")))
  expect_true(is_bbt_model(bbtcomp(ss, lrope = F, hyper_prior = 1, scale = 2.0)))
  expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d")))
})




