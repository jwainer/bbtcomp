


test_that("some examples", {
  options(mc.cores = parallel::detectCores(logical = FALSE))
  m1 <- bbtcomp(ss, lrope = T)
  expect_true(is_bbt_model(m1))
  expect_equal(dim(table_pwin(m1)),c(10,5))
  expect_s3_class(plot_pwin(m1),"ggplot")
  expect_s3_class(plot_ppc(m1),"bayesplot_grid" )
  expect_equal(dim(table_ppc(m1)),c(4,2))

  expect_true(is_bbt_model(bbtcomp(ss, lrope = T, paired = F)))
  expect_true(is_bbt_model(bbtcomp(ss, lrope_value = 0.2, deal_with_ties = "f")))
  expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d")))
  expect_true(is_bbt_model(bbtcomp(ss, lrope = F, hyper_prior = 1, scale = 2.0)))
  expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d")))
})




