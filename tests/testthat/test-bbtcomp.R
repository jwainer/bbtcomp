


test_that("some examples", {
  (expect_true(is_bbt_model(bbtcomp(ss, lrope = T))))
  (expect_true(is_bbt_model(bbtcomp(ss, lrope = T, paired = F))))

  (expect_true(is_bbt_model(bbtcomp(ss, lrope_value = 0.2, deal_with_ties = "f"))))
  (expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d"))))
  (expect_true(is_bbt_model(bbtcomp(ss, use_log_lik = T))))
  (expect_true(is_bbt_model(bbtcomp(ss, lrope = F, hyper_prior = 1, scale = 2.0))))
  (expect_true(is_bbt_model(bbtcomp(ss, deal_with_ties = "d", use_log_lik = T))))
})




