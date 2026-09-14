

test_that("no ties ssmeanx all the same",{
  w1 = make_wintable(ssmeanx)
  w2 = make_wintable(ssmeanx, lrope=F)
  w3 = make_wintable(ssmeanx, lrope=T, paired = FALSE)
  w4 = make_wintable(ssmeanx, lrope=T, paired = TRUE)
  w5 = make_wintable(ssmeanx, lrope=F, deal_with_ties = "s")
  w6 = make_wintable(ssmeanx, lrope=F, deal_with_ties = "r")
  w7 = make_wintable(ssmeanx, lrope=F, deal_with_ties = "f")
  w8 = make_wintable(ssmeanx, lrope=F, deal_with_ties = "d")
  expect_equal(w1$table, w2$table)
  expect_equal(w1$table, w3$table)
  expect_equal(w1$table, w4$table)
  expect_equal(w1$table, w5$table)
  expect_equal(w1$table, w6$table)
  expect_equal(w1$table, w7$table)
  expect_equal(w1$table, w8$table)
})

test_that("tabmean and tabsd",{
  expect_equal(make_wintable(amean, sssd[,-1], lrope=T, dbcol=0, paired=F)$table,
               make_wintable(ss, lrope=T, paired = F)$table)
  expect_equal(make_wintable(ssmean[,-1], sssd[,-1], lrope=T, dbcol=0,, paired=F)$table,
               make_wintable(ss, lrope=T, paired = F)$table)
  expect_equal(make_wintable(ssmean[,-1], sssd[,-1], dbcol=0,lrope=T, lrope_value = 0.3, paired=F)$table,
               make_wintable(ss, lrope=T, paired = F, lrope_value = 0.3)$table)
  expect_equal(make_wintable(ssmean[,-1], sssd[,-1], dbcol=0,lrope=T, deal_with_ties = "d", paired=F)$table,
               make_wintable(ss, lrope=T, paired = F, deal_with_ties = "d")$table)
  expect_equal(make_wintable(amean, sssd[,-1], dbcol=0,lrope=T, deal_with_ties = "f", paired=F)$table,
               make_wintable(ss, lrope=T, paired = F, deal_with_ties = "f")$table)
})

test_that("different variations all should work",{
  expect_length(make_wintable(ss, lrope = F, paired = T),6)
  expect_length(make_wintable(ss, lrope = T, paired = F),6)
  expect_length(make_wintable(ss, deal_with_ties = "s"),6)
  expect_length(make_wintable(ss, deal_with_ties = "f"),6)
  expect_length(make_wintable(ss, deal_with_ties = "d"),6)
  expect_length(make_wintable(ss, deal_with_ties = "r", lrope = F),6)
})

