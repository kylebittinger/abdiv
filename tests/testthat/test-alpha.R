context("alpha diversity")

test_that("Alpha diversity values are consistent with QIIME2 tests", {
  x_qiime2 <- c(0, 1, 1, 4, 2, 5, 2, 4, 1, 2)
  expect_equal(berger_parker_d(x_qiime2), 5 / 22)
  expect_equal(brillouin_d(c(1, 2, 0, 0, 3, 1)), 0.8628935302)
  expect_equal(dominance(c(1, 0, 2, 5, 2)), 0.340, tol=0.001)
  # scikit-bio enspie implemented as invsimpson
  expect_equal(invsimpson(c(1, 1, 1, 1, 1, 1)), 6)
  expect_equal(invsimpson(c(13, 13, 13, 13)), 4)
  x_enspie1 <- c(1, 41, 0, 0, 12, 13)
  expect_equal(invsimpson(x_enspie1), 2.250, tol=0.001)
  expect_equal(invsimpson(x_enspie1), 1 / dominance(x_enspie1))
  x_enspie2 <- c(1, 0, 2, 5, 2)
  expect_equal(invsimpson(x_enspie2), 1 / dominance(x_enspie2))
  #expect_equal(fisher_alpha(c(4, 3, 4, 0, 1, 0, 2)), 2.782, tol=0.001)
  #expect_equal(fisher_alpha(c(1, 6, 1, 0, 1, 0, 5)), 2.782, tol=0.001)
  #expect_equal(fisher_alpha(c(61, 0, 0, 1)), 0.395, tol=0.001)
  #expect_equal(fisher_alpha(c(999, 0, 10)), 0.240, tol=0.001)
  x_heip <- c(1, 2, 3, 1)
  expect_equal(heip_e(x_heip), (exp(shannon(x_heip)) - 1) / 3)
  expect_equal(heip_e(c(500, 300, 200)), 0.90, tol=0.01)
  expect_equal(heip_e(c(500, 299, 200, 1)), 0.61, tol=0.01)
  x_kempton <- c(
    2, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6, 7, 7, 9, 9, 11, 14,
    15, 15, 20, 29, 33, 34, 36, 37, 53, 57, 138, 146, 170)
  expect_equal(kempton_taylor_q(x_kempton), 14 / log(34 / 4))
  expect_equal(kempton_taylor_q(sample(x_kempton)), 14 / log(34 / 4))
  expect_equal(margalef(x_qiime2), 8 / log(22))
  expect_equal(mcintosh_d(c(1, 2, 3)), 0.636, tol=0.001)
  expect_equal(mcintosh_e(c(1, 2, 3, 1)), sqrt(15 / 19))
  expect_equal(menhinick(x_qiime2), 9 / sqrt(22))
  expect_equal(richness(c(4, 3, 4, 0, 1, 0, 2)), 5)
  expect_equal(richness(x_qiime2), 9)
  x_pielou <- c(1, 2, 3, 1)
  expect_equal(pielou_e(x_pielou), shannon(x_pielou) / log(4))
  expect_equal(pielou_e(x_qiime2), 0.925, tol=0.001)
  expect_equal(pielou_e(c(1, 1)), 1)
  expect_equal(pielou_e(c(1, 1, 196, 1, 1)), 0.078, tol=0.001)
  expect_equal(pielou_e(c(0, 0, 200, 0, 0)), NaN)
  expect_equal(shannon(5, base = 2), 0)
  expect_equal(shannon(c(5, 5), base = 2), 1)
  expect_equal(shannon(c(1, 1, 1, 1, 0), base = 2), 2)
  expect_equal(simpson(c(1, 0, 2, 5, 2)), 0.66)
  expect_equal(simpson(5), 0)
  expect_equal(simpson_e(c(1, 1, 1, 1, 1, 1, 1)), 1)
  expect_equal(simpson_e(c(500, 400, 600, 500)), 0.980, tol=0.001)
  expect_equal(strong(c(1, 2, 3, 1)), 0.214, tol=0.001)
})

test_that("Alpha diversity values are consistent with vegan", {
  bci19 <- c(
    0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 13, 0, 0, 2, 1, 13, 1, 0, 2,
    0, 0, 0, 2, 4, 0, 1, 0, 4, 0, 0, 1, 2, 0, 8, 3, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 5, 8, 4, 3, 1, 1, 0, 0, 1, 3, 1, 0, 1,
    1, 1, 1, 0, 0, 2, 1, 1, 6, 42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3,
    0, 0, 3, 0, 0, 2, 3, 1, 1, 6, 0, 3, 0, 7, 0, 24, 0, 1, 0, 2,
    0, 1, 0, 0, 0, 0, 0, 0, 3, 1, 1, 7, 1, 4, 1, 3, 4, 2, 0, 3, 1,
    4, 0, 1, 0, 0, 1, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 2, 0, 1, 1,
    24, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 4, 0, 0, 10, 5, 5, 0,
    2, 5, 0, 0, 0, 1, 2, 0, 2, 0, 1, 0, 0, 6, 0, 1, 1, 3, 0, 1, 0,
    0, 0, 1, 6, 4, 0, 0, 0, 9, 2, 0, 0, 0, 0, 9, 0, 0, 0, 0, 3, 1,
    0, 1, 44, 0, 1, 0, 3, 0, 4, 0, 12, 1, 0, 0, 1, 3, 1, 3, 0, 0)
  expect_equal(richness(bci19), vegan::specnumber(bci19))
  expect_equal(shannon(bci19), vegan::diversity(bci19, index = "shannon"))
  expect_equal(simpson(bci19), vegan::diversity(bci19, index = "simpson"))
  expect_equal(invsimpson(bci19), vegan::diversity(bci19, index = "invsimpson"))
})

test_that("Alpha diversity values are consistent with Mothur", {
  mothur_cluster_cts <- function (x) {
    rep(seq_along(x), times=x)
  }
  # https://www.mothur.org/wiki/Sobs
  amazonian10 <- mothur_cluster_cts(c(34, 13, 3, 2, 0, 0, 3))
  expect_equal(richness(amazonian10), 55)
  # https://www.mothur.org/wiki/Bergerparker
  amazonian03 <- mothur_cluster_cts(c(75, 6, 1, 2))
  expect_equal(berger_parker_d(amazonian03), 4 / 98)
  # https://www.mothur.org/wiki/Shannon
  expect_equal(shannon(amazonian03), 4.3532944387)
  # https://www.mothur.org/wiki/Simpson
  # dominance does not match, we don't use the ubiased version
  expect_equal(dominance(amazonian03), 0.0145772595) # 0.004418262 in Mothur
  # no examplee for invsimpson
  # https://www.mothur.org/wiki/Qstat
  # kempton_taylor_q does not match, need to work on this
  expect_equal(kempton_taylor_q(amazonian10), 38.9527661040) # 33.90 in Mothur
})

test_that("Empty count vectors give correct values", {
  x_empty <- c(0, 0, 0, 0)
  expect_equal(berger_parker_d(x_empty), NaN)
  expect_equal(simpson(x_empty), NaN)
  expect_equal(dominance(x_empty), NaN)
  expect_equal(invsimpson(x_empty), NaN)
  expect_equal(simpson_e(x_empty), NaN)
  expect_equal(kempton_taylor_q(x_empty), NaN)
  expect_equal(margalef(x_empty), NaN)
  expect_equal(mcintosh_d(x_empty), NaN)
  expect_equal(mcintosh_e(x_empty), NaN)
  expect_equal(menhinick(x_empty), NaN)
  expect_equal(richness(x_empty), 0) # Only richness gives an answer here
  expect_equal(shannon(x_empty), NaN)
  expect_equal(brillouin_d(x_empty), NaN)
  expect_equal(heip_e(x_empty), NaN)
  expect_equal(pielou_e(x_empty), NaN)
  expect_equal(strong(x_empty), NaN)
})

test_that("Count vector with one observation gives correct values", {
  x_single <- c(0, 0, 0, 1)
  expect_equal(berger_parker_d(x_single), 1) # max = 1, sum = 1, max / sum = 1
  expect_equal(simpson(x_single), 0) # 1 - D = 1 - 1 = 0
  expect_equal(dominance(x_single), 1) # D = sum of squared proportions = 1
  expect_equal(invsimpson(x_single), 1)# 1 / D = 1 / 1 = 1
  expect_equal(simpson_e(x_single), 1) # 1 / (D * S) = 1 / (1 * 1) = 1
  expect_equal(kempton_taylor_q(x_single), 0) # Depends on length of vector
  expect_equal(margalef(x_single), NaN) # See comments in function
  expect_equal(mcintosh_d(x_single), NaN) # n = 1, u = 1, (1 - 1)/(1 - 1)
  expect_equal(mcintosh_e(x_single), 1) # 1 / (1 + 1 - 1) = 1
  expect_equal(menhinick(x_single), 1) # 1 / sqrt(1) = 1
  expect_equal(richness(x_single), 1)
  expect_equal(shannon(x_single), 0) # 1 * log(1) = 0
  expect_equal(brillouin_d(x_single), 0) # 1 * log(1 / 1) = 1 * 0 = 0
  expect_equal(heip_e(x_single), NaN) # exp(H) = exp(0) = 1, (1 - 1) / (1 - 1)
  expect_equal(pielou_e(x_single), NaN) # 0 / log(1) = 0 / 0 = NaN
  expect_equal(strong(x_single), 0) # 1 / 1 - 1 / 1 = 0
})

test_that("Alpha diversity values are correct for simple example", {
  x_simple <- c(0, 1, 1, 4, 3, 2)
  expect_equal(berger_parker_d(x_simple), 4 / 11)
  expect_equal(brillouin_d(x_simple), 1.076, tol=0.001)
  expect_equal(dominance(x_simple), 0.256, tol=0.001)
  expect_equal(heip_e(x_simple), 0.835, tol=0.001)
  expect_equal(invsimpson(x_simple), 3.903, tol=0.001)
  expect_equal(kempton_taylor_q(x_simple), 1.820, tol=0.001)
  expect_equal(margalef(x_simple), 1.668, tol=0.001)
  expect_equal(mcintosh_d(x_simple), 0.707, tol=0.01)
  expect_equal(mcintosh_e(x_simple), 0.765, tol=0.001)
  expect_equal(menhinick(x_simple), 1.508, tol=0.001)
  expect_equal(pielou_e(x_simple), 0.912, tol=0.001)
  expect_equal(richness(x_simple), 5)
  expect_equal(shannon(x_simple), 1.468, tol=0.001)
  expect_equal(simpson(x_simple), 0.744, tol=0.001)
  expect_equal(simpson_e(x_simple), 0.781, tol=0.001)
  expect_equal(strong(x_simple), 0.236, tol=0.001)
})

