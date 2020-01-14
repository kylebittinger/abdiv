context("beta diversity components")

test_that("Components match examples from betapart", {
  # Examples run with betapart 1.5.1
  AL <- c(
    1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
    1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0,
    1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0,
    0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0,
    1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0,
    1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0,
    1, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0,
    0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1,
    1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0)
  AT <- c(
    1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0,
    0, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0,
    1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 1,
    0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1,
    1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
    1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    1, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0,
    0, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0,
    1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0,
    1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0,
    1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1,
    1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1,
    1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0)
  expect_equal(sorenson_nestedness(AL, AT), 0.1080251572)
  expect_equal(sorenson_turnover(AL, AT), 0.2893081761)
  expect_equal(sorenson(AL, AT), 0.3973333333)

  expect_equal(jaccard_nestedness(AL, AT), 0.1199218023)
  expect_equal(jaccard_turnover(AL, AT), 0.4487804878)
  expect_equal(jaccard(AL, AT), 0.5687022901)

  BCI1 <- c(
    0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 25, 0, 0, 0, 1,
    13, 2, 0, 0, 6, 0, 0, 4, 5, 0, 0, 0, 1, 0, 0,
    2, 2, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0,
    0, 0, 0, 2, 12, 8, 0, 2, 0, 0, 0, 2, 0, 0, 1,
    1, 2, 0, 0, 0, 0, 0, 0, 0, 3, 14, 0, 0, 0, 1,
    0, 0, 0, 0, 1, 0, 4, 0, 3, 1, 0, 2, 6, 0, 1, 10,
    0, 5, 0, 4, 0, 21, 0, 0, 0, 2, 0, 0, 0, 0, 0,
    0, 3, 0, 2, 0, 0, 6, 1, 1, 0, 0, 0, 0, 0, 0, 7,
    1, 0, 4, 0, 1, 2, 0, 2, 0, 0, 1, 1, 0, 0, 0, 1,
    1, 0, 0, 0, 1, 22, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    3, 2, 0, 24, 5, 0, 5, 0, 13, 5, 2, 11, 0, 0, 0,
    1, 11, 0, 3, 0, 0, 0, 0, 14, 3, 0, 1, 15, 0, 1,
    0, 1, 2, 1, 3, 1, 0, 1, 1, 9, 6, 0, 1, 1, 0, 5,
    0, 1, 0, 0, 3, 0, 0, 0, 18, 0, 0, 2, 1, 0, 1,
    0, 17, 4, 0, 0, 1, 3, 0, 2, 0, 0)
  BCI2 <- c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 26, 0, 0, 0, 0,
    12, 0, 0, 2, 0, 1, 0, 5, 2, 0, 2, 0, 1, 0, 0,
    1, 0, 0, 5, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 3, 14, 6, 0, 2, 0, 0, 0, 2, 3, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 0, 0, 2, 36, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 6, 16, 0, 5, 5,
    0, 9, 0, 5, 0, 14, 0, 2, 0, 4, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 2, 0, 10, 0, 0, 1, 1, 0, 0, 0, 0,
    7, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 3, 1, 1, 16, 3, 0, 7, 0, 12, 4, 0, 8, 0, 0,
    0, 0, 12, 0, 2, 0, 0, 0, 0, 6, 2, 0, 0, 22, 0,
    1, 0, 1, 0, 2, 3, 4, 1, 0, 2, 5, 1, 0, 0, 0, 0,
    7, 0, 1, 0, 1, 1, 0, 0, 1, 27, 0, 0, 0, 1, 1,
    5, 0, 12, 3, 0, 0, 0, 4, 0, 2, 0, 0)
  expect_equal(bray_curtis_balanced(BCI1, BCI2), 0.2597701149)
  expect_equal(bray_curtis_gradient(BCI1, BCI2), 0.0108980617)
  expect_equal(bray_curtis(BCI1, BCI2), 0.2706681767)

  expect_equal(ruzicka_balanced(BCI1, BCI2), 0.4124087591)
  expect_equal(ruzicka_gradient(BCI1, BCI2), 0.0136161963)
  expect_equal(ruzicka(BCI1, BCI2), 0.4260249554)
})

test_that("Sorenson components match 2010 paper", {
  A1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  A2 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
  A3 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  # A1, A2, and A3 are completely nested
  expect_equal(sorenson_turnover(A1, A2), 0)
  expect_equal(sorenson_nestedness(A1, A2), 0.5)
  expect_equal(sorenson_turnover(A2, A3), 0)
  expect_equal(sorenson_nestedness(A2, A3), 1 / 3)

  B1 <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  B2 <- c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0)
  B3 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1)
  # B1, B2, and B3 have no nesting
  expect_equal(sorenson_turnover(B1, B2), 0.5)
  expect_equal(sorenson_nestedness(B1, B2), 0)
  expect_equal(sorenson_turnover(B2, B3), 0.5)
  expect_equal(sorenson_nestedness(B2, B3), 0)

  # Panels C and D are not relevant for pairwise diversity
})

test_that("Sorenson and Jaccard components match 2012 paper", {
  A1 <- c(0, 1, 1)
  A2 <- c(1, 1, 0)
  expect_equal(sorenson_turnover(A1, A2), 1 / 2)
  expect_equal(sorenson_nestedness(A1, A2), 0)
  expect_equal(jaccard_turnover(A1, A2), 2 / 3)
  expect_equal(jaccard_nestedness(A1, A2), 0)

  B1 <- c(0, 1, 1, 1)
  B2 <- c(1, 1, 0, 0)
  expect_equal(sorenson_turnover(B1, B2), 1 / 2)
  expect_equal(sorenson_nestedness(B1, B2), 3 / 5 - 1 / 2)
  expect_equal(jaccard_turnover(B1, B2), 2 / 3)
  expect_equal(jaccard_nestedness(B1, B2), 3 / 4 - 2 / 3)

  C1 <- c(0, 1, 1, 1, 1)
  C2 <- c(1, 1, 0, 0, 0)
  expect_equal(sorenson_turnover(C1, C2), 1 / 2)
  expect_equal(sorenson_nestedness(C1, C2), 4 / 6 - 1 / 2)
  expect_equal(jaccard_turnover(C1, C2), 2 / 3)
  expect_equal(jaccard_nestedness(C1, C2), 4 / 5 - 2 / 3)

  D1 <- c(0, 1, 1, 1, 1, 1)
  D2 <- c(1, 1, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(D1, D2), 1 / 2)
  expect_equal(sorenson_nestedness(D1, D2), 5 / 7 - 1 / 2)
  expect_equal(jaccard_turnover(D1, D2), 2 / 3)
  expect_equal(jaccard_nestedness(D1, D2), 5 / 6 - 2 / 3)

  E1 <- c(0, 1, 1, 1, 1, 1, 1)
  E2 <- c(1, 1, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(E1, E2), 1 / 2)
  expect_equal(sorenson_nestedness(E1, E2), 6 / 8 - 1 / 2)
  expect_equal(jaccard_turnover(E1, E2), 2 / 3)
  expect_equal(jaccard_nestedness(E1, E2), 6 / 7 - 2 / 3)

  F1 <- c(0, 1, 1, 1, 1, 1, 1, 1)
  F2 <- c(1, 1, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(F1, F2), 1 / 2)
  expect_equal(sorenson_nestedness(F1, F2), 7 / 9 - 1 / 2)
  expect_equal(jaccard_turnover(F1, F2), 2 / 3)
  expect_equal(jaccard_nestedness(F1, F2), 7 / 8 - 2 / 3)

  G1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  G2 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(G1, G2), 2 / 6)
  expect_equal(sorenson_nestedness(G1, G2), 8 / 12 - 2 / 6)
  expect_equal(jaccard_turnover(G1, G2), 2 / 4)
  expect_equal(jaccard_nestedness(G1, G2), 8 / 10 - 2 / 4)

  H1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  H2 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(H1, H2), 2 / 8)
  expect_equal(sorenson_nestedness(H1, H2), 7 / 13 - 2 / 8)
  expect_equal(jaccard_turnover(H1, H2), 2 / 5)
  expect_equal(jaccard_nestedness(H1, H2), 7 / 10 - 2 / 5)

  I1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  I2 <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(I1, I2), 2 / 10)
  expect_equal(sorenson_nestedness(I1, I2), 6 / 14 - 2 / 10)
  expect_equal(jaccard_turnover(I1, I2), 2 / 6)
  expect_equal(jaccard_nestedness(I1, I2), 6 / 10 - 2 / 6)

  J1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  J2 <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0)
  expect_equal(sorenson_turnover(J1, J2), 2 / 12)
  expect_equal(sorenson_nestedness(J1, J2), 5 / 15 - 2 / 12)
  expect_equal(jaccard_turnover(J1, J2), 2 / 7)
  expect_equal(jaccard_nestedness(J1, J2), 5 / 10 - 2 / 7)

  K1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  K2 <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
  expect_equal(sorenson_turnover(K1, K2), 2 / 14)
  expect_equal(sorenson_nestedness(K1, K2), 4 / 16 - 2 / 14)
  expect_equal(jaccard_turnover(K1, K2), 2 / 8)
  expect_equal(jaccard_nestedness(K1, K2), 4 / 10 - 2 / 8)

  L1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  L2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0)
  expect_equal(sorenson_turnover(L1, L2), 2 / 16)
  expect_equal(sorenson_nestedness(L1, L2), 3 / 17 - 2 / 16)
  expect_equal(jaccard_turnover(L1, L2), 2 / 9)
  expect_equal(jaccard_nestedness(L1, L2), 3 / 10 - 2 / 9)
})

test_that("Bray-Curtis components match figure 1 in 2017 paper", {
  # No abundance gradient, only balanced variation
  A1 <- c(30, 5, 15, 20)
  A2 <- c(20, 30, 5, 15)
  expect_equal(bray_curtis_balanced(A1, A2), 50 / 140)
  expect_equal(bray_curtis_gradient(A1, A2), 0)

  # Only abundance gradient
  B1 <- c(30, 30, 30, 30)
  B2 <- c(20, 20, 20, 20)
  expect_equal(bray_curtis_balanced(B1, B2), 0)
  expect_equal(bray_curtis_gradient(B1, B2), 40 / 200)
})

test_that("Phylogenetic components match Leprieur 2012", {
  .a <- c(1, 1, 1, 0, 0, 0, 0, 0)
  .b <- c(0, 1, 1, 1, 0, 0, 0, 0)
  .c <- c(0, 1, 1, 1, 1, 0, 0, 0)
  .d <- c(0, 1, 1, 1, 1, 1, 0, 0)
  .e <- c(0, 1, 1, 1, 1, 1, 1, 0)
  .f <- c(0, 1, 1, 1, 1, 1, 1, 1)

  # The values for total diversity and turnover were obtained by counting
  # branches on the tree. I did not reduce the fractions, in case someone
  # wants to re-check by hand.

  # Values for nestedness were obtained by subtraction, and here I did reduce
  # the fractions for convenience.

  expect_equal(phylosor(.a, .b, leprieur_tree), 2 / 12)
  expect_equal(phylosor_turnover(.a, .b, leprieur_tree), 2 / 12)
  expect_equal(phylosor_nestedness(.a, .b, leprieur_tree), 0)
  expect_equal(unweighted_unifrac(.a, .b, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_turnover(.a, .b, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_nestedness(.a, .b, leprieur_tree), 0)

  expect_equal(phylosor(.a, .c, leprieur_tree), 5 / 15)
  expect_equal(phylosor_turnover(.a, .c, leprieur_tree), 2 / 12)
  expect_equal(phylosor_nestedness(.a, .c, leprieur_tree), 2 / 12)
  expect_equal(unweighted_unifrac(.a, .c, leprieur_tree), 1 / 2)
  expect_equal(unweighted_unifrac_turnover(.a, .c, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_nestedness(.a, .c, leprieur_tree), 3 / 14)

  expect_equal(phylosor(.a, .d, leprieur_tree), 6 / 16)
  expect_equal(phylosor_turnover(.a, .d, leprieur_tree), 2 / 12)
  expect_equal(phylosor_nestedness(.a, .d, leprieur_tree), 5 / 24)
  expect_equal(unweighted_unifrac(.a, .d, leprieur_tree), 6 / 11)
  expect_equal(unweighted_unifrac_turnover(.a, .d, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_nestedness(.a, .d, leprieur_tree), 20 / 77)

  expect_equal(phylosor(.a, .e, leprieur_tree), 8 / 18)
  expect_equal(phylosor_turnover(.a, .e, leprieur_tree), 2 / 12)
  expect_equal(phylosor_nestedness(.a, .e, leprieur_tree), 5 / 18)
  expect_equal(unweighted_unifrac(.a, .e, leprieur_tree), 8 / 13)
  expect_equal(unweighted_unifrac_turnover(.a, .e, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_nestedness(.a, .e, leprieur_tree), 30 / 91)

  expect_equal(phylosor(.a, .f, leprieur_tree), 9 / 19)
  expect_equal(phylosor_turnover(.a, .f, leprieur_tree), 2 / 12)
  expect_equal(phylosor_nestedness(.a, .f, leprieur_tree), 35 / 114)
  expect_equal(unweighted_unifrac(.a, .f, leprieur_tree), 9 / 14)
  expect_equal(unweighted_unifrac_turnover(.a, .f, leprieur_tree), 2 / 7)
  expect_equal(unweighted_unifrac_nestedness(.a, .f, leprieur_tree), 5 / 14)

  expect_equal(phylosor(.b, .c, leprieur_tree), 3 / 15)
  expect_equal(phylosor_turnover(.b, .c, leprieur_tree), 0)
  expect_equal(phylosor_nestedness(.b, .c, leprieur_tree), 3 / 15)
  expect_equal(unweighted_unifrac(.b, .c, leprieur_tree), 3 / 9)
  expect_equal(unweighted_unifrac_turnover(.b, .c, leprieur_tree), 0)
  expect_equal(unweighted_unifrac_nestedness(.b, .c, leprieur_tree), 3 / 9)

  expect_equal(phylosor(.b, .d, leprieur_tree), 4 / 16)
  expect_equal(phylosor_turnover(.b, .d, leprieur_tree), 0)
  expect_equal(phylosor_nestedness(.b, .d, leprieur_tree), 4 / 16)
  expect_equal(unweighted_unifrac(.b, .d, leprieur_tree), 4 / 10)
  expect_equal(unweighted_unifrac_turnover(.b, .d, leprieur_tree), 0)
  expect_equal(unweighted_unifrac_nestedness(.b, .d, leprieur_tree), 4 / 10)

  expect_equal(phylosor(.b, .e, leprieur_tree), 6 / 18)
  expect_equal(phylosor_turnover(.b, .e, leprieur_tree), 0)
  expect_equal(phylosor_nestedness(.b, .e, leprieur_tree), 6 / 18)
  expect_equal(unweighted_unifrac(.b, .e, leprieur_tree), 6 / 12)
  expect_equal(unweighted_unifrac_turnover(.b, .e, leprieur_tree), 0)
  expect_equal(unweighted_unifrac_nestedness(.b, .e, leprieur_tree), 6 / 12)

  expect_equal(phylosor(.b, .f, leprieur_tree), 7 / 19)
  expect_equal(phylosor_turnover(.b, .f, leprieur_tree), 0)
  expect_equal(phylosor_nestedness(.b, .f, leprieur_tree), 7 / 19)
  expect_equal(unweighted_unifrac(.b, .f, leprieur_tree), 7 / 13)
  expect_equal(unweighted_unifrac_turnover(.b, .f, leprieur_tree), 0)
  expect_equal(unweighted_unifrac_nestedness(.b, .f, leprieur_tree), 7 / 13)
})

test_that("Beta diversity components are correct for empty vectors", {
  x <- c(0, 0, 0, 0, 0)
  expect_identical(jaccard_nestedness(x, x), NaN)
  expect_identical(jaccard_turnover(x, x), NaN)
  expect_identical(sorenson_nestedness(x, x), NaN)
  expect_identical(sorenson_turnover(x, x), NaN)
  expect_identical(bray_curtis_balanced(x, x), NaN)
  expect_identical(ruzicka_gradient(x, x), NaN)
  expect_identical(ruzicka_balanced(x, x), NaN)
  expect_identical(unweighted_unifrac_nestedness(x, x, faith_tree), NaN)
  expect_identical(unweighted_unifrac_turnover(x, x, faith_tree), NaN)
  expect_identical(phylosor_nestedness(x, x, faith_tree), NaN)
  expect_identical(phylosor_turnover(x, x, faith_tree), NaN)

  y <- c(1, 1, 1, 1, 1)
  expect_identical(jaccard_nestedness(x, y), NaN)
  expect_identical(jaccard_turnover(x, y), NaN)
  expect_identical(sorenson_nestedness(x, y), NaN)
  expect_identical(sorenson_turnover(x, y), NaN)
  expect_identical(bray_curtis_balanced(x, y), NaN)
  expect_identical(ruzicka_gradient(x, y), NaN)
  expect_identical(ruzicka_balanced(x, y), NaN)
  expect_identical(unweighted_unifrac_nestedness(x, y, faith_tree), NaN)
  expect_identical(unweighted_unifrac_turnover(x, y, faith_tree), NaN)
  expect_identical(phylosor_nestedness(x, y, faith_tree), NaN)
  expect_identical(phylosor_turnover(x, y, faith_tree), NaN)
})
