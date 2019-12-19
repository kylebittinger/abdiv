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
  expect_equal(sorenson_nestedness_component(AL, AT), 0.1080251572)
  expect_equal(sorenson_turnover_component(AL, AT), 0.2893081761)
  expect_equal(sorenson(AL, AT), 0.3973333333)

  expect_equal(jaccard_nestedness_component(AL, AT), 0.1199218023)
  expect_equal(jaccard_turnover_component(AL, AT), 0.4487804878)
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
  expect_equal(bray_curtis_balanced_component(BCI1, BCI2), 0.2597701149)
  expect_equal(bray_curtis_gradient_component(BCI1, BCI2), 0.0108980617)
  expect_equal(bray_curtis(BCI1, BCI2), 0.2706681767)

  expect_equal(ruzicka_balanced_component(BCI1, BCI2), 0.4124087591)
  expect_equal(ruzicka_gradient_component(BCI1, BCI2), 0.0136161963)
  expect_equal(ruzicka(BCI1, BCI2), 0.4260249554)
})

test_that("Sorenson components match 2010 paper", {
  A1 <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  A2 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0)
  A3 <- c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  # A1, A2, and A3 are completely nested
  expect_equal(sorenson_turnover_component(A1, A2), 0)
  expect_equal(sorenson_nestedness_component(A1, A2), 0.5)
  expect_equal(sorenson_turnover_component(A2, A3), 0)
  expect_equal(sorenson_nestedness_component(A2, A3), 1 / 3)

  B1 <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  B2 <- c(1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0)
  B3 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1)
  # B1, B2, and B3 have no nesting
  expect_equal(sorenson_turnover_component(B1, B2), 0.5)
  expect_equal(sorenson_nestedness_component(B1, B2), 0)
  expect_equal(sorenson_turnover_component(B2, B3), 0.5)
  expect_equal(sorenson_nestedness_component(B2, B3), 0)

  # Panels C and D are not relevant for pairwise diversity
})

test_that("Sorenson and Jaccard components match 2012 paper", {
  A1 <- c(0, 1, 1)
  A2 <- c(1, 1, 0)
  expect_equal(sorenson_turnover_component(A1, A2), 1 / 2)
  expect_equal(sorenson_nestedness_component(A1, A2), 0)
  expect_equal(jaccard_turnover_component(A1, A2), 2 / 3)
  expect_equal(jaccard_nestedness_component(A1, A2), 0)

  B1 <- c(0, 1, 1, 1)
  B2 <- c(1, 1, 0, 0)
  expect_equal(sorenson_turnover_component(B1, B2), 1 / 2)
  expect_equal(sorenson_nestedness_component(B1, B2), 3 / 5 - 1 / 2)
  expect_equal(jaccard_turnover_component(B1, B2), 2 / 3)
  expect_equal(jaccard_nestedness_component(B1, B2), 3 / 4 - 2 / 3)

  C1 <- c(0, 1, 1, 1, 1)
  C2 <- c(1, 1, 0, 0, 0)
  expect_equal(sorenson_turnover_component(C1, C2), 1 / 2)
  expect_equal(sorenson_nestedness_component(C1, C2), 4 / 6 - 1 / 2)
  expect_equal(jaccard_turnover_component(C1, C2), 2 / 3)
  expect_equal(jaccard_nestedness_component(C1, C2), 4 / 5 - 2 / 3)

  D1 <- c(0, 1, 1, 1, 1, 1)
  D2 <- c(1, 1, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(D1, D2), 1 / 2)
  expect_equal(sorenson_nestedness_component(D1, D2), 5 / 7 - 1 / 2)
  expect_equal(jaccard_turnover_component(D1, D2), 2 / 3)
  expect_equal(jaccard_nestedness_component(D1, D2), 5 / 6 - 2 / 3)

  E1 <- c(0, 1, 1, 1, 1, 1, 1)
  E2 <- c(1, 1, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(E1, E2), 1 / 2)
  expect_equal(sorenson_nestedness_component(E1, E2), 6 / 8 - 1 / 2)
  expect_equal(jaccard_turnover_component(E1, E2), 2 / 3)
  expect_equal(jaccard_nestedness_component(E1, E2), 6 / 7 - 2 / 3)

  F1 <- c(0, 1, 1, 1, 1, 1, 1, 1)
  F2 <- c(1, 1, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(F1, F2), 1 / 2)
  expect_equal(sorenson_nestedness_component(F1, F2), 7 / 9 - 1 / 2)
  expect_equal(jaccard_turnover_component(F1, F2), 2 / 3)
  expect_equal(jaccard_nestedness_component(F1, F2), 7 / 8 - 2 / 3)

  G1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  G2 <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(G1, G2), 2 / 6)
  expect_equal(sorenson_nestedness_component(G1, G2), 8 / 12 - 2 / 6)
  expect_equal(jaccard_turnover_component(G1, G2), 2 / 4)
  expect_equal(jaccard_nestedness_component(G1, G2), 8 / 10 - 2 / 4)

  H1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  H2 <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(H1, H2), 2 / 8)
  expect_equal(sorenson_nestedness_component(H1, H2), 7 / 13 - 2 / 8)
  expect_equal(jaccard_turnover_component(H1, H2), 2 / 5)
  expect_equal(jaccard_nestedness_component(H1, H2), 7 / 10 - 2 / 5)

  I1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  I2 <- c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(I1, I2), 2 / 10)
  expect_equal(sorenson_nestedness_component(I1, I2), 6 / 14 - 2 / 10)
  expect_equal(jaccard_turnover_component(I1, I2), 2 / 6)
  expect_equal(jaccard_nestedness_component(I1, I2), 6 / 10 - 2 / 6)

  J1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  J2 <- c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0)
  expect_equal(sorenson_turnover_component(J1, J2), 2 / 12)
  expect_equal(sorenson_nestedness_component(J1, J2), 5 / 15 - 2 / 12)
  expect_equal(jaccard_turnover_component(J1, J2), 2 / 7)
  expect_equal(jaccard_nestedness_component(J1, J2), 5 / 10 - 2 / 7)

  K1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  K2 <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0)
  expect_equal(sorenson_turnover_component(K1, K2), 2 / 14)
  expect_equal(sorenson_nestedness_component(K1, K2), 4 / 16 - 2 / 14)
  expect_equal(jaccard_turnover_component(K1, K2), 2 / 8)
  expect_equal(jaccard_nestedness_component(K1, K2), 4 / 10 - 2 / 8)

  L1 <- c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  L2 <- c(1, 1, 1, 1, 1, 1, 1, 1, 0, 0)
  expect_equal(sorenson_turnover_component(L1, L2), 2 / 16)
  expect_equal(sorenson_nestedness_component(L1, L2), 3 / 17 - 2 / 16)
  expect_equal(jaccard_turnover_component(L1, L2), 2 / 9)
  expect_equal(jaccard_nestedness_component(L1, L2), 3 / 10 - 2 / 9)
})

test_that("Bray-Curtis components match figure 1 in 2017 paper", {
  A1 <- c(30, 5, 15, 20)
  A2 <- c(20, 30, 5, 15)
  A3 <- c(15, 20, 30, 5)
  A4 <- c(5, 15, 20, 30)
  B1 <- c(30, 30, 30, 30)
  B2 <- c(20, 20, 20, 20)
  B3 <- c(15, 15, 15, 15)
  B4 <- c(5, 5, 5, 5)
  C1 <- c(30, 20, 15, 5)
  C2 <- c(20, 30, 20, 20)
  C3 <- c(15, 15, 30, 15)
  C4 <- c(5, 5, 5, 30)
  obs <- usedist::dist_make(fig1, bray_curtis_balanced_component)
  obs2 <- usedist::dist_make(fig1, bray_curtis)
  expect_equal(bray_curtis_balanced_component(B1, B2), 0)
})
