context("unifrac")

t_no_branch_lengths <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 9L, 8L, 6L,
  10L, 10L, 7L, 8L, 9L, 1L, 2L, 3L, 10L, 4L, 5L), .Dim = c(9L,
  2L)), Nnode = 5L, tip.label = c("OTU1", "OTU2", "OTU3", "OTU4",
  "OTU5")), .Names = c("edge", "Nnode", "tip.label"), class = "phylo", order = "cladewise")

skbio_t1 <- structure(list(edge = structure(c(6L, 7L, 8L, 9L, 10L, 10L, 9L,
  7L, 11L, 11L, 7L, 8L, 9L, 10L, 1L, 2L, 3L, 11L, 4L, 5L), .Dim = c(10L,
  2L)), edge.length = c(0, 0, 1, 0.5, 0.5, 0.5, 1, 1.25, 0.75,
  0.75), Nnode = 6L, node.label = c("root", "", "", "", "", ""),
  tip.label = c("OTU1", "OTU2", "OTU3", "OTU4", "OTU5")), .Names = c("edge",
  "edge.length", "Nnode", "node.label", "tip.label"), class = "phylo",
  order = "cladewise")

skbio_minimal <- structure(list(edge = structure(c(3L, 3L, 1L, 2L), .Dim = c(2L,
  2L)), edge.length = c(0.25, 0.25), Nnode = 1L, node.label = "root",
  tip.label = c("OTU1", "OTU2")), class = "phylo", order = "cladewise")

skbio_unobs_root <- structure(list(edge = structure(c(5L, 6L, 6L, 5L, 7L, 7L,
  6L, 1L, 2L, 7L, 3L, 4L), .Dim = c(6L, 2L)), edge.length = c(0.3,
  0.1, 0.2, 1.1, 0.5, 0.7), Nnode = 3L, node.label = c("root",
  "", ""), tip.label = c("OTU1", "OTU2", "OTU3", "OTU4")), class = "phylo",
  order = "cladewise")

unifrac_tree <- structure(list(edge = structure(c(4L, 5L, 5L, 4L, 5L, 1L, 2L,
  3L), .Dim = c(4L, 2L)), edge.length = c(0.1, 0.3, 0.3, 0.4),
  Nnode = 2L, tip.label = c("D", "E", "C")), class = "phylo", order = "cladewise")

b0 <- c(1, 3, 0, 1, 0)
b1 <- c(0, 2, 0, 4, 4)
b2 <- c(0, 0, 6, 2, 1)
b3 <- c(0, 0, 1, 1, 1)
b4 <- c(5, 3, 5, 0, 0)
b5 <- c(0, 0, 0, 3, 5)

sk_cts <- rbind(b0=setNames(b0, paste0("OTU", 1:5)), b1, b2, b3, b4, b5)

test_that("Scikit-bio tests are satisfied", {
  expect_equal(faith_pd(c(0, 0, 0, 0, 0), skbio_t1), 0)

  expect_equal(faith_pd(b0, skbio_t1), 4.5)
  expect_equal(faith_pd(b1, skbio_t1), 4.75)
  expect_equal(faith_pd(b2, skbio_t1), 4.75)
  expect_equal(faith_pd(b3, skbio_t1), 4.75)

  expect_equal(faith_pd(c(1, 0), skbio_minimal), 0.25)

  expect_equal(faith_pd(c(1, 1, 0, 0), skbio_unobs_root), 0.6)
  expect_equal(faith_pd(c(0, 0, 1, 1), skbio_unobs_root), 2.3)

  expect_equal(unweighted_unifrac(b0, b1, skbio_t1), 0.238095238095)
  expect_equal(unweighted_unifrac(b0, b2, skbio_t1), 0.52)
  expect_equal(unweighted_unifrac(b0, b3, skbio_t1), 0.52)
  expect_equal(unweighted_unifrac(b0, b4, skbio_t1), 0.545454545455)
  expect_equal(unweighted_unifrac(b0, b5, skbio_t1), 0.619047619048)
  expect_equal(unweighted_unifrac(b1, b2, skbio_t1), 0.347826086957)
  expect_equal(unweighted_unifrac(b1, b3, skbio_t1), 0.347826086957)
  expect_equal(unweighted_unifrac(b1, b4, skbio_t1), 0.68)
  expect_equal(unweighted_unifrac(b1, b5, skbio_t1), 0.421052631579)
  expect_equal(unweighted_unifrac(b2, b3, skbio_t1), 0)
  expect_equal(unweighted_unifrac(b2, b4, skbio_t1), 0.68)
  expect_equal(unweighted_unifrac(b2, b5, skbio_t1), 0.421052631579)
  expect_equal(unweighted_unifrac(b3, b4, skbio_t1), 0.68)
  expect_equal(unweighted_unifrac(b3, b5, skbio_t1), 0.421052631579)
  expect_equal(unweighted_unifrac(b4, b5, skbio_t1), 1)
  expect_equal(weighted_unifrac(b0, b1, skbio_t1), 2.4)
  expect_equal(weighted_unifrac(b0, b2, skbio_t1), 1.86666666667)
  expect_equal(weighted_unifrac(b0, b3, skbio_t1), 2.53333333333)
  expect_equal(weighted_unifrac(b0, b4, skbio_t1), 1.35384615385)
  expect_equal(weighted_unifrac(b0, b5, skbio_t1), 3.2)
  expect_equal(weighted_unifrac(b1, b2, skbio_t1), 2.26666666667)
  expect_equal(weighted_unifrac(b1, b3, skbio_t1), 0.933333333333)
  expect_equal(weighted_unifrac(b1, b4, skbio_t1), 3.2)
  expect_equal(weighted_unifrac(b1, b5, skbio_t1), 0.8375)
  expect_equal(weighted_unifrac(b2, b3, skbio_t1), 1.33333333333)
  expect_equal(weighted_unifrac(b2, b4, skbio_t1), 1.89743589744)
  expect_equal(weighted_unifrac(b2, b5, skbio_t1), 2.666666666667)
  expect_equal(weighted_unifrac(b3, b4, skbio_t1), 2.666666666667)
  expect_equal(weighted_unifrac(b3, b5, skbio_t1), 1.333333333333)
  expect_equal(weighted_unifrac(b4, b5, skbio_t1), 4.0)
  expect_equal(weighted_normalized_unifrac(b0, b1, skbio_t1), 0.6)
  expect_equal(weighted_normalized_unifrac(b0, b2, skbio_t1), 0.466666666667)
  expect_equal(weighted_normalized_unifrac(b0, b3, skbio_t1), 0.633333333333)
  expect_equal(weighted_normalized_unifrac(b0, b4, skbio_t1), 0.338461538462)
  expect_equal(weighted_normalized_unifrac(b0, b5, skbio_t1), 0.8)
  expect_equal(weighted_normalized_unifrac(b1, b2, skbio_t1), 0.566666666667)
  expect_equal(weighted_normalized_unifrac(b1, b3, skbio_t1), 0.233333333333)
  expect_equal(weighted_normalized_unifrac(b1, b4, skbio_t1), 0.8)
  expect_equal(weighted_normalized_unifrac(b1, b5, skbio_t1), 0.209375)
  expect_equal(weighted_normalized_unifrac(b2, b3, skbio_t1), 0.333333333333)
  expect_equal(weighted_normalized_unifrac(b2, b4, skbio_t1), 0.474358974359)
  expect_equal(weighted_normalized_unifrac(b2, b5, skbio_t1), 0.666666666667)
  expect_equal(weighted_normalized_unifrac(b3, b4, skbio_t1), 0.666666666667)
  expect_equal(weighted_normalized_unifrac(b3, b5, skbio_t1), 0.333333333333)
  expect_equal(weighted_normalized_unifrac(b4, b5, skbio_t1), 1.0)

  expect_error(unweighted_unifrac(b0, b1, t_no_branch_lengths))
})

test_that("Implementation is consistent with GUniFrac package", {
  # Distances generated with GUniFrac version 1.1
  expect_equal(generalized_unifrac(b0, b1, skbio_t1, 0), 0.647619047619048)
  expect_equal(generalized_unifrac(b0, b2, skbio_t1, 0), 0.590861244019139)
  expect_equal(generalized_unifrac(b0, b3, skbio_t1, 0), 0.723574660633484)
  expect_equal(generalized_unifrac(b0, b4, skbio_t1, 0), 0.646626447541779)
  expect_equal(generalized_unifrac(b0, b5, skbio_t1, 0), 0.821256038647343)
  expect_equal(generalized_unifrac(b1, b2, skbio_t1, 0), 0.641976726709297)
  expect_equal(generalized_unifrac(b1, b3, skbio_t1, 0), 0.434782608695652)
  expect_equal(generalized_unifrac(b1, b4, skbio_t1, 0), 0.833135669362084)
  expect_equal(generalized_unifrac(b1, b5, skbio_t1, 0), 0.490045596551042)
  expect_equal(generalized_unifrac(b2, b3, skbio_t1, 0), 0.33859649122807)
  expect_equal(generalized_unifrac(b2, b4, skbio_t1, 0), 0.754926829268293)
  expect_equal(generalized_unifrac(b2, b5, skbio_t1, 0), 0.703251657005612)
  expect_equal(generalized_unifrac(b3, b4, skbio_t1, 0), 0.771428571428571)
  expect_equal(generalized_unifrac(b3, b5, skbio_t1, 0), 0.531027056131377)
  expect_equal(generalized_unifrac(b4, b5, skbio_t1, 0), 1)
  expect_equal(generalized_unifrac(b0, b1, skbio_t1, 0.5), 0.617769602045319)
  expect_equal(generalized_unifrac(b0, b2, skbio_t1, 0.5), 0.525472704687848)
  expect_equal(generalized_unifrac(b0, b3, skbio_t1, 0.5), 0.675072475858197)
  expect_equal(generalized_unifrac(b0, b4, skbio_t1, 0.5), 0.481341923061457)
  expect_equal(generalized_unifrac(b0, b5, skbio_t1, 0.5), 0.80995546052754)
  expect_equal(generalized_unifrac(b1, b2, skbio_t1, 0.5), 0.598521685917017)
  expect_equal(generalized_unifrac(b1, b3, skbio_t1, 0.5), 0.320989734713416)
  expect_equal(generalized_unifrac(b1, b4, skbio_t1, 0.5), 0.817047318097553)
  expect_equal(generalized_unifrac(b1, b5, skbio_t1, 0.5), 0.31789725272092)
  expect_equal(generalized_unifrac(b2, b3, skbio_t1, 0.5), 0.335375036602275)
  expect_equal(generalized_unifrac(b2, b4, skbio_t1, 0.5), 0.61103756582007)
  expect_equal(generalized_unifrac(b2, b5, skbio_t1, 0.5), 0.686639667598045)
  expect_equal(generalized_unifrac(b3, b4, skbio_t1, 0.5), 0.71763688469824)
  expect_equal(generalized_unifrac(b3, b5, skbio_t1, 0.5), 0.420437490178619)
  expect_equal(generalized_unifrac(b4, b5, skbio_t1, 0.5), 1)
  expect_equal(variance_adjusted_unifrac(b0, b1, skbio_t1), 0.61086482600333)
  expect_equal(variance_adjusted_unifrac(b0, b2, skbio_t1), 0.477692393409754)
  expect_equal(variance_adjusted_unifrac(b0, b3, skbio_t1), 0.654812113390644)
  expect_equal(variance_adjusted_unifrac(b0, b4, skbio_t1), 0.316461536812682)
  expect_equal(variance_adjusted_unifrac(b0, b5, skbio_t1), 0.803815008981085)
  expect_equal(variance_adjusted_unifrac(b1, b2, skbio_t1), 0.581412810553451)
  expect_equal(variance_adjusted_unifrac(b1, b3, skbio_t1), 0.281183644615338)
  expect_equal(variance_adjusted_unifrac(b1, b4, skbio_t1), 0.806808747363346)
  expect_equal(variance_adjusted_unifrac(b1, b5, skbio_t1), 0.21897829192221)
  expect_equal(variance_adjusted_unifrac(b2, b3, skbio_t1), 0.335791781867113)
  expect_equal(variance_adjusted_unifrac(b2, b4, skbio_t1), 0.487740375942116)
  expect_equal(variance_adjusted_unifrac(b2, b5, skbio_t1), 0.664433327183042)
  expect_equal(variance_adjusted_unifrac(b3, b4, skbio_t1), 0.703053361713065)
  expect_equal(variance_adjusted_unifrac(b3, b5, skbio_t1), 0.352712330815847)
  expect_equal(variance_adjusted_unifrac(b4, b5, skbio_t1), 1.0)
})

test_that("Methods from 'Expanding the UniFrac Toolbox' paper work", {
  expect_equal(information_unifrac(b0, b1, skbio_t1), 0.1394037899)
})
