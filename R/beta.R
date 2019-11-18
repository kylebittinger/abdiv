#' All beta-diversity measures
#'
#' We exclude the phylogenetic measures and functions that return multiple
#' values.
#' @export
beta_diversities <- c(
  "euclidean", "rms_distance", "chord", "clark_coefficient_of_divergence",
  "geodisc_metric", "manhattan", "mean_character_difference",
  "modified_mean_character_difference", "canberra", "chebyshev", "correlation",
  "cosine", "bray_curtis", "hellinger", "kulczynski", "kulczynski_cody",
  "kulczynski_mothur", "kulczynski_scipy", "rogers_tanimoto", "russel_rao",
  "sokal_michener", "sokal_sneath", "yule", "gower", "alt_gower", "minkowski",
  "morisita", "cao", "millar", "morisita_horn", "jaccard", "sorenson",
  "whittaker", "hamming")

#' Euclidean distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in dist() function with method = "euclidean".
#'   \item Equivalent to vegdist() with method = "euclidean".
#'   \item Equivalent to euclidean() function in scipy.spatial.distance.
#'   \item Equivalent to D_1 in Legendre & Legendre.
#'   \item Equivalent to D_18 in Legendre & Legendre after transformation to
#'     relative abundance.
#' }
#' @export
euclidean <- function (x, y) {
  sqrt(sum((y - x) ^ 2))
}

#' Kullback-Leibler divergence
#'
#' @param x,y Numeric vectors representing probabilities
#'
#' @details
#' Kullback-Leibler divergence is a non-symmetric measure of difference between
#' two porbability vectors. In general, KL(x, y) is not equal to KL(y, x).
#'
#' Because this measure is defined for probabilities, the vectors x and y are
#' normalized in the function so they sum to 1.
#'
#' The Kullback-Leibler divergence is not defined when y_i == 0 but x_i > 0. In
#' this case, the function returns NaN.
#' @export
kullback_leibler_divergence <- function (x, y) {
  x <- x / sum(x)
  y <- y / sum(y)
  y0_but_not_x0 <- (y == 0) & (!(x == 0))
  if (any(y0_but_not_x0)) {
    return(NaN)
  }
  terms <- x * log(x / y)
  sum(ifelse(x > 0, terms, 0))
}

#' Root-mean-square distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_XXXX in Legendre & Legendre.
#' }
#' @export
rms_distance <- function (x, y) {
  sqrt(mean((y - x) ^ 2))
}

#' Chord distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_3 in Legendre & Legendre.
#' }
#' @export
chord <- function (x, y) {
  x <- x / sqrt(sum(x ^ 2))
  y <- y / sqrt(sum(y ^ 2))
  euclidean(x, y)
}

#' Clark's coefficient of divergence
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_11 in Legendre & Legendre.
#' }
#' @export
clark_coefficient_of_divergence <- function (x, y) {
  keep <- (x > 0) | (y > 0)
  x <- x[keep]
  y <- y[keep]
  sqrt(sum(((x - y) / (x + y)) ^ 2) / length(x))
}

#' Geodisc metric
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_4 in Legendre & Legendre.
#' }
#' @export
geodisc_metric <- function (x, y) {
  acos(1 - chord(x, y) / 2)
}

#' Manhattan or city block distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in dist() function with method = "manhattan".
#'   \item Equivalent to vegdist() with method = "manhattan".
#'   \item Equivalent to cityblock() function in scipy.spatial.distance.
#'   \item Equivalent to D_7 in Legendre & Legendre.
#'   \item Whittaker's index of assiciation (D_9 in Legendre & Legendre) is the
#' Manhattan distance computed after transforming to proportions and dividing
#' by 2.
#' }
#' @export
manhattan <- function (x, y) {
  sum(abs(y - x))
}

#' Mean character difference
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_8 in Legendre & Legendre.
#'   \item For binary data, equivalent to 1 - S_1 in Legendre & Legendre.
#' }
#' @export
mean_character_difference <- function (x, y) {
  manhattan(x, y) / length(x)
}

#' Modified mean character difference
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_19 in Legendre & Legendre.
#'   \item For binary data, equivalent to Jaccard distance.
#' }
#' @export
modified_mean_character_difference <- function (x, y) {
  pp <- sum((x > 0) | (y > 0))
  manhattan(x, y) / pp
}

#' Canberra distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in dist() function with method = "canberra".
#'   \item Equivalent to vegdist() with method = "canberra", multiplied by the number of nonzero entries.
#'   \item Equivalent to canberra() function in scipy.spatial.distance.
#'   \item Equivalent to D_10 in Legendre & Legendre.
#' }
#' @export
canberra <- function (x, y) {
  numerator <- abs(x - y)
  denominator <- abs(x + y)
  keep <- denominator != 0
  sum(numerator[keep] / denominator[keep])
}

#' Chebyshev distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to chebyshev() function in scipy.spatial.distance.
#' }
#' @export
chebyshev <- function (x, y) {
  max(abs(x - y))
}

#' Correlation distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to correlation() function in scipy.spatial.distance.
#' }
#' @export
correlation <- function (x, y) {
  x_centered <- x - mean(x)
  y_centered <- y - mean(y)
  cosine(x_centered, y_centered)
}

#' Cosine distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to cosine() function in scipy.spatial.distance.
#' }
#' @export
cosine <- function (x, y) {
  1 - mean(x * y) / sqrt(mean(x * x) * mean(y * y))
}

#' Bray-Curtis distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "bray".
#'   \item Equivalent to D_14 = 1 - s_17 in Legendre & Legendre.
#' }
#' @export
bray_curtis <- function (x, y) {
  sum(abs(x - y)) / sum(x + y)
}

#' Hellinger distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to D_17 in Legendre & Legendre.
#' }
#' @export
hellinger <- function (x, y) {
  x <- x / sum(x)
  y <- y / sum(y)
  chord(sqrt(x), sqrt(y))
}

#' Kulczynski distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "kulczynski".
#' }
#' @export
kulczynski <- function (x, y) {
  min_sum <- sum(pmin(x, y))
  x_sum <- sum(x)
  y_sum <- sum(y)
  1 - 0.5 * ((min_sum / x_sum) + (min_sum / y_sum))
}

#' Kulczynski-Cody distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to 1 - S_13 in Legendre & Legendre.
#' }
#' @export
kulczynski_cody <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  1 - 0.5 * (a / sum(x) + a / sum(y))
}

#' Kulczynski distance (Mothur implementation)
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to 1 - S_12 in Legendre & Legendre.
#' }
#' @export
kulczynski_mothur <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  1 - a / (sum(x) + sum(y) - 2 * a)
}

#' Kulczynski distance (scipy implementation)
#'
#' @param x,y Numeric vectors
#'
#' @export
kulczynski_scipy <- function (x, y) {
  n <- length(x)
  with(scipy_coefficients(x, y), {
    (cTF + cFT - cTT + n) / (cTF + cFT + n)
  })
}

scipy_coefficients <- function (x, y) {
  x <- x > 0
  y <- y > 0
  list(
    cTT = sum(x & y),
    cFT = sum((!x) & y),
    cTF = sum(x & (!y)),
    cFF = sum((!x) & (!y)))
}

#' Rogers-Tanimoto distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to rogerstanimoto() function in scipy.spatial.distance.
#'   \item Equivalent to 1 - S_2 in Legendre & Legendre
#' }
#' @export
rogers_tanimoto <- function (x, y) {
  with(scipy_coefficients(x, y), {
    R <- 2 * (cTF + cFT)
    R / (cTT + cFF + R)
  })
}

#' Russel-Rao distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to russelrao() function in scipy.spatial.distance.
#'   \item Equivalent to 1 - S_11 in Legendre & Legendre
#' }
#' @export
russel_rao <- function (x, y) {
  x <- x > 0
  y <- y > 0
  cTT <- sum(x & y)
  n <- length(x)
  (n - cTT) / n
}

#' Sokal-Michener distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to... XXXXX
#' }
#' @export
sokal_michener <- function (x, y) {
  with(scipy_coefficients(x, y), {
    R <- 2 * (cTF + cFT)
    S <- cFF + cTT
    R / (S + R)
  })
}

#' Sokal-Sneath distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to sokalsneath() function in scipy.spatial.distance.
#'   \item Equivalent to 1 - S_10 in Legendre & Legendre.
#' }
#' @export
sokal_sneath <- function (x, y) {
  with(scipy_coefficients(x, y), {
    R <- 2 * (cTF + cFT)
    R / (cTT + R)
  })
}

#' Yule distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to yule() function in scipy.spatial.distance.
#'   \item Equivalent to 1 - S, where S is the Yule coefficient in Legendre & Legendre.
#' }
#' @export
yule <- function (x, y) {
  with(scipy_coefficients(x, y), {
    2 * cTF * cFT / (cTT * cFF + cTF * cFT)
  })
}

make_range_scale_fcn <- function (x, y) {
  xy_min <- pmin(x, y)
  xy_max <- pmax(x, y)
  xy_range <- xy_max - xy_min
  xy_range_is_zero <- xy_range == 0
  function (z) {
    ifelse(xy_range_is_zero, 0, (z - xy_min) / xy_range)
  }
}

#' Gower distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "gower".
#' }
#' @export
gower <- function (x, y) {
  range_scale <- make_range_scale_fcn(x, y)
  x <- range_scale(x)
  y <- range_scale(y)
  sum(abs(x - y)) / length(x)
}

#' Alternate Gower distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "altGower".
#' }
#' @export
alt_gower <- function (x, y) {
  sum(abs(x - y)) / sum((x > 0) | (y > 0))
}

#' Minkowski distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in dist() function with method = "minkowski".
#'   \item Equivalent to minkowski() function in scipy.spatial.distance.
#'   \item Equivalent to D_6 in Legendre & Legendre.
#' }
#' @export
minkowski <- function (x, y, p = 1) {
  sum(abs(x - y) ^ p) ^ (1 / p)
}

#' Morisita index of dissimilarity
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "morisita".
#' }
#' @export
morisita <- function (x, y) {
  # BUG IN VEGAN?!?!?!
  # Documentation says denominator of lambda should be sum(x) * sum(x - 1)
  lambda_x <- sum(x * (x - 1)) / (sum(x) * (sum(x) - 1))
  lambda_y <- sum(y * (y - 1)) / (sum(y) * (sum(y) - 1))
  xy_term <- sum(x * y) / (sum(x) * sum(y))
  1 - 2 * xy_term / (lambda_x + lambda_y)
}

# TODO: Mountford
# TODO: Raup
# TODO: Chao-Jaccard, vegan uses some correction from the paper
# DONT DO: mahalanobis, in R stats module
# Mahalanobis equivalent to D_5 in Legendre & Legendre.


#' Cao index of dissimilarity
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "cao".
#' }
#' @export
cao <- function (x, y) {
  keep <- (x > 0) | (y > 0)
  x <- x[keep]
  y <- y[keep]
  # using truncation at 0.1, just like vegan
  x_trunc <- ifelse(x > 0.1, x, 0.1)
  y_trunc <- ifelse(y > 0.1, y, 0.1)
  s <- length(x)
  n <- x_trunc + y_trunc
  t1 <- log(n) - log(2)
  t2 <- x_trunc * log(y_trunc)
  t3 <- y_trunc * log(x_trunc)
  sum(t1 - (t2 + t3) / n) / s
}

#' Binomial index of dissimilarity
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "binomial".
#' }
#' @export
millar <- function (x, y) {
  keep <- (x > 0) | (y > 0)
  x <- x[keep]
  y <- y[keep]
  n <- x + y
  t1 <- ifelse(x > 0, x * log(x / n), 0)
  t2 <- ifelse(y > 0, y * log(y / n), 0)
  t3 <- n * log(2)
  sum((t1 + t2 + t3) / n)
}

#' Morisita-Horn index of dissimilarity
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "horn".
#' }
#' @export
morisita_horn <- function (x, y) {
  lambda_x <- sum(x ^ 2) / (sum(x) ^ 2)
  lambda_y <- sum(y ^ 2) / (sum(y) ^ 2)
  xy_term <- sum(x * y) / (sum(x) * sum(y))
  1 - 2 * xy_term / (lambda_x + lambda_y)
}

#' Jaccard index of dissimilarity
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in dist() function with method = "binary".
#'   \item Equivalent to vegdist() with method = "jaccard" and binary = TRUE.
#'   \item Equivalent to jaccard() function in scipy.spatial.distance.
#'   \item Equivalent to 1 - S7 in Legendre & Legendre.
#' }
#' @export
jaccard <- function (x, y) {
  x <- x > 0
  y <- y > 0
  sum(xor(x, y)) / sum(x | y)
}

koleff_abc <- function (x, y) {
  x <- x > 0
  y <- y > 0
  list(c = sum(x & (!y)), b = sum(y & (!x)), a = sum(x & y))
}

#' Sorenson or Dice index of dissimilarity
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to dice() function in scipy.spatial.distance.
#'   \item Equivalent to D_13 = 1 - S_8 in Legendre & Legendre.
#' }
#' @export
sorenson <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  bc / (2 * a + bc)
}

#' Whittaker index of dissimilarity
#' @export
whittaker <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  (a + bc) / (((2 * a) + bc) / 2)
}

#' Hamming distance
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to hamming() function in scipy.spatial.distance
#'   \item For binary data, equivalent to 1 - S_1 in Legendre & Legendre.
#' }
#' @export
hamming <- function (x, y) {
  sum(x != y) / length(x)
}

# Legendre & Legendre notes:
# S_3, S_4, S_5, S_6 not implemented
# Hamann and Pearson's phi not implemented
# S_9 = 3a / (3a + b + c) not implemented
# S_14 as a distance is proportional to sqrt of chord or Hellinger for binary data
# S_26 is not implemented
# Gower's coefficient, S_15 and S_16 not implemented. I think they are equivalent to Hamming distance.
# Skip rest of Q-mode similarity for now
# Pearson coefficient of racial likeness, D_12, not implemented
# Chi-square distances, D_15 and D_16, not implemented (they need the full matrix)
# Need to note how Bray-Curtis distance is related to almost everything else.

