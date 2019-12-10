#' All beta-diversity measures
#'
#' We exclude the phylogenetic measures and functions that return multiple
#' values.
#' @export
beta_diversities <- c(
  "euclidean", "rms_distance", "chord", "clark_coefficient_of_divergence",
  "geodesic_metric", "manhattan", "mean_character_difference",
  "modified_mean_character_difference", "canberra", "chebyshev", "correlation",
  "cosine", "bray_curtis", "hellinger", "kulczynski", "kulczynski_cody",
  "kulczynski_mothur", "kulczynski_scipy", "rogers_tanimoto", "russel_rao",
  "sokal_michener", "sokal_sneath", "yule", "gower", "alt_gower", "minkowski",
  "morisita", "cao", "millar", "morisita_horn", "jaccard", "sorenson",
  "whittaker", "hamming")

#' Euclidean distance and related measures
#'
#' These distance and diversity measures are mathematically similar to the
#' Euclidean distance between two vectors.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' For vectors \code{x} and \code{y}, the Euclidean distance is defined as
#' \deqn{d(x, y) = \sqrt{\Sigma_i^n (x_i - y_i) ^ 2}.}
#' Relation of \code{euclidean()} to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "euclidean"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "euclidean"}.
#'   \item Equivalent to the \code{euclidean()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{D_1}{D_1} in Legendre & Legendre.
#'   \item Equivalent to the \emph{distance between species profiles},
#'     \eqn{D_{18}}{D_18} in Legendre & Legendre if \code{x} and \code{y} are
#'     transformed to relative abundance.
#' }
#'
#' The \emph{root-mean-square} distance or \emph{average} distance is similar
#' to Euclidean distance. As the name implies, it is computed as the square
#' root of the mean of the squared differences between elements of \code{x}
#' and \code{y}:
#' \deqn{d(x, y) = \sqrt{\frac{1}{n} \Sigma_i^n (x_i - y_i) ^ 2}.}
#' Relation of \code{rms_distance()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_2}{D_2} in Legendre & Legendre.
#' }
#'
#' The \emph{chord} distance is the Euclidean distance after scaling each
#' vector by its root sum of squares, \eqn{\hat{x} = \sqrt{\Sigma_i x_i^2}}.
#' The chord distance between any two vectors ranges from 0 to
#' \eqn{\sqrt{2}}{sqrt(2)}. Relation of \code{chord()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_3}{D_3} in Legendre & Legendre.
#' }
#'
#' The \emph{Hellinger} distance is equal to the chord distance computed after
#' a square-root transformation. Relation of \code{hellinger()} to other
#' definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{17}}{D_17} in Legendre & Legendre.
#' }
#'
#' The \emph{geodesic metric} is a transformed version of the chord distance.
#' \deqn{d(x, y) = \textrm{arccos}(1 - \frac{d_c^2(x, y)}{2}),} where \eqn{d_c}
#' is the chord distance. It gives the length of the arc on a hypersphere
#' between the vectors, if the vectors are normalized to unit length. Relation
#' of \code{geodesic_metric()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_4}{D_4} in Legendre & Legendre.
#' }
#' @export
euclidean <- function (x, y) {
  sqrt(sum((y - x) ^ 2))
}

#' @rdname euclidean
#' @export
rms_distance <- function (x, y) {
  sqrt(mean((y - x) ^ 2))
}

#' @rdname euclidean
#' @export
chord <- function (x, y) {
  x <- x / sqrt(sum(x ^ 2))
  y <- y / sqrt(sum(y ^ 2))
  euclidean(x, y)
}

#' @rdname euclidean
#' @export
hellinger <- function (x, y) {
  x <- x / sum(x)
  y <- y / sum(y)
  chord(sqrt(x), sqrt(y))
}

#' @rdname euclidean
#' @export
geodesic_metric <- function (x, y) {
  acos(1 - chord(x, y) / 2)
}

#' Clark's coefficient of divergence
#'
#' @param x,y Numeric vectors
#'
#' Relation of \code{clark_coefficient_of_divergence()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{11}}{D_11} in Legendre & Legendre.
#' }
#' @export
clark_coefficient_of_divergence <- function (x, y) {
  keep <- (x > 0) | (y > 0)
  x <- x[keep]
  y <- y[keep]
  sqrt(sum(((x - y) / (x + y)) ^ 2) / length(x))
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

#' Manhattan or city block distance
#'
#' The Manhattan distance is the sum of absolute differences between the
#' elements of two vectors.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "manhattan"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "manhattan"}.
#'   \item Equivalent to the \code{cityblock()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{D_7}{D_7} in Legendre & Legendre.
#'   \item Whittaker's index of assiciation (\eqn{D_9}{D_9} in Legendre &
#'     Legendre) is the Manhattan distance computed after transforming to
#'     proportions and dividing by 2.
#' }
#' @export
manhattan <- function (x, y) {
  sum(abs(y - x))
}

#' Mean character difference
#'
#' The mean character difference is the Manhattan distance divided by the
#' vector length.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_8}{D_8} in Legendre & Legendre.
#'   \item For binary data, equivalent to \eqn{1 - S_1}{1 - S_1} in Legendre
#'     & Legendre, where \eqn{S_1}{S_1} is the simple matching coefficient.
#' }
#' @export
mean_character_difference <- function (x, y) {
  manhattan(x, y) / length(x)
}

#' Modified mean character difference
#'
#' The modified mean character difference is the Manhattan distance divided by
#' the number of elements where the vectors are not both equal to zero.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{19}}{D_19} in Legendre & Legendre.
#'   \item For binary data, equivalent to the Jaccard distance.
#' }
#' @export
modified_mean_character_difference <- function (x, y) {
  pp <- sum((x > 0) | (y > 0))
  manhattan(x, y) / pp
}

#' Canberra distance
#'
#' The sum over all elements of the absolute difference in abundance divided by
#' the sum.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "canberra"}.
#'   \item Equivalent to the \code{vegdist()} function with
#'     \code{method = "canberra"}, multiplied by the number of nonzero entries.
#'   \item Equivalent to the \code{canberra()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{D_{10}}{D_10} in Legendre & Legendre.
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
#' The Bray-Curtis distance is the sum of absolute differences between the
#' vectors divided by the sum of both vectors.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "bray"}.
#'   \item Equivalent to \eqn{D_{14} = 1 - S_{17}}{D_14 = 1 - S_17} in
#'     Legendre & Legendre.
#' }
#' @export
bray_curtis <- function (x, y) {
  sum(abs(x - y)) / sum(x + y)
}

#' Kulczynski distance
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "kulczynski"}.
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
#'   \item Equivalent to \eqn{1 - S_{13}}{1 - S_13} in Legendre & Legendre.
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
#'   \item Equivalent to \eqn{1 - S_{12}}{1 - S_12} in Legendre & Legendre.
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
#'   \item Equivalent to the \code{rogerstanimoto()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{1 - S_2}{1 - S_2} in Legendre & Legendre.
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
#'   \item Equivalent to the \code{russelrao()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \code{1 - S_{11}}{1 - S_11} in Legendre & Legendre.
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
#'   \item Equivalent to \code{sokalsneath()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{1 - S_{10}}{1 - S_10} in Legendre & Legendre.
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
#'   \item Equivalent to the \code{yule()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{1 - S}, where \eqn{S} is the Yule coefficient
#'     in Legendre & Legendre.
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
#' The Minkowski metric is a generalized form of Euclidean (p=2) and Manhattan
#' (p=1) distance.
#'
#' @param x,y Numeric vectors
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "minkowski"}.
#'   \item Equivalent to the \code{minkowski()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{D_6}{D_6} in Legendre & Legendre.
#' }
#' The default value of \code{p = 1} makes this distance equal to the Manhattan
#' distance.
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
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "binary"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "jaccard"}
#'     and \code{binary = TRUE}.
#'   \item Equivalent to the \code{jaccard()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{1 - S_7}{1 - S_7} in Legendre & Legendre.
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
#' After transforming to presence/absence, the average fraction of species
#' found in one sample but not the other.
#'
#' @details
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{dice()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to \eqn{D_{13} = 1 - S_8}{D_13 = 1 - S_8} in Legendre &
#'     Legendre.
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
#'   \item Equivalent to the \code{hamming()} function in
#'     \code{scipy.spatial.distance}.
#'   \item For binary data, equivalent to \code{1 - S_1}{1 - S_1} in Legendre
#'     & Legendre.
#' }
#' @export
hamming <- function (x, y) {
  sum(x != y) / length(x)
}

# Euclidean-like:
# euclidean, rms_distance, chord, geodesic_metric, hellinger

# Manhattan-like:
# manhattan, mean_character_difference, modified_mean_character_difference

# Canberra-like: clark_coefficient_of_divergence

# Vegan notes:
# TODO: Mountford
# TODO: Raup
# TODO: Chao-Jaccard, vegan uses some correction from the paper

# Legendre & Legendre notes:
# D_1 implemented as euclidean
# D_2 (average distance) implemented as rms_distance
# D_3 implemented as chord
# D_4 implemented as geodesic_metric
# D_5 (Mahalanobis distance) not implemented, available in R stats module
# D_6 implemented as minkowski
# D_7 implemented as manhattan
# D_8 implemented as mean_character_difference
# D_19 implemented as modified_mean_character_difference
# D_9 (Whittaker's index of association) not implemented, noted under manhattan
# D_10 implemented as canberra
# D_11 implemented as clark_coefficient_of_divergence
# D_12 (Pearson coefficient of racial likeness) not implemented
# D_18 (distance between species profiles) not implemented, noted under euclidean
# D_15 (Chi-square metric) not implemented, needs full matrix
# D_16 (Chi-square distance) not implemented, needs full matrix
# D_17 implemented as hellinger
# D_13 implemented as sorenson
# D_14 implemented as bray_curtis
#   Need to note how Bray-Curtis distance is related to almost everything else.
# S_1 (simple matching coefficient) is the mean_character_difference
# S_2 (coefficient of Rogers & Tanimoto) implemented as rogers_tanimoto
# S_3, S_4, S_5, S_6 (Sokal and Sneath) not implemented
# Hamann coefficient not implemented
# Yule coefficient implemented as yule
# Pearson's phi not implemented
# S_7 implemented as jaccard
# S_8 implemented as sorenson
# S_9 = 3a / (3a + b + c) not implemented
# S_10 implemented as sokal_sneath
# S_11 implemented as russel_rao
# S_12 implemented as kulczynski_mothur
# S_13 implemented as kulczynski_cody
#   The Kulczynski distances are a mess!
# S_14 (name unclear) not implemented. It should be noted that S_14 as a
#   distance is proportional to sqrt of chord or Hellinger for binary data.
# S_26 (name unclear) not implemented
# S_15 (Gower coefficient) not implemented. Equivalent to Hamming distance?
# S_16 (Estabrook and Rogers) not implemented
# S_17 (Steinhaus coefficient) implemented as bray_curtis
# S_18 (yet another Kulcynski coefficient) is equivalent to kulczynski, I think
# S_19 (modified Gower) not implemented
# S_20 (Legendre & Chodorowski) not implemented, also similar to Gower
# S_21 (Chi-square similarity) not implemented, see D_15
# S_22 (Goodall 1) not implemented, needs full matrix
# S_23 (Goodall 2) not implemented, needs full matrix
# Raup-Crick not implemented, available as raupcrick() in vegan
# S_27 (Raup-Crick p-value) not implemented
# Skip rest of Q-mode similarity for now

