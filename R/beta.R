#' Euclidean and related distances
#'
#' These distance and diversity measures are mathematically similar to the
#' Euclidean distance between two vectors.
#'
#' @param x,y Numeric vectors
#' @return The distance between \code{x} and \code{y}. The chord distance,
#'   Hellinger distance, and geodesic metric are not defined if all elements
#'   of either vector are zero. We return \code{NaN} in this case.
#' @details
#' For vectors \code{x} and \code{y}, the Euclidean distance is defined as
#' \deqn{d(x, y) = \sqrt{\sum_i (x_i - y_i) ^ 2}.}
#' Relation of \code{euclidean()} to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "euclidean"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "euclidean"}.
#'   \item Equivalent to the \code{euclidean()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to the \code{structeuclidean} calculator in Mothur, to
#'     \code{speciesprofile} if \code{x} and \code{y} are transformed to
#'     relative abundance, and to \code{memeuclidean} if \code{x} and \code{y}
#'     are transformed to presence/absence.
#'   \item Equivalent to \eqn{D_1}{D_1} in Legendre & Legendre.
#'   \item Equivalent to the \emph{distance between species profiles},
#'     \eqn{D_{18}}{D_18} in Legendre & Legendre if \code{x} and \code{y} are
#'     transformed to relative abundance.
#' }
#'
#' The \emph{root-mean-square} distance or \emph{average} distance is similar
#' to Euclidean distance. As the name implies, it is computed as the square
#' root of the mean of the squared differences between elements of \code{x}
#' and \code{y}: \deqn{d(x, y) = \sqrt{\frac{1}{n} \sum_i^n (x_i - y_i) ^ 2}.}
#' Relation of \code{rms_distance()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_2}{D_2} in Legendre & Legendre.
#' }
#'
#' The \emph{chord} distance is the Euclidean distance after scaling each
#' vector by its root sum of squares, \eqn{\sqrt{\sum_i x_i^2}}. The chord
#' distance between any two vectors ranges from 0 to \eqn{\sqrt{2}}{sqrt(2)}.
#' Relation of \code{chord()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_3}{D_3} in Legendre & Legendre.
#' }
#'
#' The \emph{Hellinger} distance is equal to the chord distance computed after
#' a square-root transformation. Relation of \code{hellinger()} to other
#' definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{17}}{D_17} in Legendre & Legendre.
#'   \item Equivalent to the \code{hellinger} calculator in Mothur.
#' }
#'
#' The \emph{geodesic metric} is a transformed version of the chord distance.
#' \deqn{d(x, y) = \textrm{arccos} \left(1 - \frac{d_c^2(x, y)}{2} \right),}
#' where \eqn{d_c} is the chord distance. It gives the length of the arc on a
#' hypersphere between the vectors, if the vectors are normalized to unit
#' length. Relation of \code{geodesic_metric()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_4}{D_4} in Legendre & Legendre.
#' }
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' y <- c(10, 2, 0, 1, 1, 0)
#' euclidean(x, y)
#' # The "distance between species profiles"
#' euclidean(x / sum(x), y / sum(y))
#' rms_distance(x, y)
#' chord(x, y)
#' hellinger(x, y)
#' # Hellinger is chord distance after square root transform
#' chord(sqrt(x), sqrt(y))
#' geodesic_metric(x, y)
#'
#' # No species in common with x
#' v <- c(0, 0, 0, 5, 0, 5)
#' chord(v, x)
#' sqrt(2)
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
  acos(1 - (chord(x, y) ^ 2) / 2)
}

#' Kullback-Leibler divergence
#'
#' @param x,y Numeric vectors representing probabilities
#' @return The Kullback-Leibler divergence between \code{x} and \code{y}. We
#'   adopt the following conventions if elements of \code{x} or \code{y} are
#'   zero: \eqn{0 \log (0 / y_i) = 0}, \eqn{0 \log (0 / 0) = 0}, and
#'   \eqn{x_i \log (x_i / 0) = \infty}. As a result, if elements of \code{x} are
#'   zero, they do not contribute to the sum. If elements of \code{y} are zero
#'   where \code{x} is nonzero, the result will be \code{Inf}. If either
#'   \code{x} or \code{y} sum to zero, we are not able to compute the
#'   proportions, and we return \code{NaN}.
#' @details
#' Kullback-Leibler divergence is a non-symmetric measure of difference between
#' two probability vectors. In general, KL(x, y) is not equal to KL(y, x).
#'
#' Because this measure is defined for probabilities, the vectors x and y are
#' normalized in the function so they sum to 1.
#' @export
kullback_leibler_divergence <- function (x, y) {
  x <- x / sum(x)
  y <- y / sum(y)
  terms <- x * log(x / y)
  # By convention, 0 * log(0 / y) = 0
  x_is_defined_and_zero <- (x == 0) %in% TRUE
  sum(ifelse(x_is_defined_and_zero, 0, terms))
}

#' Manhattan and related distances
#'
#' The Manhattan or city block distance is the sum of absolute differences
#' between the elements of two vectors. The \emph{mean character} difference
#' is a closely related measure.
#'
#' @param x,y Numeric vectors
#' @return The distance between \code{x} and \code{y}. The modified mean
#'   character difference is undefined if all elements in \code{x} and \code{y}
#'   are zero, in which case we return \code{NaN}.
#' @details
#' For vectors \code{x} and \code{y}, the Manhattan distance is given by
#' \deqn{d(x, y) = \sum_i |x_i - y_i|.} Relation of \code{manhattan()} to
#' other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "manhattan"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "manhattan"}.
#'   \item Equivalent to the \code{cityblock()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to the \code{manhattan} calculator in Mothur.
#'   \item Equivalent to \eqn{D_7}{D_7} in Legendre & Legendre.
#'   \item Whittaker's index of association (\eqn{D_9}{D_9} in Legendre &
#'     Legendre) is the Manhattan distance computed after transforming to
#'     proportions and dividing by 2.
#' }
#'
#' The mean character difference is the Manhattan distance divided by the
#' length of the vectors. It was proposed by Cain and Harrison in 1958.
#' Relation of \code{mean_character_difference()} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_8} in Legendre & Legendre.
#'   \item For binary data, equivalent to \eqn{1 - S_1} in Legendre & Legendre,
#'     where \eqn{S_1} is the simple matching coefficient.
#' }
#'
#' The modified mean character difference is the Manhattan distance divided by
#' the number elements where either \code{x} or \code{y} (or both) are nonzero.
#' Relation of \code{modified_mean_character_difference()} to other
#' definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{19}} in Legendre & Legendre.
#'   \item Equivalent to \code{vegdist()} with \code{method = "altGower"}.
#'   \item For binary data, it is equivalent to the Jaccard distance.
#' }
#' @references
#' Cain AJ, Harrison GA. An analysis of the taxonomist's judgment of affinity.
#' Proceedings of the Zoological Society of London 1958;131:85-98.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' y <- c(10, 2, 0, 1, 1, 0)
#' manhattan(x, y)
#' # Whittaker's index of association
#' manhattan(x / sum(x), y / sum(y)) / 2
#'
#' mean_character_difference(x, y)
#' # Simple matching coefficient for presence/absence data
#' # Should be 2 / 6
#' mean_character_difference(x > 0, y > 0)
#'
#' modified_mean_character_difference(x, y)
#' # Jaccard distance for presence/absence data
#' modified_mean_character_difference(x > 0, y > 0)
#' jaccard(x, y)
#' @export
manhattan <- function (x, y) {
  sum(abs(y - x))
}

#' @rdname manhattan
#' @export
mean_character_difference <- function (x, y) {
  manhattan(x, y) / length(x)
}

#' @rdname manhattan
#' @export
modified_mean_character_difference <- function (x, y) {
  pp <- sum((x > 0) | (y > 0))
  manhattan(x, y) / pp
}

#' Canberra and related distances
#'
#' The Canberra distance and Clark's coefficient of divergence are measures
#' that use the absolute difference over the sum for each element of the
#' vectors.
#'
#' @param x,y Numeric vectors
#' @return The Canberra distance or Clark's coefficient of divergence. If every
#'   element in \code{x} and \code{y} is zero, Clark's coefficient of
#'   divergence is undefined, and we return \code{NaN}.
#' @details
#' For vectors \code{x} and \code{y}, the Canberra distance is defined as
#' \deqn{d(x, y) = \sum_i \frac{|x_i - y_i|}{x_i + y_i}.} Elements where
#' \eqn{x_i + y_i = 0} are not included in the sum. Relation of
#' \code{canberra()} to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "canberra"}.
#'   \item Equivalent to the \code{vegdist()} function with
#'     \code{method = "canberra"}, multiplied by the number of entries where
#'     \code{x > 0}, \code{y > 0}, or both.
#'   \item Equivalent to the \code{canberra()} function in
#'     \code{scipy.spatial.distance} for positive vectors. They take the
#'     absolute value of \eqn{x_i} and \eqn{y_i} in the denominator.
#'   \item Equivalent to the \code{canberra} calculator in Mothur, multiplied
#'     by the total number of species in \code{x} and \code{y}.
#'   \item Equivalent to \eqn{D_{10}} in Legendre & Legendre.
#' }
#'
#' Clark's coefficient of divergence involves summing squares and taking a
#' square root afterwards:
#' \deqn{
#'   d(x, y) = \sqrt{
#'     \frac{1}{n} \sum_i \left( \frac{x_i - y_i}{x_i + y_i} \right)^2
#'   },}
#' where \eqn{n} is the number of elements where \code{x > 0}, \code{y > 0}, or
#' both. Relation of \code{clark_coefficient_of_divergence()} to other
#' definitions:
#' \itemize{
#'   \item Equivalent to \eqn{D_{11}}{D_11} in Legendre & Legendre.
#' }
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' y <- c(10, 2, 0, 1, 1, 0)
#' canberra(x, y)
#' clark_coefficient_of_divergence(x, y)
#' @export
canberra <- function (x, y) {
  xy_diff <- abs(x - y)
  xy_sum <- x + y
  keep <- xy_sum != 0
  sum(xy_diff[keep] / xy_sum[keep])
}

#' @rdname canberra
#' @export
clark_coefficient_of_divergence <- function (x, y) {
  xy_diff <- x - y
  xy_sum <- x + y
  keep <- xy_sum != 0
  pp <- sum(keep)
  sqrt(sum((xy_diff[keep] / xy_sum[keep]) ^ 2) / pp)
}


#' Chebyshev distance
#'
#' The Chebyshev distance is the maximum absolute difference between the vector
#' elements.
#'
#' @param x,y Numeric vectors
#' @return The Chebyshev distance between \code{x} and \code{y}.
#' @details
#' For vectors \code{x} and \code{y}, the Chebyshev distance is defined as
#' \deqn{d(x, y) = \max_i |x_i - y_i|.} Relation to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{chebyshev()} function in
#'     \code{scipy.spatial.distance}.
#' }
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' y <- c(10, 2, 0, 1, 1, 0)
#' chebyshev(x, y) # should be 5
#' @export
chebyshev <- function (x, y) {
  max(abs(x - y))
}

#' Correlation and cosine distance
#'
#' The correlation and cosine distances, which are derived from the dot
#' product of the two vectors.
#'
#' @param x,y Numeric vectors
#' @return The correlation or cosine distance. These are undefined if either
#'   \code{x} or \code{y} contain all zero elements, that is, if \eqn{|x| = 0}
#'   or \eqn{|y| = 0}. In this case, we return \code{NaN}.
#' @details
#' For vectors \code{x} and \code{y}, the cosine distance is defined as the
#' cosine of the angle between the vectors,
#' \deqn{d(x, y) = 1 - \frac{x \cdot y}{|x| |y|},} where \eqn{|x|} is the
#' magnitude or L2 norm of the vector, \eqn{|x| = \sqrt{\sum_i x_i^2}}.
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{cosine()} function in
#'     \code{scipy.spatial.distance}.
#' }
#'
#' The correlation distance is simply equal to one minus the Pearson
#' correlation between vectors. Mathematically, it is equivalent to the cosine
#' distance between the vectors after they are centered (\eqn{x - \bar{x}}).
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{correlation()} function in
#'     \code{scipy.spatial.distance}.
#'   \item Equivalent to the \code{1 - mempearson} calculator in Mothur.
#' }
#' @examples
#' x <- c(2, 0)
#' y <- c(5, 5)
#' cosine_distance(x, y)
#' # The two vectors form a 45 degree angle, or pi / 4
#' 1 - cos(pi / 4)
#'
#' v <- c(3.5, 0.1, 1.4)
#' w <- c(3.3, 0.5, 0.9)
#' correlation_distance(v, w)
#' 1 - cor(v, w)
#' @export
correlation_distance <- function (x, y) {
  x <- x - mean(x)
  y <- y - mean(y)
  cosine_distance(x, y)
}

#' @rdname correlation_distance
#' @export
cosine_distance <- function (x, y) {
  xy_dot <- sum(x * y)
  x_norm <- sqrt(sum(x ^ 2))
  y_norm <- sqrt(sum(y ^ 2))
  1 - xy_dot / (x_norm * y_norm)
}

#' Bray-Curtis distance
#'
#' The Bray-Curtis distance is the Manhattan distance divided by the sum of
#' both vectors.
#'
#' @param x,y Numeric vectors
#' @return The Bray-Curtis distance between \code{x} and \code{y}. The
#'   Bray-Curtis distance is undefined if the sum of all elements in \code{x}
#'   and \code{y} is zero, in which case we return \code{NaN}.
#' @details
#' For two vectors \code{x} and \code{y}, the Bray-Curtis distance is defined
#' as \deqn{d(x, y) = \frac{\sum_i |x_i - y_i|}{\sum_i x_i + y_i}.} The
#' Bray-Curtis distance is connected to many other distance measures in this
#' package; we try to list some of the more important connections here. Relation
#' to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "bray"}.
#'   \item Equivalent to the \code{braycurtis()} function in
#'     \code{scipy.spatial.distance} for positive vectors. They take the
#'     absolute value of \eqn{x_i + y_i} in the denominator.
#'   \item Equivalent to the \code{braycurtis} and \code{odum} calculators in
#'     Mothur.
#'   \item Equivalent to \eqn{D_{14} = 1 - S_{17}}{D_14 = 1 - S_17} in
#'     Legendre & Legendre.
#'   \item The Bray-Curtis distance on proportions is equal to half the
#'     Manhattan distance.
#'   \item The Bray-Curtis distance on presence/absence vectors is equal to the
#'     Sorenson index of dissimilarity.
#' }
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' y <- c(10, 2, 0, 1, 1, 0)
#' bray_curtis(x, y)
#'
#' # For proportions, equal to half the Manhattan distance
#' bray_curtis(x / sum(x), y / sum(y))
#' manhattan(x / sum(x), y / sum(y)) / 2
#' @export
bray_curtis <- function (x, y) {
  sum(abs(x - y)) / sum(x + y)
}

#' Weighted Kulczynski distance
#'
#' The quantitative version of the second Kulczynski index
#'
#' @param x,y Numeric vectors
#' @return The weighted Kulczynski distance between \code{x} and \code{y}. The
#'   distance is undefined if the sum of \code{x} or the sum of \code{y} is
#'   zero, in which case we return \code{NaN}.
#' @details
#' The quantitative version of the second Kulczynski index is defined as
#' \deqn{
#'   d(x, y) = 1 - \frac{1}{2} \left (
#'     \frac{\sum_i \min{(x_i, y_i)}}{\sum_i x_i} +
#'     \frac{\sum_i \min{(x_i, y_i)}}{\sum_i y_i}
#'   \right ).
#' }
#' Relation of \code{weighted_kulczynski_second()} to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "kulczynski"}.
#'   \item Equivalent to \code{structkulczynski} in Mothur.
#'   \item Equivalent to \eqn{1 - S_{18}} in Legendre & Legendre.
#' }
#' @export
weighted_kulczynski_second <- function (x, y) {
  xy_min <- sum(pmin(x, y))
  1 - (1 / 2) * (xy_min / sum(x) + xy_min / sum(y))
}

#' Minkowski distance
#'
#' The Minkowski metric is a generalized form of Euclidean (p=2) and Manhattan
#' (p=1) distance.
#'
#' @param x,y Numeric vectors.
#' @param p Exponent parameter, a single number greater than zero.
#' @return The Minkowski distance between \code{x} and \code{y}.
#' @details
#' For vectors \code{x} and \code{y}, the Minkowski distance is defined as
#' \deqn{d(x, y) = \left( \sum_i |x_i - y_i|^p \right)^{1/p}.} Relation to
#' other definitions:
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
  stopifnot(p > 0, length(p) == 1)
  sum(abs(x - y) ^ p) ^ (1 / p)
}

#' The Morisita index and Horn-Morisita index
#'
#' The Morisita and the Horn-Morisita indices measure the probability that
#' individuals drawn one from each vector will belong to different species,
#' relative to drawing from each vector separately. The Morisita index is
#' formulated for count data only, whereas the Horn-Morisita index can be
#' used with transformed counts or proportions.
#'
#' @param x,y Numeric vectors
#' @return The Morisita or Horn-Morisita index between \code{x} and \code{y}.
#'   Both are undefined if \code{x} or \code{y} have no nonzero elements, in
#'   which case we return \code{NaN}.
#' @details
#' For two vectors \code{x} and \code{y}, the Morisita index of dissimilarity
#' is
#' \deqn{d(x,y) = 1 - \frac{2 \sum_i x_i y_i}{(\lambda_x + \lambda_y) N_x N_y},}
#' where \deqn{\lambda_x = \frac{\sum_i x_i (x_i - 1)}{N_x (N_x - 1)}} and
#' \eqn{N_x = \sum_i x_i} The formula for \eqn{\lambda_x} is the unbiased
#' estimate for the probability of drawing two individuals of the same species
#' from \code{x}, without replacement. The correction for sampling without
#' replacement only makes sense for species count data.
#'
#' Relation of \code{morisita()} to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "morisita"}.
#' }
#'
#' Horn (1966) reformulated the index to use the equations for sampling with
#' replacement in \eqn{\lambda_x} and \eqn{\lambda_y}:
#' \deqn{\lambda_x = \frac{\sum_i x_i^2}{N_x^2}} With this modification,
#' the index is valid for proportions or transformed count data.
#'
#' Relation of \code{horn_morisita()} to other definitions:
#' \itemize{
#'   \item Equivalent to \code{vegdist()} with \code{method = "horn"}.
#'   \item Equivalent to the \code{morisitahorn} calculator in Mothur.
#' }
#' @references
#' Mosrisita M. Measuring of interspecific association and similarity between
#' communities. Memoirs of the Faculty of Science, Kyushu Univ., Series E
#' (Biology). 1959;3:65-80.
#'
#' Horn HS. Measurement of "Overlap" in Comparative Ecological Studies. The
#' American Naturalist, 1966;100(914):419-424.
#' @export
morisita <- function (x, y) {
  # Vegan docs not consistent with the paper, but implementation is correct
  Nx <- sum(x)
  Ny <- sum(y)
  # Eqn 1
  lambda_x <- sum(x * (x - 1)) / (Nx * (Nx - 1))
  # Eqn 2
  lambda_y <- sum(y * (y - 1)) / (Ny * (Ny - 1))
  # Eqn 3
  1 - 2 * sum(x * y) / ((lambda_x + lambda_y) * Nx * Ny)
}

#' @rdname morisita
#' @export
horn_morisita <- function (x, y) {
  Nx <- sum(x)
  Ny <- sum(y)
  # Horn defines new lambda estimates in middle of page 420
  lambda_x <- sum(x ^ 2) / (Nx ^ 2)
  lambda_y <- sum(y ^ 2) / (Ny ^ 2)
  1 - 2 * sum(x * y) / ((lambda_x + lambda_y) * Nx * Ny)
}


#' Binomial deviance and CY index of dissimilarity
#'
#' The binomial deviance dissimilarity and the CY (or Cao) index of
#' dissimilarity were created to compare species counts at sites with moderate
#' to large differences.
#'
#' @param x,y Numeric vectors
#' @param base Base of the logarithm
#' @param min_value Replacement for zero or near-zero values. Values less than
#'   \code{min_value} are replaced with \code{min_value}.
#' @return The Binomial deviance or CY index of dissimilarity. The CY index is
#'   undefined if all elements of \code{x} and \code{y} are zero, in which case
#'   we return \code{NaN}.
#' @details
#' Both of these measures were designed to be used with whole-numbered counts,
#' and may not make sense for comparing normalized vectors or vectors of
#' species proportions.
#'
#' For two vectors \code{x} and \code{y}, the binomial deviance dissimilarity
#' is
#' \deqn{
#'   d(x,y) = \sum_i{
#'     \frac{1}{n_i}
#'     \left (
#'       x_i \log{\frac{x_i}{n_i}} +
#'       y_i \log{\frac{y_i}{n_i}} -
#'       (x_i + y_i) log{2}
#'     \right )
#'   },
#' }
#' where \eqn{n_i = x_i + y_i}. This value is the weighted average of the
#' deviance for each species, under a binomial model where the expected counts
#' are \eqn{n_i / 2} at each site. It was proposed by Anderson and Millar in
#' 2004. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "binomial".
#' }
#'
#' The CY index was proposed by Cao, Williams, and Bark in 1997. For two
#' vectors \code{x} and \code{y}, the CY index is
#' \deqn{
#'   d(x,y) = \frac{1}{N} \sum_i
#'   \left (
#'     \frac{
#'       (x_i + y_i) \log_{10} ( \frac{x_i + y_i}{2} ) -
#'       x_i \log_{10}(y_i) - y_i \log_{10}(x_i)
#'     }{
#'       x_i + y_i
#'     }
#'   \right ),
#' }
#' where \eqn{N} is the total number of species in vectors \eqn{x} and \eqn{y}.
#' Double zeros are not considered in the measure.
#'
#' When either \eqn{x_i} or \eqn{y_i} are zero, they need to be replaced by
#' another value in the CY index to avoid infinities. Cao suggested replacing
#' zero values with \eqn{0.1}, which is one log lower than the minimum value
#' for whole-numbered counts. Here, we use a \code{min_value} argument to allow
#' the user set a lower limit on the values. For vectors of species counts,
#' this function follows the formulation of Cao by default.
#'
#' Relation of the CY index to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{vegdist()} function with
#'     \code{method = "cao"}, if \code{base = exp(1)}.
#' }
#' @references
#' Anderson MJ, Millar RB. Spatial variation and effects of habitat on
#' temperate reef fish assemblages in northeastern New Zealand. Journal of
#' Experimental Marine Biology and Ecology 2004;305:191–221.
#'
#' Cao Y, Williams WP, Bark AW. Similarity measure bias in river benthic
#' Aufwuchs community analysis. Water Environment Research 1997;69(1):95-106.
#' @export
binomial_deviance <- function (x, y) {
  x_is_defined_and_zero <- (x == 0) %in% TRUE
  y_is_defined_and_zero <- (y == 0) %in% TRUE
  n <- x + y
  n_is_defined_and_zero <- x_is_defined_and_zero & y_is_defined_and_zero
  # Formula at top of page 199
  t1 <- ifelse(x_is_defined_and_zero, 0, x * log(x / n))
  t2 <- ifelse(y_is_defined_and_zero, 0, y * log(y / n))
  t3 <- n * log(2)
  # BUG IN VEGAN?!?!?!
  # Formula in paper subtracts t3, vegan function adds this term
  terms <- ifelse(n_is_defined_and_zero, 0, (1 / n) * (t1 + t2 + t3))
  sum(terms)
}

#' @rdname binomial_deviance
#' @export
cy_dissimilarity <- function (x, y, base = 10, min_value = 0.1) {
  # Remove double zeros
  keep <- (x > 0) | (y > 0)
  x <- x[keep]
  y <- y[keep]
  # Substitute individual zeros with 0.1, just like vegan
  x <- ifelse(x > min_value, x, min_value)
  y <- ifelse(y > min_value, y, min_value)
  N <- length(x)
  xy_sum <- x + y
  t1 <- xy_sum * log(xy_sum / 2, base = base)
  t2 <- x * log(y, base = base)
  t3 <- y * log(x, base = base)
  # Formula 12 in the Cao paper
  (1 / N) * sum((t1 - t2 - t3) / xy_sum)
}

#' Ruzicka or weighted Jaccard distance
#'
#' @param x,y Numeric vectors.
#' @return The Ruzicka distance between \code{x} and \code{y}. The distance is
#'   not defined if all elements in \code{x} and \code{y} are zero, and we
#'   return \code{NaN} in this case.
#' @details
#' For vectors \code{x} and \code{y}, the Ruzicka distance is defined as
#' \deqn{d(x, y) = 1 - \frac{\sum_i \min(x, y)}{\sum_i \max(x, y)}.} Relation
#' to other definitions:
#' \itemize{
#'   \item Equivalent to vegdist() with method = "jaccard".
#'   \item Related to the Bray-Curtis distance,
#'     \eqn{d_r = 2 d_{bc} / (1 + d_{bc})}.
#' }
#' @export
ruzicka <- function (x, y) {
  1 - sum(pmin(x, y)) / sum(pmax(x, y))
}

#' Beta diversity for presence/absence data
#'
#' These functions transform the input vectors to binary or presence/absence
#' format, then compute a distance or dissimilarity.
#'
#' @param x,y Numeric vectors
#' @return The dissimilarity between \code{x} and \code{y}, based on
#'   presence/absence. The Jaccard, Sorenson, Sokal-Sneath, Yule, and both
#'   Kulczynski dissimilarities are not defined if both \code{x} and \code{y}
#'   have no nonzero elements. In addition, the second Kulczynski index and the
#'   Yule index of dissimilarity are not defined if one of the vectors has no
#'   nonzero elements. We return \code{NaN} for undefined values.
#' @details
#' Many of these indices are covered in Koleff et al. (2003), so we adopt their
#' notation. For two vectors \code{x} and \code{y}, we define three quantities:
#' \itemize{
#'   \item \eqn{a} is the number of species that are present in both \code{x}
#'     and \code{y},
#'   \item \eqn{b} is the number of species that are present in \code{y} but
#'     not \code{x},
#'   \item \eqn{c} is the number of species that are present in \code{x} but
#'     not \code{y}, and
#'   \item \eqn{d} is the number of species absent in both vectors.
#' }
#' The quantity \eqn{d} is seldom used in ecology, for good reason. For
#' details, please see the discussion on the "double zero problem," in section
#' 2 of chapter 7.2 in Legendre & Legendre.
#'
#' The \emph{Jaccard} index of dissimilarity is \eqn{1 - a / (a + b + c)}, or
#' one minus the proportion of shared species, counting over both samples
#' together. Relation of \code{jaccard()} to other definitions:
#' \itemize{
#'   \item Equivalent to R's built-in \code{dist()} function with
#'     \code{method = "binary"}.
#'   \item Equivalent to \code{vegdist()} with \code{method = "jaccard"}
#'     and \code{binary = TRUE}.
#'   \item Equivalent to the \code{jaccard()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to \eqn{1 - S_7} in Legendre & Legendre.
#'   \item Equivalent to \eqn{1 - \beta_j}, as well as \eqn{\beta_{cc}}, and
#'     \eqn{\beta_g} in Koleff (2003).
#' }
#'
#' The \emph{\enc{Sørenson}{Sorenson}} or \emph{Dice} index of dissimilarity is
#' \eqn{1 - 2a / (2a + b + c)}, or one minus the average proportion of shared
#' species, counting over each sample individually. Relation of
#' \code{sorenson()} to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{dice()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to the \code{sorclass} calculator in Mothur, and to
#'     \code{1 - whittaker}.
#'   \item Equivalent to \eqn{D_{13} = 1 - S_8}{D_13 = 1 - S_8} in Legendre &
#'     Legendre.
#'   \item Equivalent to \eqn{1 - \beta_{sor}} in Koleff (2003). Also
#'     equivalent to Whittaker's beta diversity
#'     (the second definition, \eqn{\beta_w = (S / \bar{a}) - 1}), as well as
#'     \eqn{\beta_{-1}}, \eqn{\beta_t}, \eqn{\beta_{me}}, and
#'     \eqn{\beta_{hk}}.
#' }
#'
#' I have not been able to track down the original reference for the first and
#' second Kulczynski indices, but we have good formulas from Legendre &
#' Legendre. The \emph{first Kulczynski index} is \eqn{1 - a / (b + c)}, or
#' one minus the ratio of shared to unshared species.
#'
#' Relation of \code{kulczynski_first} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{1 - S_{12}}{1 - S_12} in Legendre & Legendre.
#'   \item Equivalent to the \code{kulczynski} calculator in Mothur.
#' }
#'
#' Some people refer to the \emph{second Kulczynski index} as the
#' Kulczynski-Cody index. It is defined as one minus the average proportion of
#' shared species in each vector,
#' \deqn{
#'   d = 1 - \frac{1}{2} \left ( \frac{a}{a + b} + \frac{a}{a + c} \right ).
#' }
#' Relation of \code{kulczynski_second} to other definitions:
#' \itemize{
#'   \item Equivalent to \eqn{1 - S_{13}}{1 - S_13} in Legendre & Legendre.
#'   \item Equivalent to the \code{kulczynskicody} calculator in Mothur.
#'   \item Equivalent to one minus the Kulczynski similarity in Hayek (1994).
#'   \item Equivalent to \code{vegdist()} with \code{method = "kulczynski"} and
#'     \code{binary = TRUE}.
#' }
#'
#' The \emph{Rogers-Tanimoto} distance is defined as
#' \eqn{(2b + 2c) / (a + 2b + 2c + d)}. Relation of \code{rogers_tanimoto()}
#' to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{rogerstanimoto()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to \eqn{1 - S_2}{1 - S_2} in Legendre & Legendre.
#' }
#'
#' The \emph{Russel-Rao} distance is defined
#' \eqn{(b + c + d) / (a + b + c + d)}, or the fraction of elements not present
#' in both vectors, counting double absences. Relation of \code{russel_rao()} to
#' other definitions:
#' \itemize{
#'   \item Equivalent to the \code{russelrao()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to \eqn{1 - S_{11}}{1 - S_11} in Legendre & Legendre.
#' }
#'
#' The \emph{Sokal-Michener} distance is defined as
#' \eqn{(2b + 2c) / (a + 2b + 2c + d)}. Relation of \code{sokal_michener()} to
#' other definitions:
#' \itemize{
#'   \item Equivalent to the \code{sokalmichener()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#' }
#'
#' The \emph{Sokal-Sneath} distance is defined as
#' \eqn{(2b + 2c) / (a + 2b + 2c)}. Relation of \code{sokal_sneath()} to other
#' definitions:
#' \itemize{
#'   \item Equivalent to the \code{sokalsneath()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to the \code{anderberg} calculator in Mothur.
#'   \item Equivalent to \eqn{1 - S_{10}}{1 - S_10} in Legendre & Legendre.
#' }
#'
#' The \emph{Yule} dissimilarity is defined as \eqn{2bc / (ad + bc)}. Relation
#' of \code{yule_dissimilarity()} to other definitions:
#' \itemize{
#'   \item Equivalent to the \code{yule()} function in
#'     \code{scipy.spatial.distance}, except that we always convert vectors to
#'     presence/absence.
#'   \item Equivalent to \eqn{1 - S}, where \eqn{S} is the Yule coefficient
#'     in Legendre & Legendre.
#' }
#' @export
jaccard <- function (x, y) {
  x <- x > 0
  y <- y > 0
  sum(xor(x, y)) / sum(x | y)
}

#' @rdname jaccard
#' @export
sorenson <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  bc / (2 * a + bc)
}

#' @rdname jaccard
#' @export
kulczynski_first <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  1 - a / bc
}

#' @rdname jaccard
#' @export
kulczynski_second <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  1 - 0.5 * (a / sum(x) + a / sum(y))
}

#' @rdname jaccard
#' @export
rogers_tanimoto <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  d <- sum((!x) & (!y))
  2 * bc / (a + d + 2 * bc)
}

#' @rdname jaccard
#' @export
russel_rao <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  n <- length(x)
  (n - a) / n
}

#' @rdname jaccard
#' @export
sokal_michener <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  d <- sum((!x) & (!y))
  2 * bc / (a + d + 2 * bc)
}

#' @rdname jaccard
#' @export
sokal_sneath <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  bc <- sum(xor(x, y))
  2 * bc / (a + 2 * bc)
}

#' @rdname jaccard
#' @export
yule_dissimilarity <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  b <- sum((!x) & y)
  c <- sum(x & (!y))
  d <- sum((!x) & (!y))
  2 * b * c / (a * d + b * c)
}

#' Hamming distance
#'
#' The Hamming distance is the number of positions where the values are
#' different.
#'
#' @param x,y Numeric vectors.
#' @return The Hamming distance between \code{x} and \code{y}.
#' @details
#' For vectors \code{x} and \code{y}, the Hamming distance is defined as
#' \deqn{d(x, y) = \sum_i [x_i \neq y_i],} where the quantity in the brackets
#' is 1 if the elements are not equal, and zero if the elements are equal.
#' Relation to other definitions:
#' \itemize{
#'   \item The \code{hamming()} function in \code{scipy.spatial.distance}
#'     divides the result by the vector length. Our function is equivalent to
#'     the SciPy version multiplied by the vector length.
#'   \item Equivalent to the \code{hamming} calculator in Mothur for
#'     presence/absence vectors.
#' }
#' @export
hamming <- function (x, y) {
  sum(x != y)
}

#' Chao's abundance-weighted indices
#'
#' These indices were first developed and introduced as probabilistic,
#' abundance-weighted versions of the Jaccard and Sorenson indices. The
#' abundance information is summarized by two quantities, named \eqn{U} and
#' \eqn{V} in the literature. Each quantity indicates the total abundance of
#' species that are common to both samples: \eqn{U} is the total abundance of
#' such species in the first sample, and \eqn{V} the total in the second sample.
#'

#' @param x,y Numeric vectors
#' @return The dissimilarity between \code{x} and \code{y}
#' @details
#' In addition to defining the quantities \eqn{U} and \eqn{V}, the articles
#' provide formulas for estimating the value \eqn{U} and \eqn{V} if more
#' sampling was to be performed. In other software packages where these indices
#' are implemented, it is the estimated quantities, \eqn{U_{est}} and
#' \eqn{V_{est}} that are used to calculate the distance. For example, the
#' \code{chaodist} function in Vegan uses the estimators. In this package, we
#' use the observed values of \eqn{U} and \eqn{V}, rather than the estimated
#' values. Thus, the distances returned here may differ from other software
#' packages that use the estimators rather than the observed values of \eqn{U}
#' and \eqn{V}. In general, the estimators increase \eqn{U} and \eqn{V},
#' resulting in greater estimated abundance of shared species and thus a
#' smaller distance between samples.
#'
#' The \emph{abundance-weighted Jaccard} distance is \eqn{UV / (U + V - UV)}.
#' Relation of \code{abundance_jaccard()} to other definitions:
#' \itemize{
#'   \item Differs \code{vegdist()} with \code{method = "chao"}. We use the
#'     observed values, not the estimators.
#'   \item Differs from the \code{jabund} calculator in Mothur. We use the
#'     observed values, not the estimators.
#' }
#'
#' The \emph{abundance-weighted Sorenson} distance is \eqn{2UV / (U + V)}.
#' Relation of \code{abundance_sorenson()} to other definitions:
#' \itemize{
#'   \item Differs from the \code{sabund} calculator in Mothur. We use the
#'     observed values, not the estimators.
#' }
#' @references
#' Chao A, Chazdon RL, Colwell R, Shen TJ. A new statistical approach for
#' assessing similarity of species composition with incidence and abundance
#' data. Ecology Letters 2005;8:148-159.
#' @export
abundance_jaccard <- function (x, y) {
  shared <- (x > 0) & (y > 0)
  x <- x / sum(x)
  y <- y / sum(y)
  u <- sum(x[shared])
  v <- sum(y[shared])
  1 - u * v / (u + v - u * v)
}

#' @rdname abundance_jaccard
#' @export
abundance_sorenson <- function (x, y) {
  shared <- (x > 0) & (y > 0)
  x <- x / sum(x)
  y <- y / sum(y)
  u <- sum(x[shared])
  v <- sum(y[shared])
  1 - 2 * u * v / (u + v)
}

# Scipy notes:
# # Distance functions for continuous vectors
# braycurtis implemented as bray_curtis
# canberra implemented
# chebyshev implemented
# cityblock implemented as manhattan
# correlation implemented as correlation_distance
# cosine implemented
# euclidean implemented
# jensenshannon not implemented
# mahalanobis not implemented, available in R stats module
# minkowski implemented
# seuclidean not implemented
# sqeuclidean not implemented
# wminkowski not implemented
# # Distance functions for presence/absence vectors
# *** We convert the vectors to presence/absence, so our results do not match
# *** the SciPy tests when elements are greater than 1. The SciPy
# *** implementation is arguably wrong here.
# dice implemented as sorenson
# hamming implemented differently
# jaccard implemented
# kulsinski not implemented, not sure measure is correct in SciPy
# rogerstanimoto implemented as rogers_tanimoto
# russelrao implemented as russel_rao
# sokalmichener implemented as sokal_michener
# sokalsneath implemented as sokal_sneath
# yule implemented as yule_dissimilarity

# Mothur notes
# https://www.mothur.org/wiki/Calculators
# ## Community membership
# anderberg implemented as sokal_sneath
# hamming implemented as hamming for presence/absence vectors
# jclass implemented as jaccard (Mothur docs give formula for similarity)
# jest (Jaccard with Chao estimates) not implemented
# kulczynski implemented as kulczynski_first
# kulczynskicody implemented as kulczynski_second
# lennon not implemented, TODO
# memchi2 not implemented, needs full matrix
# memchord not implemented, Mothur docs seem off
#   If the terms in the denominator were squared, this would be the chord
#   distance for presence/absence data. Must check source code.
# memeuclidean equivalent to euclidean on presence/absence vectors
# mempearson equivalent to 1 - correlation_distance
# ochiai not implemented, TODO
# sorclass implemented as sorenson
# sorest not implemented
# whittaker equivalent to 1 - sorenson
# ## Community structure
# braycurtis implemented as bray_curtis
# canberra implemented, see note
# gower not implemented, needs full matrix
# hellinger implemented
# jabund (Chao-Jaccard or abundance-Jaccard) not implemented, TODO
# manhattan implemented
# morisitahorn implemented as horn_morisita
# odum is equivalent to bray_curtis
# soergel not implemented, TODO
# sorabund (Chao-Sorenson or abundance-Sorenson) not implemented, TODO
# spearman not implemented
# speciesprofile equivalent to Euclidean on proportions, noted
# structchi2 not implemented, needs full matrix
# structchord implemented as chord
# structeuclidean implemented as euclidean
# structkulczynski implemented as weighted_kulczynski_second
# structpearson implemented as correlation_distance
# thetan (community Jaccard index of Smith aka the Θn of Yue) not implemented,
#   TODO
# thetayc (Yue & Clayton measure) not implemented, TODO

# Koleff notes:
#  1. \beta_w (Whittaker's beta diversity) implemented as sorenson
#  2. \beta_{-1} implemented as sorenson
#  3. \beta_c not implemented
#  4. \beta_{wb} (Weiher & Boylen) not implemented
#  5. \beta_r (Routledge) not implemented
#  6. \beta_I not implemented
#  7. \beta_e not implemented
#  8. \beta_t implemented as sorenson
#  9. \beta_{me} implemented as sorenson
# 10. \beta_j implemented as jaccard
# 11. \beta_{sor} implemented as sorenson
# 12. \beta_m (Magurran) not implemented
# 13. \beta_{-2} not implemented
# 14. \beta_{co} not implemented
# 15. \beta_{cc} implemented as jaccard
# 16. \beta_g implemented as jaccard
# 17. \beta_{-3} not implemented
# 18. \beta_1 not implemented
# 19. (Williams) not implemented
# 20. \beta_{hk} (Harte & Kinzig) implemented as sorenson
# 21. \beta_{rlb} not implemented
# 22. \beta_{sim} not implemented
# 23. \beta_{gl} not implemented
# 24. \beta_z not implemented

# Vegan notes:
# # TODO: document these functions better
# # Methods in vegdist, listed in order of the METHODS vector
#  1. manhattan implemented
#     manhattan with binary = TRUE not implemented
#  2. euclidean implemented
#     euclidean with binary = TRUE not implemented
#  3. canberra implemented
#     canberra with binary = TRUE not implemented
#  4. bray implemented as bray_curtis
#     bray with binary = TRUE not implemented
#  5. kulczynski implemented as weighted_kulczynski_second
#     kulczynski with binary = TRUE implemented as kulczynski_second
#  6. gower not implemented, needs full matrix
#     gower with binary = TRUE not implemented
#  7. morisita implemented
#     morisita with binary = TRUE can't be calculated
#  8. horn implemented as horn_morisita
#     horn with binary = TRUE not implemented
#  9. mountford not imlemented, TODO
# 10. jaccard implemented as ruzicka
#     jaccard with binary = TRUE implemented as jaccard
# 11. raup no implemented, TODO
# 12. binomial implemented as millar (from Anderson & Millar 2004)
# 13. chao implemented as abundance_jaccard, but we implement the *observed*
#     distance, not the estimated distance
#.    chao with binary = TRUE not implemented
# 14. altGower implemented as modified_mean_character_difference
#     altGower with binary = TRUE not implemented
# 15. cao implemented
#     cao with binary = TRUE not implemented
# 16. mahalanobis not implemented, available in R stats module
#     mahalanobis with binary = TRUE not implemented
# 17. clark not implemented?
#     clark with binary = TRUE not implemented
# 18. chisq not implemented, needs full matrix
# 19. chord implemented (must test)
#     chord with binary = TRUE not implemented
# 20. hellinger implemented (must test)
# 21. aitchison not implemented, TODO
# 22. robust.aitchison not implemented, TODO


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
# D_18 (distance between species profiles) not implemented,
#   noted under euclidean
# D_15 (Chi-square metric) not implemented, needs full matrix
# D_16 (Chi-square distance) not implemented, needs full matrix
# D_17 implemented as hellinger
# D_13 implemented as sorenson
# D_14 implemented as bray_curtis
# S_1 (simple matching coefficient) is the mean_character_difference
# S_2 (coefficient of Rogers & Tanimoto) implemented as rogers_tanimoto
# S_3, S_4, S_5, S_6 (Sokal and Sneath) not implemented
# Hamann coefficient not implemented
# Yule coefficient implemented as yule_dissimilarity
# Pearson's phi not implemented
# S_7 implemented as jaccard
# S_8 implemented as sorenson
# S_9 = 3a / (3a + b + c) not implemented
# S_10 implemented as sokal_sneath
# S_11 implemented as russel_rao
# S_12 implemented as kulczynski_first
# S_13 implemented as kulczynski_second
# S_14 (name unclear) not implemented. It should be noted that S_14 as a
#   distance is proportional to sqrt of chord or Hellinger for binary data.
# S_26 (name unclear) not implemented
# S_15 (Gower coefficient) not implemented. Equivalent to Hamming distance?
# S_16 (Estabrook and Rogers) not implemented
# S_17 (Steinhaus coefficient) implemented as bray_curtis
# S_18 (quantitative Kulcynski) implemented as weighted_kulczynski_second
# S_19 (modified Gower) not implemented
# S_20 (Legendre & Chodorowski) not implemented, also similar to Gower
# S_21 (Chi-square similarity) not implemented, see D_15
# S_22 (Goodall 1) not implemented, needs full matrix
# S_23 (Goodall 2) not implemented, needs full matrix
# Raup-Crick not implemented, available as raupcrick() in vegan
# S_27 (Raup-Crick p-value) not implemented
# Skip rest of Q-mode similarity for now

