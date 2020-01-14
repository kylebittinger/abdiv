#' Berger-Parker dominance
#'
#' The Berger-Parker dominance is the proportion of the most abundant species.
#'
#' @param x A numeric vector of species counts or proportions.
#' @return The Berger-Parker dominance, \eqn{0 < D_{BP} \leq 1}. If the vector
#'   sums to zero, the Berger-Parker dominance is undefined, and we return
#'   \code{NaN}.
#' @details
#' \itemize{
#'   \item Equivalent to \code{berger_parker_d()} in
#'     \code{skbio.diversity.alpha}.
#'   \item Equivalent to the \code{bergerparker} calculator in Mothur.
#' }
#' @references
#' Berger WH, Parker FL. Diversity of Planktonic Foraminifera in Deep-Sea
#' Sediments. Science. 1970;168(3937):1345-1347.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' berger_parker_d(x) # 15 / 28
#' @export
berger_parker_d <- function (x) {
  max(x) / sum(x)
}

#' Simpson's index and related measures
#'
#' These measures are based on the sum of squared species proportions. The
#' function \code{dominance()} gives this quantity, \code{simpson()} gives one
#' minus this quantity, \code{invsimpson()} gives the reciprocal of the
#' quantity, and \code{simpson_e} gives the reciprocal divided by the number
#' of species.
#'
#' @param x A numeric vector of species counts or proportions.
#' @return The value of the dominance (\eqn{0 < D \leq 1}), Simpson index, or
#'   inverse Simpson index. The dominance is undefined if the vector sums to
#'   zero, in which case we return \code{NaN}.
#' @details
#' For a vector of species counts \code{x}, the dominance index is defined as
#' \deqn{D = \sum_i p_i^2,} where \eqn{p_i} is the species proportion,
#' \eqn{p_i = x_i / N}, and \eqn{N} is the total number of counts. This is
#' equal to the probability of selecting two individuals from the same species,
#' with replacement. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{dominance()} in \code{skbio.diversity.alpha}.
#'   \item Similar to the \code{simpson} calculator in Mothur. They use the
#'     unbiased estimate \eqn{p_i = x_i (x_i - 1) / (N (N -1))}.
#' }
#'
#' Simpson's index is defined here as \eqn{1 - D}, or the probability of
#' selecting two individuals from different species, with replacement. Relation
#' to other definitions:
#' \itemize{
#'   \item Equivalent to \code{diversity()} in \code{vegan} with
#'     \code{index = "simpson"}.
#'   \item Equivalent to \code{simpson()} in \code{skbio.diversity.alpha}.
#' }
#'
#' The inverse Simpson index is \eqn{1/D}. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{diversity()} in \code{vegan} with
#'     \code{index = "invsimpson"}.
#'   \item Equivalent to \code{enspie()} in \code{skbio.diversity.alpha}.
#'   \item Similar to the \code{invsimpson} calculator in Mothur. They use
#'     the unbiased estimate \eqn{p_i = x_i (x_i - 1) / (N (N -1))}.
#' }
#'
#' Simpson's evenness index is the inverse Simpson index divided by the
#' number of species observed, \eqn{1 / (D S)}. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{simpson_e()} in \code{skbio.diversity.alpha}.
#' }
#'
#' Please be warned that the naming conventions vary between sources. For
#' example Wikipedia calls \eqn{D} the Simpson index and \eqn{1 - D} the
#' Gini-Simpson index. We have followed the convention from \code{vegan}, to
#' avoid confusion within the \code{R} ecosystem.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' dominance(x)
#'
#' # Simpson is 1 - D
#' simpson(x)
#' 1 - dominance(x)
#'
#' # Inverse Simpson is 1/D
#' invsimpson(x)
#' 1 / dominance(x)
#'
#' # Simpson's evenness is 1 / (D * S)
#' simpson_e(x)
#' 1 / (dominance(x) * richness(x))
#' @export
simpson <- function (x) {
  p <- x / sum(x)
  1 - sum(p ** 2)
}

#' @rdname simpson
#' @export
dominance <- function (x) {
  p <- x / sum(x)
  sum(p ** 2)
}

#' @rdname simpson
#' @export
invsimpson <- function (x) {
  p <- x / sum(x)
  1 / sum(p ** 2)
}

#' @rdname simpson
#' @export
simpson_e <- function (x) {
  p <- x / sum(x)
  D <- sum(p ** 2)
  S <- sum(x > 0)
  1 / (D * S)
}

#' Kempton-Taylor Q index
#'
#' The Kempton-Taylor Q index is designed to measure species in the middle of
#' the abundance distribution.
#'
#' @param x A numeric vector of species counts or proportions.
#' @param lower_quantile,upper_quantile Lower and upper quantiles of the
#'   abundance distribution. Default values are the ones suggested by Kempton
#'   and Taylor.
#' @return The Kempton-Taylor Q index, \eqn{Q < 0}. If the vector sums to zero,
#'   we cannot compute the quantiles, and this index is undefined. In that
#'   case, we return \code{NaN}.
#' @details
#' For a vector of species counts \code{x}, the Kempton-Taylor Q statistic is
#' equal to the slope of the cumulative abundance curve across a specified
#' quantile range. The cumulative abundance curve is the plot of the number of
#' species against the log-abundance.
#'
#' Kempton and Taylor originally defined the index as
#' \deqn{Q = \frac{\frac{1}{2}S}{\log{R_2} - \log{R_1}},} where \eqn{S} is the
#' total number of species observed, \eqn{R_1} is the abundance at the lower
#' quantile, and \eqn{R_2} is the abundance at the upper quantile. However,
#' this definition only holds if one uses the interquartile range. Because we
#' allow the user to adjust the upper and lower quantiles, we have to find the
#' number of species at these abundance values. Here, we follow the
#' implementation in \code{scikit-bio} and round inwards to find the quantile
#' values, taking the number of species and log-abundance values at these data
#' points exactly.
#'
#' \itemize{
#'   \item Equivalent to \code{kempton_taylor_q()} in
#'     \code{skbio.diversity.alpha}.
#'   \item Similar to the \code{qstat} calculator in Mothur. Our implementation
#'     differs slightly, and this difference affects the result.
#' }
#' @references
#' Kempton RA, Taylor LR. Models and statistics for species diversity. Nature.
#' 1976;262:818-820.
#' @export
kempton_taylor_q <- function (x, lower_quantile=0.25, upper_quantile=0.75) {
  # I'm sure there is a better way to do this with R's quantile function,
  # but not sure how to guarantee that the result always replicates the one
  # obtained via this algorithm.
  # !!!! The scikit-bio implementation seems to not match the paper here.
  # Must check this and make sure the answers line up with their results.
  n <- length(x)
  lower_idx <- ceiling(n * lower_quantile) + 1
  upper_idx <- floor(n * upper_quantile) + 1
  x_sorted <- sort(x)
  x_upper <- x_sorted[upper_idx]
  x_lower <- x_sorted[lower_idx]
  (upper_idx - lower_idx) / log(x_upper / x_lower)
}

#' Margalef's richness index
#'
#' @param x A numeric vector of species counts.
#' @return The value of Margalef's index, \eqn{D \geq 0}. This index is
#'   undefined when the total number of counts is 1 or 0, in which case we
#'   return \code{NaN}.
#' @details
#' For a vector \code{x} of species counts, Margalef's index is
#' \deqn{D = \frac{S -1}{\log N},} where \eqn{S} is the total number of species
#' observed and \eqn{N} is the total number of counts.
#'
#' This index is appropriate only for raw counts, not transformed counts or
#' proportions.
#'
#' Equivalent to \code{margalef()} in \code{skbio.diversity.alpha}.
#' @references
#' Margalef R. Information theory in ecology. General Systems 3. 1958;36-71.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' margalef(x)
#' @export
margalef <- function (x) {
  # Margalef is based on the slope of the species-area curve as proposed by
  # Gleason (1922), where the number of species increases with the log of the
  # area. Here, the number of individuals, n, is a stand-in for the area. The
  # equation is:
  #   s = s0 + z * log(N)
  # When the number of individuals, n, is one, then z * log(N) is zero. In this
  # case, s is one by definition, so s0 must be one.  Making this substitution
  # and solving for z, we get
  #   z = (s - 1) / log(N)
  # When we only observe one individual, we don't know the slope. Numerically,
  # we encounter zero divided by zero, so we expect NaN as a result. This is
  # what R returns.
  # When we observe no individuals, we also dont't know the slope. Numerically,
  # we encounter negative one over negative infinity, which R interprets as
  # zero. However, if we return to the original equation, we see that z must be
  # a number that multiplies negative infinity to produce negative one. Zero
  # does not fulfill this role, and is therefore not an acceptable answer.
  #   0 = 1 + z * (-Inf)
  # The quantity z is undefined at n = 1, and we also consider it to be
  # undefined at n = 0. Therefore, we take special care to return NaN when
  # n = 0.
  s <- sum(x > 0)
  n <- sum(x)
  if (n == 0) {
    NaN
  } else {
    (s - 1) / log(n)
  }
}

#' McIntosh dominance index D
#' @param x A numeric vector of species counts.
#' @return The McIntosh dominance index, \eqn{0 \leq D < 1}. The index is undefined
#'   when the total number of counts is 1 or 0, in which case we return
#'   \code{NaN}.
#' @details
#' For a vector \code{x} of raw species counts, the McIntosh dominance index is
#' defined as \deqn{D = \frac{N - U}{N - \sqrt{N}},} where \eqn{N} is the total
#' number of counts and \eqn{U = \sqrt{\sum_i x_i^2}}.
#'
#' This index is appropriate only for raw counts, not transformed counts or
#' proportions.
#'
#' Equivalent to \code{mcintosh_d()} in \code{skbio.diversity.alpha}.
#' @references
#' McIntosh RP. An index of diversity and the relation of certain concepts to
#' diversity. Ecology. 1967;48:1115-1126.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' mcintosh_d(x)
#' @export
mcintosh_d <- function (x) {
  # Equation 4
  n <- sum(x)
  u <- sqrt(sum(x ^ 2))
  (n - u) / (n - sqrt(n))
}

#' McIntosh's evenness measure E
#' @param x A numeric vector of species counts.
#' @return McIntosh's evenness measure, \eqn{0 < E \leq 1}.  The index is
#'   undefined when the total number of counts is 0, in which case we return
#'   \code{NaN}.
#' @details
#' For a vector \code{x} of raw species counts, the McIntosh evenness measure
#' is \deqn{E = \frac{\sqrt{\sum_i x_i^2}}{\sqrt{(N - S + 1)^2 + S - 1},}}
#' where \eqn{N} is the total number of counts and \eqn{S} is the total
#' number of species observed.
#'
#' This index is appropriate only for raw counts, not transformed counts or
#' proportions.
#'
#' Equivalent to \code{mcintosh_e()} in \code{skbio.diversity.alpha}.
#' @references
#' Heip C, Engels P. Comparing Species Diversity and Evenness Indices. J. Mar.
#' Bioi. Ass. U.K. 1974;54:559-563.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' mcintosh_e(x)
#' @export
mcintosh_e <- function (x) {
  n <- sum(x)
  s <- sum(x > 0)
  numerator <- sqrt(sum(x ^ 2))
  denominator <- sqrt((n - s + 1) ^ 2 + s - 1)
  numerator / denominator
}

#' Menhinick's richness index
#' @param x A numeric vector of species counts.
#' @return Menhinick's richness index, \eqn{R > 0}. The index is undefined when
#'   the total number of counts is 0, in which case we return \code{NaN}.
#' @details
#' For a vector \code{x} of raw species counts, the Menhinick's richness index
#' is \eqn{\frac{S}{\sqrt{N}}}, where \eqn{N} is the total number
#' of counts and \eqn{S} is the total number of species observed.
#'
#' This index is appropriate only for raw counts, not transformed counts or
#' proportions.
#'
#' Equivalent to \code{menhinick()} in \code{skbio.diversity.alpha}.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' menhinick(x)
#' @export
menhinick <- function (x) {
  n <- sum(x)
  s <- sum(x > 0)
  s / sqrt(n)
}

#' Richness or number of observed species
#' @param x A numeric vector of species counts or proportions.
#' @return The number of species observed, \eqn{R \geq 0}.
#' @details The richness is simply the number of nonzero elements in \code{x}.
#' Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{observed_otus()} in \code{skbio.diversity.alpha}.
#'   \item Equivalent to \code{specnumber} in \code{vegan}.
#'   \item Equivalent to the \code{sobs} calculator in Mothur.
#' }
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' richness(x) # 4
#' @export
richness <- function (x) {
  sum(x > 0)
}

#' Shannon diversity and related measures
#'
#' The Shannon index of diversity
#'
#' @param x A numeric vector of species counts or proportions.
#' @param base Base of the logarithm to use in the calculation.
#' @return The Shannon diversity, \eqn{H \geq 0}, or related quantity. The
#'   value of \eqn{H} is undefined if \code{x} sums to zero, and we return
#'   \code{NaN} in this case.  Heip's evenness measure and Pielou's Evenness
#'   index are undefined if only one element of \code{x} is nonzero, and again
#'   we return \code{NaN} if this is the case.
#' @details
#' The Shannon index of diversity or Shannon information entropy has deep roots
#' in information theory. It is defined as \deqn{H = - \sum_i p_i \log{p_i},}
#' where \eqn{p_i} is the species proportion. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{diversity()} in \code{vegan} with
#'     \code{index = "shannon"}.
#'   \item Equivalent to \code{shannon()} in \code{skbio.diversity.alpha}.
#' }
#'
#' The Brillouin index (Brillouin 1956) is similar to Shannon's index, but
#' accounts for sampling without replacement. For a vector of species counts
#' \code{x}, the Brillouin index is
#' \deqn{
#'   \frac{1}{N}\log{\frac{N!}{\prod_i x_i!}} =
#'   \frac{\log{N!} - \sum_i \log{x_i!}}{N}
#' } where \eqn{N} is the total number of counts. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{brillouin_d()} in \code{skbio.diversity.alpha}.
#'   \item Equivalent to the \code{shannon} calculator in Mothur.
#' }
#'
#' The Brillouin index accounts for the total number of individuals sampled,
#' and should be used on raw count data, not proportions.
#'
#' Heip's evenness measure is \deqn{\frac{e^H - 1}{S - 1},} where \eqn{S} is
#' the total number of species observed. Relation to other definitions:
#' \itemize{
#'   \item Equivalent to \code{heip_e()} in \code{skbio.diversity.alpha}.
#' }
#'
#' Pielou's Evenness index \eqn{J = H / \log{S}}. Relation to other
#' definitions:
#' \itemize{
#'   \item Equivalent to \code{peilou_e()} in \code{skbio.diversity.alpha}.
#' }
#'
#' @references
#' Brillouin L. Science and Information Theory. 1956;Academic Press, New York.
#'
#' Pielou EC. The Measurement of Diversity in Different Types of Biological
#' Collections. Journal of Theoretical Biology. 1966;13:131-144.
#' @examples
#' x <- c(15, 6, 4, 0, 3, 0)
#' shannon(x)
#'
#' # Using a different base is the same as dividing by the log of that base
#' shannon(x, base = 10)
#' shannon(x) / log(10)
#'
#' brillouin_d(x)
#'
#' # Brillouin index should be almost identical to Shannon index for large N
#' brillouin_d(10000 * x)
#' shannon(10000 * x)
#'
#' heip_e(x)
#' (exp(shannon(x)) - 1) / (richness(x) - 1)
#'
#' pielou_e(x)
#' shannon(x) / log(richness(x))
#' @export
shannon <- function (x, base=exp(1)) {
  p <- x / sum(x)
  # By convention, 0 * log(0) = 0
  p_is_defined_and_zero <- (p == 0) %in% TRUE
  p_logp <- ifelse(p_is_defined_and_zero, 0, p * log(p, base=base))
  -sum(p_logp)
}

#' @rdname shannon
#' @export
brillouin_d <- function (x) {
  n <- sum(x)
  nz <- x[x > 0]
  (lfactorial(n) - sum(lfactorial(nz))) / n
}

#' @rdname shannon
#' @export
heip_e <- function (x) {
  s <- sum(x > 0)
  h <- shannon(x)
  (exp(h) - 1) / (s - 1)
}

#' @rdname shannon
#' @export
pielou_e <- function (x) {
  h <- shannon(x)
  s <- sum(x > 0)
  h / log(s)
}

#' Strong's dominance index
#'
#' Strong's dominance index measures the maximum departure between the observed
#' proportions and a perfectly even community.
#' @param x A numeric vector of species counts.
#' @return Strong's dominance index, \eqn{0 \leq D_W < 1}. The index is
#'   undefined if \code{x} sums to 0, and we return \code{NaN} in this case.
#' @details
#' Strong's dominance index is defined as
#' \deqn{D_W = \max_i \left [ \frac{b_i}{N} - \frac{i}{S} \right ],} where
#' \eqn{b_i} is the abundance of the \eqn{i}th species, ordered from smallest
#' to largest, \eqn{N} is the total number of counts, and \eqn{S} is the number
#' of species observed.
#'
#' Equivalent to \code{strong()} in \code{skbio.diversity.alpha}.
#' @references
#' Strong WL. Assessing species abundance uneveness within and between plant
#' communities. Community Ecology. 2002;3:237-246.
#' @examples
#' x <- c(9, 0, 1, 2, 5, 2, 1, 1, 0, 7, 2, 1, 0, 1, 1)
#' strong(x)
#' @export
strong <- function (x) {
  n <- sum(x)
  s <- sum(x > 0)
  sorted_sum <- cumsum(sort(x, decreasing = TRUE))
  idx <- seq_along(x)
  max((sorted_sum / n) - (idx / s))
}

# Scikit-bio notes
# ace not implemented, TODO
# berger_parker_d implemented
# brillouin_d implemented
# chao1 not implemented, out of scope
# dominance implemented
# doubles not implemented, out of scope
# enspie implemented as invsimpson
# etsy_ci not implemented, out of scope
# faith_pd implemented
# fisher_alpha not implemented, TODO
# gini_index not implemented, TODO
# goods_coverage not implemented, out of scope
# heip_e implemented
# kempton_taylor_q implemented, but needs work
# lladser_ci not implemented, out of scope
# lladser_pe not implemented
# margalef implemented
# mcintosh_d implemented
# mcintosh_e implemented
# menhinick implemented
# michaelis_menten_fit not implemented
# observed_otus implemented as richness
# osd (obs. OTUs, singletons, doubletons) not implemented, out of scope
# pielou_e implemented
# robbins not implemented, out of scope
# shannon implemented
# simpson implemented
# simpson_e implemented
# singles not implemented, out of scope
# strong implemented

# Mothur notes
# ## Community richness
# sobs implemented as richness
# chao not implemented, out of scope
# ace not implemented, out of scope
# jack not implemented, out of scope
# bootstrap not implemented, out of scope
# ## Community evenness
# simpsoneven need to look at implementation, TODO
# shannoneven need to look at implementation, TODO
# heip need to look at implementation, TODO
# smithwilson need to look at implementation, TODO
# ## Community diversity
# bergerparker implemented as berger_parker_d
# shannon implemented
# npshannon not implemented, TODO
# simpson implemented (biased version) as dominance
# invsimpson implemented (biased version)
# coverage (aka Good's coverage) not implemented, out of scope
# qstat implemented as kempton_taylor_q, implementation does not match
# ## Estimates of number of additional OTUs
# boneh not implemented, out of scope
# efron not implemented, out of scope
# shen not implemented, out of scope
# solow not implemented, out of scope
