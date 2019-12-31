#' Diversity measures implemented
#'
#' The diversity functions offered in \code{abdiv} are organized based on the
#' function signature.
#'
#' @format
#' The following character vectors are provided:
#' \describe{
#'   \item{\code{alpha_diversities}}{
#'     All non-phylogenetic alpha diversity measures. These functions take a
#'     single numeric vector as an argument.}
#'   \item{\code{beta_diversities}}{
#'     All non-phylogenetic beta diversity measures. These functions take two
#'     numeric vectors as arguments.}
#'   \item{\code{phylogenetic_alpha_diversities}}{
#'     There is only one phylogenetic alpha diversity measure implemented, but
#'     we use the plural to be consistent with the other vectors. This function
#'     takes a numeric vector, a phylogenetic tree object, and optionally a
#'     character vector of species labels.}
#'   \item{\code{phylogenetic_beta_diversities}}{
#'     Phylogenetic measures of beta diversity. These functions take two
#'     numeric vectors, a phylogenetic tree object, and optionally a
#'     character vector of species labels.}
#' }
#' @name diversity_measures
NULL

#' @rdname diversity_measures
#' @export
alpha_diversities <- sort(c(
  "berger_parker_d", "brillouin_d", "dominance", "heip_e", "invsimpson",
  "kempton_taylor_q", "margalef", "mcintosh_d", "mcintosh_e", "menhinick",
  "pielou_e", "richness", "shannon", "simpson", "simpson_e", "strong"))

#' @rdname diversity_measures
#' @export
beta_diversities <- sort(c(
  # From beta.R
  "euclidean", "rms_distance", "chord", "hellinger", "geodesic_metric",
  "kullback_leibler_divergence", "manhattan", "mean_character_difference",
  "modified_mean_character_difference", "canberra",
  "clark_coefficient_of_divergence", "chebyshev", "correlation_distance",
  "cosine_distance", "bray_curtis", "weighted_kulczynski_second", "minkowski",
  "morisita", "horn_morisita", "binomial_deviance", "cy_dissimilarity",
  "ruzicka", "jaccard", "sorenson", "kulczynski_first", "kulczynski_second",
  "rogers_tanimoto", "russel_rao", "sokal_michener", "sokal_sneath",
  "yule_dissimilarity", "hamming",
  # From beta_components.R
  "jaccard_turnover_component", "jaccard_nestedness_component",
  "sorenson_turnover_component", "sorenson_nestedness_component",
  "bray_curtis_balanced_component", "bray_curtis_gradient_component",
  "ruzicka_balanced_component", "ruzicka_gradient_component"))

#' @rdname diversity_measures
#' @export
phylogenetic_alpha_diversities <- c("faith_pd")

#' @rdname diversity_measures
#' @export
phylogenetic_beta_diversities <- sort(c(
  # From unifrac.R
  "unweighted_unifrac", "weighted_unifrac", "weighted_normalized_unifrac",
  "variance_adjusted_unifrac", "generalized_unifrac", "information_unifrac",
  # From beta_components.R
  "unweighted_unifrac_turnover_component",
  "unweighted_unifrac_nestedness_component",
  "phylosor_turnover_component", "phylosor_nestedness_component"))
