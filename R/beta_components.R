#' Nestedness and turnover components for presence/absence data
#' @param x,y Numeric vectors
#' @references
#' Baselga A. Partitioning the turnover and nestedness components of beta
#' diversity. Global Ecol. Biogeogr. 2010;19:134-143.
#'
#' Baselga A. The relationship between species replacement, dissimilarity
#' derived from nestedness, and nestedness. Global Ecol. Biogeogr.
#' 2012;21:1223–1232.
#' @name jaccard_components
NULL

#' @rdname jaccard_components
#' @export
jaccard_turnover_component <- function (x, y) {
    x <- x > 0
    y <- y > 0
    a <- sum(x & y)
    b = sum((!x) & y)
    c = sum(x & (!y))
    2 * min(b, c) / (a + 2 * min(b, c))
}

#' @rdname jaccard_components
#' @export
jaccard_nestedness_component <- function (x, y) {
  x <- x > 0
  y <- y > 0
  jaccard(x, y) - jaccard_turnover_component(x, y)
}

#' @rdname jaccard_components
#' @export
sorenson_turnover_component <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  b = sum((!x) & y)
  c = sum(x & (!y))
  min(b, c) / (a + min(b, c))
}

#' @rdname jaccard_components
#' @export
sorenson_nestedness_component <- function (x, y) {
  x <- x > 0
  y <- y > 0
  sorenson(x, y) - sorenson_turnover_component(x, y)
}

#' Balanced variation and abundance gradient components for abundance data
#' @param x,y Numeric vectors
#' @references
#' Baselga A. Separating the two components of abundance-based dissimilarity:
#' balanced changes in abundance vs.abundance gradients. Methods in Ecology
#' and Evolution. 2013;4:552–557.
#'
#' Baselga A. Partitioning abundance-based multiple-site dissimilarity into
#' components: balanced variation in abundance and abundance gradients.
#' Methods in Ecology and Evolution. 2017;8:799–808.
#' @name bray_curtis_components

#' @rdname bray_curtis_components
#' @export
bray_curtis_balanced_component <- function (x, y) {
  # Check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  min(B, C) / (A + min(B, C))
}

#' @rdname bray_curtis_components
#' @export
bray_curtis_gradient_component <- function (x, y) {
  bray_curtis(x, y) - bray_curtis_balanced_component(x, y)
}

#' @rdname bray_curtis_components
#' @export
ruzicka_balanced_component <- function (x, y) {
  # Double check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  2 * min(B, C) / (A +  2* min(B, C))
}

#' @rdname bray_curtis_components
#' @export
ruzicka_gradient_component <- function (x, y) {
  ruzicka(x, y) - ruzicka_balanced_component(x, y)
}

#' Nestedness and turnover components of unweighted UniFrac distance
#'
#' @param x,y Numeric vectors of species counts or proportions.
#' @param tree A phylogenetic tree object.
#' @param xy_labels A character vector of species labels for \code{x} and
#'   \code{y}.
#' @return The nestedness or turnover component of the UniFrac distance
#'   between communities \code{x} and \code{y}.
#' @details
#' Leprieur et al. (2012) showed that measures of phylogenetic beta diversity
#' could be partitioned into nestedness and turnover components, following the
#' approach of Baselga (2010) for Sorenson dissimilarity.
#' @references
#' Leprieur F, Albouy C, De Bortoli J, Cowman PF, Bellwood DR, Mouillot D.
#' Quantifying phylogenetic beta diversity: distinguishing between "true"
#' turnover of lineages and phylogenetic diversity gradients. PLoS One.
#' 2012;7(8):e42760. 10.1371/journal.pone.0042760
#' @name unifrac_components
NULL

#' @rdname unifrac_components
#' @export
unweighted_unifrac_turnover_component <- function (x, y, tree, xy_labels = NULL) {
  xy <- (x > 0) | (y > 0)
  pd_tot <- faith_pd(xy, tree, xy_labels)
  pd_x <- faith_pd(x, tree, xy_labels)
  pd_y <- faith_pd(y, tree, xy_labels)
  pd_min <- min(pd_tot - pd_x, pd_tot - pd_y)
  2 * pd_min / (pd_x + pd_y - pd_tot + 2 * pd_min)
}

#' @rdname unifrac_components
#' @export
unweighted_unifrac_nestedness_component <- function (x, y, tree, xy_labels = NULL) {
  unweighted_unifrac(x, y, tree, xy_labels) -
    unweighted_unifrac_turnover_component(x, y, tree, xy_labels)
}
