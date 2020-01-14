#' Nestedness and turnover components for presence/absence data
#' @param x,y Numeric vectors
#' @return The nestedness or turnover component of distance between \code{x}
#'   and \code{y}. This quantity is undefined when either \code{x} or \code{y}
#'   have no observations, in which case we return \code{NaN}.
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
jaccard_turnover <- function (x, y) {
    x <- x > 0
    y <- y > 0
    a <- sum(x & y)
    b <- sum((!x) & y)
    c <- sum(x & (!y))
    2 * min(b, c) / (a + 2 * min(b, c))
}

#' @rdname jaccard_components
#' @export
jaccard_nestedness <- function (x, y) {
  x <- x > 0
  y <- y > 0
  jaccard(x, y) - jaccard_turnover(x, y)
}

#' @rdname jaccard_components
#' @export
sorenson_turnover <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  b <- sum((!x) & y)
  c <- sum(x & (!y))
  min(b, c) / (a + min(b, c))
}

#' @rdname jaccard_components
#' @export
sorenson_nestedness <- function (x, y) {
  x <- x > 0
  y <- y > 0
  sorenson(x, y) - sorenson_turnover(x, y)
}

#' Balanced variation and abundance gradient components for abundance data
#' @param x,y Numeric vectors
#' @return The balanced variation or abundance gradient component of distance
#'   between \code{x} and \code{y}. This quantity is undefined when either
#'   \code{x} or \code{y} have all elements equal to zero, in which case we
#'   return \code{NaN}.
#' @references
#' Baselga A. Separating the two components of abundance-based dissimilarity:
#' balanced changes in abundance vs.abundance gradients. Methods in Ecology
#' and Evolution. 2013;4:552–557.
#'
#' Baselga A. Partitioning abundance-based multiple-site dissimilarity into
#' components: balanced variation in abundance and abundance gradients.
#' Methods in Ecology and Evolution. 2017;8:799–808.
#' @name bray_curtis_components
NULL

#' @rdname bray_curtis_components
#' @export
bray_curtis_balanced <- function (x, y) {
  # Check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  min(B, C) / (A + min(B, C))
}

#' @rdname bray_curtis_components
#' @export
bray_curtis_gradient <- function (x, y) {
  bray_curtis(x, y) - bray_curtis_balanced(x, y)
}

#' @rdname bray_curtis_components
#' @export
ruzicka_balanced <- function (x, y) {
  # Double check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  2 * min(B, C) / (A +  2* min(B, C))
}

#' @rdname bray_curtis_components
#' @export
ruzicka_gradient <- function (x, y) {
  ruzicka(x, y) - ruzicka_balanced(x, y)
}

#' Nestedness and turnover components of unweighted UniFrac distance
#'
#' @param x,y Numeric vectors of species counts or proportions.
#' @param tree A phylogenetic tree object.
#' @param xy_labels A character vector of species labels for \code{x} and
#'   \code{y}.
#' @return The nestedness or turnover component of the UniFrac distance
#'   between communities \code{x} and \code{y}. This quantity is undefined
#'   when either \code{x} or \code{y} have all elements equal to zero, in
#'   which case we return \code{NaN}.
#' @details
#' Leprieur et al. (2012) showed that measures of phylogenetic beta diversity
#' could be partitioned into nestedness and turnover components, following the
#' approach of Baselga (2010) for Sorenson dissimilarity.
#' @references
#' Baselga A. Partitioning the turnover and nestedness components of beta
#' diversity. Global Ecol. Biogeogr. 2010;19:134-143.
#'
#' Leprieur F, Albouy C, De Bortoli J, Cowman PF, Bellwood DR, Mouillot D.
#' Quantifying phylogenetic beta diversity: distinguishing between "true"
#' turnover of lineages and phylogenetic diversity gradients. PLoS One.
#' 2012;7(8):e42760. 10.1371/journal.pone.0042760
#' @examples
#' # Vectors x and y have turnover but no nestedness
#' x <- c(1, 1, 1, 0, 0, 0, 0, 0)
#' y <- c(0, 1, 1, 1, 0, 0, 0, 0)
#'
#' unweighted_unifrac(x, y, leprieur_tree)
#' unweighted_unifrac_turnover(x, y, leprieur_tree)
#' unweighted_unifrac_nestedness(x, y, leprieur_tree)
#'
#' phylosor(x, y, leprieur_tree)
#' phylosor_turnover(x, y, leprieur_tree)
#' phylosor_nestedness(x, y, leprieur_tree)
#'
#' # Vectors y and z have nestedness but no turnover
#' z <- c(0, 1, 1, 1, 1, 1, 1, 1)
#'
#' unweighted_unifrac(y, z, leprieur_tree)
#' unweighted_unifrac_turnover(y, z, leprieur_tree)
#' unweighted_unifrac_nestedness(y, z, leprieur_tree)
#'
#' phylosor(y, z, leprieur_tree)
#' phylosor_turnover(y, z, leprieur_tree)
#' phylosor_nestedness(y, z, leprieur_tree)
#' @name unifrac_components
NULL

#' @rdname unifrac_components
#' @export
unweighted_unifrac_turnover <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  py <- get_branch_abundances(em, y) > 0
  # Re-write formula 17 with a, b, and c from formulas 12-14.
  .a <- sum(b[px & py])
  .b <- sum(b[px & (!py)])
  .c <- sum(b[(!px) & py])
  2 * min(.b, .c) / (.a + 2 * min(.b, .c))
}

#' @rdname unifrac_components
#' @export
unweighted_unifrac_nestedness <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  py <- get_branch_abundances(em, y) > 0
  .a <- sum(b[px & py])
  .b <- sum(b[px & (!py)])
  .c <- sum(b[(!px) & py])
  unifrac_total <- (.b + .c) / (.a + .b + .c)
  unifrac_turnover <- 2 * min(.b, .c) / (.a + 2 * min(.b, .c))
  unifrac_total - unifrac_turnover
}

#' @rdname unifrac_components
#' @export
phylosor_turnover <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  py <- get_branch_abundances(em, y) > 0
  .a <- sum(b[px & py])
  .b <- sum(b[px & (!py)])
  .c <- sum(b[(!px) & py])
  min(.b, .c) / (.a + min(.b, .c))
}

#' @rdname unifrac_components
#' @export
phylosor_nestedness <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  py <- get_branch_abundances(em, y) > 0
  .a <- sum(b[px & py])
  .b <- sum(b[px & (!py)])
  .c <- sum(b[(!px) & py])
  phylosor_total <- (.b + .c) / (2 * .a + .b + .c)
  phylosor_turnover <- min(.b, .c) / (.a + min(.b, .c))
  phylosor_total - phylosor_turnover
}

#' Example data for phylogenetic nestedness and turnover components
#'
#' This tree was used in Figure 2 of Leprieur et al. (2005) to demonstrate
#' the nestedness and turnover components of phylogenetic beta diversity.
#'
#' @usage
#' leprieur_tree
#'
#' @format \code{leprieur_tree} is a phylogenetic tree with 8 tips, labeled
#' a-h. It was created with the \code{ape} library. All edges (branches) in the
#' tree are of length 1.
#'
#' @source
#' Leprieur F, Albouy C, De Bortoli J, Cowman PF, Bellwood DR, Mouillot D.
#' Quantifying phylogenetic beta diversity: distinguishing between "true"
#' turnover of lineages and phylogenetic diversity gradients. PLoS One.
#' 2012;7(8):e42760. 10.1371/journal.pone.0042760
"leprieur_tree"
