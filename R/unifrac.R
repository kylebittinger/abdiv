# Matrix of branches (rows) vs. leaves (columns)
# Can be used to find the nodes to add up for each branch.
make_edge_matrix <- function (tree) {
  edge_nodes <- tree$edge[,2]
  np <- ape::nodepath(tree)
  sapply(np, function (x) edge_nodes %in% x)
}

get_branch_abundances <- function (edge_matrix, node_abundances) {
  as.numeric(edge_matrix %*% node_abundances)
}

check_tree <- function (tree) {
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths")
  }
}

#' Match vector of counts to phylogenetic tree
#'
#' @param x A vector of species counts.
#' @param tree A phylogenetic tree of class \code{"phylo"}.
#' @param x_labels A vector of species labels for \code{x}.
#' @return The vector \code{x}, re-arranged to match the tree.
#' @details This function applies a couple of different methods to arrange the
#'   data in \code{x} to match a phylogenetic tree.
#'   \itemize{
#'     \item If \code{x_labels} is provided, we use this vector to match the
#'       elements of \code{x} with the tip labels in the tree.
#'     \item If \code{x_labels} is not provided and \code{x} is a named vector,
#'       we use the names to match the tip labels in the tree.
#'     \item If \code{x_labels} is not provided and \code{x} is not named, we
#'       assume that \code{x} is already in the correct order, check that the
#'       length of \code{x} matches the number of tips in the tree, and return
#'       \code{x}.
#'   }
#' @export
match_to_tree <- function (x, tree, x_labels = NULL) {
  tree_labels <- tree$tip.label
  x_labels_are_provided <- !is.null(x_labels)
  if (x_labels_are_provided) {
    if (length(x) != length(x_labels)) {
      stop("Length of x does not match length of x_labels.")
    }
    if (!all(x_labels %in% tree_labels)) {
      x_labels_not_in_tree <- setdiff(x_labels, tree_labels)
      stop(paste("x_labels not found in tree:", x_labels_not_in_tree))
    }
    idx <- match(x_labels, tree_labels)
    xout <- numeric(length(tree_labels))
    xout[idx] <- x
    return(xout)
  }
  x_has_names <- !is.null(names(x))
  if (x_has_names) {
    if (!all(names(x) %in% tree_labels)) {
      x_names_not_in_tree <- setdiff(names(x), tree_labels)
      stop(paste("Names of x not found in tree:", x_names_not_in_tree))
    }
    idx <- match(names(x), tree_labels)
    xout <- numeric(length(tree_labels))
    xout[idx] <- x
    return(xout)
  }
  if (length(x) != length(tree_labels)) {
    stop("Length of x does not match number of tips in tree.")
  }
  return(x)
}

#' Faith's phylogenetic diversity
#'
#' Faith's phylogenetic diversity gives the total branch length on a
#' phylogenetic tree that is spanned by a community.  The abundance of each
#' species in the community is not considered.
#'
#' @param x A numeric vector of species counts or proportions, or a logical
#'   vector of species presence/absence.
#' @param tree A phylogenetic tree object..
#' @param x_labels A character vector of species labels for \code{x}.
#' @references Faith DP. Conservation evaluation and phylogenetic diversity.
#'   Biol. Conserv. 1992;61:1–10. doi: 10.1016/0006-3207(92)91201-3.
#' @export
faith_pd <- function (x, tree, x_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, x_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  sum(b[px])
}

#' UniFrac distance
#'
#' The UniFrac distance is a phylogenetically-weighted distance between two
#' communities of organisms. The measure has been extended a number of
#' times to include abundance-weighted and variance-adjusted versions.
#'
#' @param x,y Numeric vectors of species counts or proportions.
#' @param tree A phylogenetic tree object.
#' @param xy_labels A character vector of species labels for \code{x} and
#'   \code{y}.
#' @return The UniFrac distance between communities \code{x} and \code{y}.
#' @details
#' These functions compute different variations of the UniFrac distance between
#' communities described by the vectors \code{x} and \code{y}. If the vectors
#' are named, the names will be automatically used to match the vectors with
#' the tree. Missing names are filled in with zero counts. If the vectors are
#' not named and \code{xy_labels} is provided, these labels will be used to
#' match the vectors with the tree. If the vectors are not named and
#' \code{xy_labels} is not provided, it is assumed that the vectors are already
#' in the correct order, and we simply check that their length matches the
#' number of tips in the tree.
#'
#' \code{unweighted_unifrac} gives the original UniFrac distance from Lozupone
#' and Knight (2005), which is the fraction of total branch length leading to
#' community \code{x} or community \code{y}, but not both. It is based on
#' species presence/absence.
#'
#' \code{weighted_unifrac} gives the abundance-weighted version of UniFrac
#' proposed by Lozupone et al. (2007). In this measure, the branch lengths of
#' the tree are multiplied by the absolute difference in species abundances
#' below each branch.
#'
#' \code{weighted_normalized_unifrac} provides a normalized version of
#' \code{weighted_unifrac}, so the distance is between 0 and 1.
#'
#' \code{variance_adjusted_unifrac} was proposed by Chang et al. (2011) to
#' adjust for the variation of weights in weighted UniFrac under random
#' sampling.
#'
#' \code{generalized_unifrac} was proposed by Chen et al. (2012) to provide a
#' unifed matematical framework for weighted and unweighted UniFrac distance.
#' It includes a parameter, \eqn{\alpha}, which can be used to adjust the
#' abundance-weighting in the distance. A value of \eqn{\alpha = 1} corresponds
#' to weighted UniFrac. A value of \eqn{\alpha = 0} corresponds to unweighted
#' UniFrac if presence/absence vectors are provided.  The authors suggest a
#' value of \eqn{\alpha = 0.5} as a compromise between weighted and unweighted
#' distances.
#'
#' \code{information_unifrac} was proposed by Wong et al. (2016) to connect
#' UniFrac distance with compositional data analysis. They also proposed a
#' "ratio UniFrac" distance, which is not yet implemented.
#' @references
#' Lozupone C, Knight R. UniFrac: a new phylogenetic method for
#' comparing microbial communities. Applied and environmental microbiology.
#' 2005;71(12):8228–8235. 10.1128/AEM.71.12.8228-8235.2005
#'
#' Lozupone CA, Hamady M, Kelley ST, Knight R. Quantitative and
#' qualitative \eqn{\beta} diversity measures lead to different insights into
#' factors that structure microbial communities. Applied and environmental
#' microbiology. 2007;73(5):1576–1585. 10.1128/AEM.01996-06
#'
#' Chang Q., et al. Variance adjusted weighted UniFrac: a powerful
#' beta diversity measure for comparing communities based on phylogeny.
#' BMC Bioinformatics. 2011;12:118. 10.1186/1471-2105-12-118
#'
#' Chen J, Bittinger K, Charlson ES, Hoffmann C, Lewis J, Wu GD,
#' et al. Associating microbiome composition with environmental covariates
#' using generalized UniFrac distances. Bioinformatics.
#' 2012;28(16):2106–2113. 10.1093/bioinformatics/bts342
#'
#' Wong RG, Wu JR, Gloor GB. Expanding the UniFrac Toolbox.
#' PLOS ONE. 2016;11(9):1–20. 10.1371/journal.pone.0161196
#' @name unifrac
NULL

#' @rdname unifrac
#' @export
unweighted_unifrac <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x) > 0
  py <- get_branch_abundances(em, y) > 0
  # Unique branch length is the sum of branch lengths where
  # either px or py is nonzero, but not both
  branch_unique <- sum(b * xor(px, py))
  # Total branch length is the sum of branch lengths where
  # px, py, or both are nonzero
  branch_total <- sum(b[px | py])
  branch_unique / branch_total
}

#' @rdname unifrac
#' @export
weighted_unifrac <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  sum(b * abs(px - py))
}

#' @rdname unifrac
#' @export
weighted_normalized_unifrac <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  sum(b * abs(px - py)) / sum(b * (px + py))
}

#' @rdname unifrac
#' @export
variance_adjusted_unifrac <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  xy_sum <- x + y
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  m <- get_branch_abundances(em, xy_sum)
  # (m_total - m) appears in the denominator; must remove double zeroes now
  m_total <- sum(xy_sum)
  keep <- (m != m_total) & (m > 0)
  b <- b[keep]
  px <- px[keep]
  py <- py[keep]
  m <- m[keep]
  # Typos in paper: this quantity should be square rooted, and the first m
  # should have an i subscript
  m_denom <- sqrt(m * (m_total - m))
  sum(b * abs(px - py) / m_denom) / sum(b * (px + py) / m_denom)
}

#' @rdname unifrac
#' @export
generalized_unifrac <- function (x, y, tree, alpha = 0.5, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  # px + py appears in the denominator; must remove double zeroes now
  keep <- (px + py) > 0
  b <- b[keep]
  px <- px[keep]
  py <- py[keep]
  p_sum <- px + py
  sum(b * (p_sum ^ (alpha - 1)) * abs(px - py)) / sum(b * (p_sum ^ alpha))
}

#' @rdname unifrac
#' @export
information_unifrac <- function (x, y, tree, xy_labels = NULL) {
  check_tree(tree)
  x <- match_to_tree(x, tree, xy_labels)
  y <- match_to_tree(y, tree, xy_labels)
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  # log(0) produces infinities; p * log(p) must be zero for formula to work
  plogp <- function (p) ifelse(p > 0, p * log(p), 0)
  sum(b * abs(plogp(px) - plogp(py))) / sum(b)
}
