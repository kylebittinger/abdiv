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

#' Faith's phylogenetic diversity
#' @export
faith_pd <- function (x, tree) {
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
  px <- get_branch_abundances(em, x) > 0
  sum(b[px])
}

#' Unweighted UniFrac distance
#' @export
unweighted_unifrac <- function (x, y, tree) {
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
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

#' Weighted normalized UniFrac distance
#' @export
weighted_normalized_unifrac <- function (x, y, tree) {
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  sum(b * abs(px - py)) / sum(b * (px + py))
}

#' Weighted UniFrac distance
#' @export
weighted_unifrac <- function (x, y, tree) {
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  sum(b * abs(px - py))
}

#' Generalized UniFrac distance
#' @references Chen et al. Bioinformatics. 2012;28(16):2106-13.
#' @export
generalized_unifrac <- function (x, y, tree, alpha = 0.5) {
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
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

#' Variance-adjusted UniFrac distance
#' @references Chang et al. BMC Bioinformatics 2011, 12:118.
#' @export
variance_adjusted_unifrac <- function (x, y, tree, alpha = 0.5) {
  xy_sum <- x + y
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
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

#' Information UniFrac distance
#'
#' From Expanding the UniFrac Toolbox
#' @export
information_unifrac <- function (x, y, tree) {
  x <- x / sum(x)
  y <- y / sum(y)
  em <- make_edge_matrix(tree)
  b <- tree$edge.length
  if (is.null(b)) {
    stop("Tree has no branch lengths")
  }
  px <- get_branch_abundances(em, x)
  py <- get_branch_abundances(em, y)
  # log(0) produces infinities; p * log(p) must be zero for formula to work
  plogp <- function (p) ifelse(p > 0, p * log(p), 0)
  sum(b * abs(plogp(px) - plogp(py))) / sum(b)
}

