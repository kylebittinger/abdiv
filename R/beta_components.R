#' Nestedness and turnover components for presence/absence data
#' @param x,y Numeric vectors
#' @references
#' Baselga A. Partitioning the turnover and nestedness components of beta
#' diversity. Global Ecol. Biogeogr. 2010;19:134-143.
#'
#' Baselga A. The relationship between species replacement, dissimilarity
#' derived from nestedness, and nestedness. Global Ecol. Biogeogr.
#' 2012;21:1223–1232.
#' @export
jaccard_turnover_component <- function (x, y) {
    x <- x > 0
    y <- y > 0
    a <- sum(x & y)
    b = sum((!x) & y)
    c = sum(x & (!y))
    2 * min(b, c) / (a + 2 * min(b, c))
}

#' @rdname jaccard_turnover_component
#' @export
jaccard_nestedness_component <- function (x, y) {
  x <- x > 0
  y <- y > 0
  jaccard(x, y) - jaccard_turnover_component(x, y)
}

#' @rdname jaccard_turnover_component
#' @export
sorenson_turnover_component <- function (x, y) {
  x <- x > 0
  y <- y > 0
  a <- sum(x & y)
  b = sum((!x) & y)
  c = sum(x & (!y))
  min(b, c) / (a + min(b, c))
}

#' @rdname jaccard_turnover_component
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
#' @export
bray_curtis_balanced_component <- function (x, y) {
  # Check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  min(B, C) / (A + min(B, C))
}

#' @rdname bray_curtis_balanced_component
#' @export
bray_curtis_gradient_component <- function (x, y) {
  bray_curtis(x, y) - bray_curtis_balanced_component(x, y)
}

#' @rdname bray_curtis_balanced_component
#' @export
ruzicka_balanced_component <- function (x, y) {
  # Double check implementation between papers
  minxy <- pmin(x, y)
  A <- sum(minxy)
  B <- sum(x - minxy)
  C <- sum(y - minxy)
  2 * min(B, C) / (A +  2* min(B, C))
}

#' @rdname bray_curtis_balanced_component
#' @export
ruzicka_gradient_component <- function (x, y) {
  ruzicka(x, y) - ruzicka_balanced_component(x, y)
}

# TODO: Phylogenetic dissimilarity
