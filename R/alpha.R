check_positive <- function (x) {
  if (any(x < 0)) {
    stop("Data in count vector is < 0")
  }
  TRUE
}

#' All alpha diversity measures
#'
#' We exclude functions that return multiple values.
#' @export
ecofuncs_alpha <- c(
    "berger_parker_d", "brillouin_d", "chao1", "dominance", "doubles",
    "enspie", "fisher_alpha", "goods_coverage", "heip_e", "invsimpson",
    "kempton_taylor_q", "margalef", "mcintosh_d", "mcintosh_e", "menhinick",
    "pielou_e", "richness", "robbins", "shannon", "simpson", "simpson_e",
    "singles", "strong")

# TODO: ace

#' Berger-Parker dominance
#' @export
berger_parker_d <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Berger-Parker dominance not defined for zero total counts")
    return(NA)
  }
  max(x) / sum(x)
}

#' Brillouin index
#' @export
brillouin_d <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Brillouin index not defined for zero total counts")
    return(NA)
  }
  n <- sum(x)
  nz <- x[x > 0]
  (lgamma(n + 1) - sum(lgamma(nz + 1))) / n
}

#' Chao1 index
#' @export
chao1 <- function (x, bias_corrected = TRUE) {
  s <- sum(x > 0)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)

  if ((!bias_corrected) & (f1 > 0) & (f2 > 0)) {
    s + (f1 ^ 2) / (f2 * 2)
  } else {
    s + f1 * (f1 - 1) / (2 * (f2 + 1))
  }
}

#' Dominance
#' @export
dominance <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Dominance not defined for zero total counts")
    return(NA)
  }
  p <- x / sum(x)
  sum(p ** 2)
}

#' Number of doubletons
#' @export
doubles <- function (x) {
  check_positive(x)
  sum(x == 2)
}

#' ENS diversity
#'
#' This is equivalent to inverse Simpson
#' @export
enspie <- function (x) {
  invsimpson(x)
}

#' Etsy's confidence interval
#'
#' The output has not been checked
#' @export
etsy_ci <- function (x, conf=0.975) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Etsy's CI not defined for zero total counts")
    return(NA)
  }
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  n <- sum(x)
  w <- (f1 * (n - f1) + 2 * n * f2) / (n ^ 3)
  z <- qnorm(conf)
  ci_center <- f1 / n
  ci_halfwidth <- z * sqrt(w)
  c(ci_center - ci_halfwidth, ci_center + ci_halfwidth)
}

# TODO: faith_pd

#' Fisher's alpha parameter
#' @export
fisher_alpha <- function (x) {
  check_positive(x)
  vegan::fisher.alpha(x)
}

# TODO: Gini-Simpson

#' Good's coverage of counts
#' @export
goods_coverage <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
   warning("Good's coverage not defined for zero total counts")
   return(NA)
  }
  f1 <- sum(x == 1)
  n <- sum(x)
  1 - (f1 / n)
}

#' Heip's evenness measure
#' @export
heip_e <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Heip's evenness not defined for zero total counts")
    return(NA)
  }
  s <- sum(x > 0)
  if (s == 1) {
    warning("Heip's evenness not defined for single species")
    return(NA)
  }
  h <- shannon(x)
  (exp(h) - 1) / (s - 1)
}

#' Inverse Simpson
#'
#' Defined to be 1 / D, like vegan.
#' @export
invsimpson <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Inverse simpson not defined for zero total counts")
    return(NA)
  }
  p <- x / sum(x)
  1 / sum(p ** 2)
}

#' Kempton-Taylor Q index
#' @export
kempton_taylor_q <- function (x, lower_quantile=0.25, upper_quantile=0.75) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Kempton-Taylor Q index not defined for zero total counts")
    return(NA)
  }
  # I'm sure there is a better way to do this with R's quantile function,
  # but not sure how to guarantee that the result always replicates the one
  # obtained via this algorithm.
  n <- length(x)
  lower_idx <- ceiling(n * lower_quantile) + 1
  upper_idx <- floor(n * upper_quantile) + 1
  x_sorted <- sort(x)
  x_upper <- x_sorted[upper_idx]
  x_lower <- x_sorted[lower_idx]
  (upper_idx - lower_idx) / log(x_upper / x_lower)
}

# TODO: lladser_pe

#' Margalef's richness index
#' @export
margalef <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Margalef's richness not defined for zero total counts")
    return(NA)
  }
  s <- sum(x > 0)
  n <- sum(x)
  (s - 1) / log(n)
}

#' McIntosh dominance index D
#' @export
mcintosh_d <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("McIntosh dominance not defined for zero total counts")
    return(NA)
  }
  n <- sum(x)
  u <- sqrt(sum(x ^ 2))
  (n - u) / (n - sqrt(n))
}

#' McIntosh's evenness measure E
#' @export
mcintosh_e <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("McIntosh evenness not defined for zero total counts")
    return(NA)
  }
  n <- sum(x)
  s <- sum(x > 0)
  numerator <- sqrt(sum(x ^ 2))
  denominator <- sqrt((n - s + 1) ^ 2 + s - 1)
  numerator / denominator
}

#' Menhinick's richness index
#' @export
menhinick <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Menhinick's richness not defined for zero total counts")
    return(NA)
  }
  n <- sum(x)
  s <- sum(x > 0)
  s / sqrt(n)
}

#TODO
michaelis_menten_fit <- function (x) {
  NA
}

#' Pielou's Evenness index J
#' @export
pielou_e <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Pielou's Evenness not defined for zero total counts")
    return(NA)
  }
  h <- shannon(x)
  s <- sum(x > 0)
  h / log(s)
}

#' Richness
#' @export
richness <- function (x) {
  check_positive(x)
  sum(x > 0)
}

#' Robbins' estimator for the probability of unobserved outcomes
#' @export
robbins <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Robbins' estimator not defined for zero total counts")
    return(NA)
  }
  f1 <- sum(x == 1)
  n <- sum(x)
  f1 / n
}

#' Shannon diversity
#' @export
shannon <- function (x, base=exp(1)) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Shannon diversity not defined zero total counts")
    return(NA)
  }
  # Remove zero values; they produce NaN's in the log function
  x <- x[x > 0]
  p <- x / sum(x)
  -sum(p * log(p, base=base))
}

#' Simpson's index
#'
#' Defined to be 1 - D, like vegan.
#' @export
simpson <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Simpson diversity not defined for zero total counts")
    return(NA)
  }
  p <- x / sum(x)
  1 - sum(p ** 2)
}

#' Simpson's evenness index
#' @export
simpson_e <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Simpson's evenness not defined for zero total counts")
    return(NA)
  }
  p <- x / sum(x)
  D <- sum(p ** 2)
  S <- sum(x > 0)
  1 / (D * S)
}

#' Number of singletons
#' @export
singles <- function (x) {
  check_positive(x)
  sum(x == 1)
}

#' Strong's dominance index
#' @export
strong <- function (x) {
  check_positive(x)
  if (sum(x) == 0) {
    warning("Strong's dominance not defined for zero total counts")
    return(NA)
  }
  n <- sum(x)
  s <- sum(x > 0)
  sorted_sum <- cumsum(sort(x, decreasing = TRUE))
  idx <- seq_along(x)
  max((sorted_sum / n) - (idx / s))
}
