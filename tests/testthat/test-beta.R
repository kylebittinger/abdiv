context("beta diversity")

test_that("Distances are computed correctly for data frame", {
  df <- data.frame(
    Observation = rep(c("A", "B", "C"), each = 5),
    Feature = paste0("F", c(1:5, 2:6, 4:8)),
    Val = 1:15,
    stringsAsFactors = FALSE)
  df_wide <- pivot_to_numeric_matrix(df, Observation, Feature, Val)
  df_wide_vals <- c(
    1, 2, 3,  4,  5,  0,  0,  0,
    0, 6, 7,  8,  9, 10,  0,  0,
    0, 0, 0, 11, 12, 13, 14, 15)
  expect_equal(df_wide, matrix(
    df_wide_vals, byrow = TRUE, nrow=3,
    dimnames = list(c("A", "B", "C"), paste0("F", 1:8))))
  expected_manhattan_dist <- structure(
    c(27, 62, 51), Size = 3L, Labels = c("A", "B", "C"),
    Diag = FALSE, Upper = FALSE, class = "dist")
  observed_manhattan_dist <- dist_long(
    df, Observation, Feature, Val, manhattan)
  expect_equal(observed_manhattan_dist, expected_manhattan_dist)
})

test_that("Distance functions are consistent with stats::dist()", {
  x1 <- c(0, 1, 2, 3, 4)
  x2 <- c(3, 0, 2, 1, 5)
  from_dist <- function (x, y, method, ...) {
    as.numeric(dist(rbind(x1, x2), method = method, ...))
  }
  expect_equal(euclidean(x1, x2), from_dist(x1, x2, "euclidean"))
  # We dont implement the "maximum" method
  expect_equal(manhattan(x1, x2), from_dist(x1, x2, "manhattan"))
  expect_equal(canberra(x1, x2), from_dist(x1, x2, "canberra"))
  expect_equal(jaccard(x1, x2), from_dist(x1, x2, "binary"))
  expect_equal(minkowski(x1, x2, 3), from_dist(x1, x2, "minkowski", p=3))
})

test_that("Distance functions are consistent with vegan::vegdist()", {
  x1 <- c(0, 0, 0, 1, 2, 3, 5)
  x2 <- c(0, 0, 3, 0, 2, 1, 1)
  x3 <- 1:5
  x4 <- 6:10
  from_vegdist <- function (x, y, method, ...) {
    as.numeric(vegan::vegdist(rbind(x, y), method = method, ...))
  }
  expect_equal(manhattan(x1, x2), from_vegdist(x1, x2, "manhattan"))
  expect_equal(euclidean(x1, x2), from_vegdist(x1, x2, "euclidean"))
  expect_equal(canberra(x1, x2), from_vegdist(x1, x2, "canberra") * sum(x1 | x2))
  expect_equal(bray_curtis(x1, x2), from_vegdist(x1, x2, "bray"))
  expect_equal(kulczynski(x1, x2), from_vegdist(x1, x2, "kulczynski"))
  expect_equal(jaccard(x1, x2), from_vegdist(x1, x2, "jaccard", binary=TRUE))
  expect_equal(gower(x1, x2), from_vegdist(x1, x2, "gower"))
  expect_equal(alt_gower(x1, x2), from_vegdist(x1, x2, "altGower"))
  expect_equal(morisita(x1, x2), from_vegdist(x1, x2, "morisita"))
  expect_equal(morisita_horn(x1, x2), from_vegdist(x1, x2, "horn"))
  expect_equal(millar(x1, x2), from_vegdist(x1, x2, "binomial"))
  expect_equal(cao(x1, x2), from_vegdist(x1, x2, "cao"))
  # We don't implement the "mahalanobis" method
})

test_that("Distance functions are consistent with Legendre & Legendre", {
  # Example under 7.19
  q1 <- c(9, 3, 7, 3, 4, 9, 5, 4, 1, 6)
  q2 <- c(2, 3, 2, 1, 2, 9, 3, 2, 1, 6)
  expect_equal(1 - hamming(q1, q2), 0.4)
  # Example under 7.33
  x1 <- c(0, 4, 8)
  x2 <- c(0, 1, 1)
  x3 <- c(1, 0, 0)
  expect_equal(euclidean(x1, x2), sqrt(sum(c(0, 3, 7) ** 2)))
  expect_equal(euclidean(x1, x3), 9)
  expect_equal(euclidean(x2, x3), sqrt(3))
  # Under 7.34
  r1 <- c(1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0)
  r2 <- c(0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0)
  expect_equal(rms_distance(r1, r2), sqrt(5 / 12))
  expect_equal(rms_distance(r1[1:8], r2[1:8]), sqrt(5 / 8))
  # Under 7.36
  expect_equal(chord(x1, x2), 0.32036449)
  expect_equal(chord(x1, x3), sqrt(2))
  expect_equal(chord(x2, x3), sqrt(2))

  expect_equal(euclidean(c(0, 4, 8), c(1, 0, 0)), 9)
  expect_equal(euclidean(c(0, 1, 1), c(1, 0, 0)), sqrt(3))

  # Under 7.57
  expect_equal(sorenson(c(1, 1, 1, 0, 0), c(0, 0, 0, 1, 1)), 1)
  expect_equal(sorenson(c(1, 1, 1, 0, 0), c(1, 1, 1, 1, 1)), 0.25)
  expect_equal(sorenson(c(0, 0, 0, 1, 1), c(1, 1, 1, 1, 1)), 3 / 7)
  # Under 7.58
  expect_equal(bray_curtis(c(2, 5, 2, 5, 3), c(3, 5, 2, 4, 3)), 2 / 34)
  expect_equal(bray_curtis(c(2, 5, 2, 5, 3), c(9, 1, 1, 1, 1)), 18 / 30)
  expect_equal(bray_curtis(c(3, 5, 2, 4, 3), c(9, 1, 1, 1, 1)), 16 / 30)

})

test_that("Kullback-Leibler divergence works", {
  px <- c(0.36, 0.48, 0.16)
  qx <- c(1/3, 1/3, 1/3)
  # Example from wikipedia
  expect_equal(kullback_leibler_divergence(px, qx), 0.0852996)
  # KL divergence is not symmetric
  expect_equal(kullback_leibler_divergence(qx, px), 0.097455)

  # Vectors should be scaled to probability distributions
  expect_equal(
    kullback_leibler_divergence(5 * px, 2 * qx),
    kullback_leibler_divergence(px, qx))

  # Whenever P(x) is zero the contribution of the corresponding term is
  # interpreted as zero.
  expect_equal(
    kullback_leibler_divergence(c(0.5, 0.5, 0), c(0.25, 0.5, 0.25)),
    0.3465736)

  # The Kullback–Leibler divergence is defined only if
  # for all Q(x)=0 implies P(x)=0 (absolute continuity).
  expect_equal(
    kullback_leibler_divergence(c(0, 0, 0.5, 0.5), c(0.25, 0.75, 0, 0)),
    NaN)
})

test_that("Distance functions are consistent with scipy.spatial.distance", {
  # Using the iris dataset, results taken from scipy
  iris1 <- c(5.1, 3.5, 1.4, 0.2)
  iris2 <- c(4.9, 3.0, 1.4, 0.2)
  expect_equal(euclidean(iris1, iris2), 0.53851648)
  expect_equal(chebyshev(iris1, iris2), 0.5)
  expect_equal(manhattan(iris1, iris2), 0.7)
  expect_equal(correlation(iris1, iris2), 0.0040013388)
  expect_equal(cosine(iris1, iris2), 0.0014208365)
  expect_equal(minkowski(iris1, iris2, 3.2), 0.50817745)
  expect_equal(minkowski(iris1, iris2, 5.8), 0.50042326)

  # Using scipy's double-inp file
  d1 <- c(
    0.827893804941075, 0.903529398447625, 0.186218899467949, 0.892115131231046,
    0.0206185911937958, 0.344063672738573, 0.153377991283033, 0.57013723000098,
    0.551002073021156, 0.1792362258426, 0.808617512087658, 0.611548718431718,
    0.0123347178716485, 0.00144164353187104, 0.404430920904569, 0.356139895949991,
    0.128198571292975, 0.866330083384748, 0.869602778629158, 0.361172737036377,
    0.528353765877262, 0.144024108809012, 0.311245722713895, 0.603128079689789,
    0.923032479274252, 0.233212188113687, 0.0319265226740344, 0.346620629499556,
    0.298868772804637, 0.0511674954204809, 0.258497583091449, 0.430202347804223,
    0.800397275171352, 0.93649319113681, 0.973709864996467, 0.471803845397223,
    0.452659168660786, 0.10564856785208, 0.588301971428541, 0.384609223767698,
    0.646150005343547, 0.101323972984882, 0.121615156165119, 0.515966892948466,
    0.845207447351023, 0.988517096224797, 0.762388307349013, 0.0229116324361543,
    0.577553098080238, 0.782069989682809, 0.823918634584297, 0.339180010526023,
    0.954631845161454, 0.37896779178677, 0.0452653339964929, 0.836678647323859,
    0.308263681104986, 0.117393682079345, 0.0763199496916944, 0.299741665072218,
    0.579520865516023, 0.394235089254201, 0.117512638329726, 0.492823251395003,
    0.942129399622595, 0.0836539105384134, 0.686805969357184, 0.358952796242944,
    0.759293942716606, 0.562384946613145, 0.211074682803205, 0.98246837046686,
    0.266123014224624, 0.616227231500712, 0.50232545366075, 0.0520285447666978,
    0.58350906688421, 0.786464211888914, 0.250401238686751, 0.672830864113599,
    0.46107935345761, 0.482050877051591, 0.972040325102227, 0.31000692852635,
    0.768101712646175, 0.0795653930600708, 0.259338963788774, 0.113785259040305,
    0.388530307328445, 0.859909466007596, 0.0521516787591828, 0.162090824857229,
    0.185923609045766, 0.624771651261048, 0.341512849552078, 0.703490336837803,
    0.603756464001957, 0.233896943442331, 0.010021048856099, 0.786605840396904)
  d2 <- c(
    0.803369411603336, 0.865326454554403, 0.746834041075404, 0.63624309199106,
    0.0512000630662547, 0.950334837263359, 0.469773260962682, 0.422130528845943,
    0.315345211983839, 0.299101484344266, 0.119066796728026, 0.348656771450934,
    0.828949364988505, 0.845481105080001, 0.91496730182119, 0.77087078371939,
    0.264015773212255, 0.210789702218961, 0.420763305505444, 0.67195002846547,
    0.145803168489306, 0.0180041273588613, 0.0840273343522001, 0.0420676015688316,
    0.137693351504131, 0.171671734102213, 0.178822072765216, 0.822431043340212,
    0.772909366686748, 0.206422362102598, 0.959209203622721, 0.8312490243755,
    0.66732893603699, 0.0463284790369077, 0.764395409835898, 0.93593415256151,
    0.191496631916303, 0.453659046940287, 0.864083601653801, 0.0394152917817546,
    0.560210199520548, 0.926380616194166, 0.155599532594482, 0.617220810295012,
    0.63355767528121, 0.0976697546036804, 0.0447579568953987, 0.3248842796105,
    0.57003771221495, 0.906696296725681, 0.545846062150568, 0.683340128558149,
    0.288724440954404, 0.131633864701683, 0.232567330524599, 0.415012196318841,
    0.383484546636606, 0.814936577396873, 0.18670038494502, 0.317032217354302,
    0.683209366268268, 0.172972851892911, 0.923655735970264, 0.915294125215009,
    0.722487998309662, 0.855792062659806, 0.534488305925164, 0.487687327444911,
    0.830827780450642, 0.391662448932221, 0.345969512227397, 0.403351249902741,
    0.655572644491301, 0.713845240938024, 0.168393731459997, 0.176938214348644,
    0.758868365517814, 0.375058989288082, 0.752517624512621, 0.60839611525383,
    0.114597230990799, 0.623961448580955, 0.13076554820659, 0.853045875067092,
    0.480160207012477, 0.0816812218986355, 0.379313962274464, 0.149698699777684,
    0.71290238783029, 0.683097923743805, 0.763537594387651, 0.182400496325123,
    0.576469584899234, 0.88651132487316, 0.5784337085544, 0.970002662875512,
    0.731820734790506, 0.385140139393671, 0.17742918511934, 0.97634232292423)
  expect_equal(chebyshev(d1, d2), 0.89084734)
  expect_equal(manhattan(d1, d2), 32.420590)
  expect_equal(correlation(d1, d2), 0.92507465)
  expect_equal(cosine(d1, d2), 0.25695885)
  expect_equal(euclidean(d1, d2), 4.0515260)
  expect_equal(minkowski(d1, d2, 3.2), 2.021505044)

  # Using scipy's boolean-inp file
  b1 <- c(
    1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0,
    1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1,
    0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0,
    0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1,
    0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1)
  b2 <- c(
    1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1,
    1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0,
    1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1)
  expect_equal(jaccard(b1, b2), 6.5714286e-01)
  expect_equal(hamming(b1, b2), 0.46)

  # Other tests
  # mahalanobis() not tested
  # matching() is deprecated
  expect_equal(jaccard(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 0.6)
  expect_equal(jaccard(c(1, 0, 1), c(1, 1, 0)), 2 / 3)
  expect_equal(yule(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 2)
  expect_equal(yule(c(1, 0, 1), c(1, 1, 0)), 2)
  expect_equal(sorenson(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 3 / 7)
  expect_equal(sorenson(c(1, 0, 1), c(1, 1, 0)), 0.5)
  expect_equal(sokal_sneath(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 0.75)
  expect_equal(sokal_sneath(c(1, 0, 1), c(1, 1, 0)), 0.8)
  expect_equal(rogers_tanimoto(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 0.75)
  expect_equal(rogers_tanimoto(c(1, 0, 1), c(1, 1, 0)), 0.8)
  expect_equal(russel_rao(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 0.6)
  expect_equal(russel_rao(c(1, 0, 1), c(1, 1, 0)), 2 / 3)
  expect_equal(canberra(3.3, 3.4), 0.01492537)
  expect_equal(minkowski(c(1, 2, 3), c(1, 1, 5), 1), 3)
  expect_equal(minkowski(c(1, 2, 3), c(1, 1, 5), 1.5), (2 ** 1.5 + 1) ** (2 / 3))
  expect_equal(minkowski(c(1, 2, 3), c(1, 1, 5), 2), sqrt(5))
  expect_equal(euclidean(c(1, 2, 3), c(1, 1, 5)), sqrt(5))
  expect_equal(cosine(c(1, 2, 3), c(1, 1, 5)), 1 - 18 / (sqrt(14) * sqrt(27)))
  expect_equal(correlation(c(1, 2, 3), c(1, 1, 5)), 0.1339746)
  expect_equal(canberra(c(1, 2, 3), c(2, 4, 6)), 1)
  expect_equal(canberra(c(1, 1, 0, 0), c(1, 0, 1, 0)), 2)
  expect_equal(bray_curtis(c(1, 2, 3), c(2, 4, 6)), 1 / 3)
  expect_equal(bray_curtis(c(1, 1, 0, 0), c(1, 0, 1, 0)), 0.5)
  expect_equal(euclidean(c(1, 1, 1), c(0, 0, 0)), sqrt(3))

  # From doctests
  expect_equal(bray_curtis(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(bray_curtis(c(1, 0, 0), c(1, 1, 0)), 1 / 3)
  expect_equal(canberra(c(1, 0, 0), c(0, 1, 0)), 2)
  expect_equal(canberra(c(1, 0, 0), c(1, 1, 0)), 1)
  expect_equal(chebyshev(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(chebyshev(c(1, 0, 0), c(1, 1, 0)), 1)
  # In scipy, Manhattan distance is called cityblock
  expect_equal(manhattan(c(1, 0, 0), c(0, 1, 0)), 2)
  expect_equal(manhattan(c(1, 0, 0), c(0, 2, 0)), 3)
  expect_equal(manhattan(c(1, 0, 0), c(1, 1, 0)), 1)
  # No doctests for correlation distance
  expect_equal(cosine(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(cosine(c(100, 0, 0), c(0, 1, 0)), 1)
  expect_equal(cosine(c(1, 0, 0), c(1, 1, 0)),  1 - 1 / (3 * sqrt(2 / 9)))
  expect_equal(euclidean(c(1, 0, 0), c(0, 1, 0)), sqrt(2))
  expect_equal(euclidean(c(1, 0, 0), c(1, 1, 0)), 1)
  # mahalanobis is not implemented
  expect_equal(minkowski(c(1, 0, 0), c(0, 1, 0), 1), 2)
  expect_equal(minkowski(c(1, 0, 0), c(0, 1, 0), 2), sqrt(2))
  expect_equal(minkowski(c(1, 0, 0), c(0, 1, 0), 3), 2 ** (1 / 3))
  expect_equal(minkowski(c(1, 0, 0), c(1, 1, 0), 1), 1)
  expect_equal(minkowski(c(1, 0, 0), c(1, 1, 0), 2), 1)
  expect_equal(minkowski(c(1, 0, 0), c(1, 1, 0), 3), 1)
  # seeuclidean is not implemented
  # sqeuclidean is trivial, checks against tests above for euclidean distance
  # wminkowski is deprecated in scipy
  # In scipy, Sorenson distance is called dice
  expect_equal(sorenson(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 3 / 7)
  expect_equal(sorenson(c(1, 0, 1), c(1, 1, 0)), 0.5)
  expect_equal(hamming(c(1, 0, 0), c(0, 1, 0)), 2 / 3)
  expect_equal(hamming(c(1, 0, 0), c(1, 1, 0)), 1 / 3)
  expect_equal(hamming(c(1, 0, 0), c(2, 0, 0)), 1 / 3)
  expect_equal(hamming(c(1, 0, 0), c(3, 0, 0)), 1 / 3)
  expect_equal(jaccard(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(jaccard(c(1, 0, 0), c(1, 1, 0)), 0.5)
  expect_equal(jaccard(c(1, 0, 0), c(1, 2, 0)), 0.5)
  expect_equal(jaccard(c(1, 0, 0), c(1, 1, 1)), 2 / 3)
  expect_equal(kulczynski_scipy(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(kulczynski_scipy(c(1, 0, 0), c(1, 1, 0)), 0.75)
  # Scipy says this is 1 / 3 based on their weird transformation
  # With a reasonable transformation, the answer should be 0.75, same as above
  expect_equal(kulczynski_scipy(c(1, 0, 0), c(2, 1, 0)), 0.75)
  # Same here. Scipy gets -0.5, should be 0.75, as above
  expect_equal(kulczynski_scipy(c(1, 0, 0), c(3, 1, 0)), 0.75)
  # Example from https://github.com/scipy/scipy/issues/2009
  expect_equal(kulczynski_scipy(c(1, 1, 0, 0), c(0, 1, 1, 0)), 5 / 6)
  expect_equal(rogers_tanimoto(c(1, 0, 0), c(0, 1, 0)), 0.8)
  expect_equal(rogers_tanimoto(c(1, 0, 0), c(1, 1, 0)), 0.5)
  # Scipy gets a value of -1 with their implementation
  # Correct value is 0; vectors have equivalent species presence/absence
  expect_equal(rogers_tanimoto(c(1, 0, 0), c(2, 0, 0)), 0)
  expect_equal(russel_rao(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(russel_rao(c(1, 0, 0), c(1, 1, 0)), 2 / 3)
  # Scipy transformation gives 1 / 3 here
  # But more resonable interpretation is that 2 of the 3 values are not
  # present in both samples, therefore the distance should be 2 / 3
  expect_equal(russel_rao(c(1, 0, 0), c(2, 0, 0)), 2 / 3)
  # I think Sokal-Michener is EXACTLY the same as Rogers-Tanimoto
  # The formula is definitely the same, and I think the implementation
  # in scipy is also the same.
  expect_equal(sokal_michener(c(1, 0, 0), c(0, 1, 0)), 0.8)
  expect_equal(sokal_michener(c(1, 0, 0), c(1, 1, 0)), 0.5)
  # Scipy gets a value of -1 with their implementation
  # Correct value is 0; vectors are equivalent in presence/absence
  expect_equal(sokal_michener(c(1, 0, 0), c(2, 0, 0)), 0)
  expect_equal(sokal_sneath(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(sokal_sneath(c(1, 0, 0), c(1, 1, 0)), 2 / 3)
  # Scipy implementation gives 0, should be 2 / 3, as above
  expect_equal(sokal_sneath(c(1, 0, 0), c(2, 1, 0)), 2 / 3)
  # Scipy implementation gives -2, should be 2 / 3, as above
  expect_equal(sokal_sneath(c(1, 0, 0), c(3, 1, 0)), 2 / 3)
  expect_equal(yule(c(1, 0, 0), c(0, 1, 0)), 2)
  expect_equal(yule(c(1, 1, 0), c(0, 1, 0)), 0)
})

test_that("Distance functions are consistent with Mothur", {
  forest <- c(
    1, 1, 1, 1, 1, 1, 3, 3, 2, 2, 1, 1, 3, 2, 1, 1, 1, 1, 2, 1,
    1, 2, 5, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
  pasture <- c(
    0, 0, 0, 1, 0, 0, 1, 0, 0, 5, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3,
    0, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 1, 1, 1, 1, 1, 7, 1,
    1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1)
  # Mothur has a different definition of Kulczynski dissimilarity
  # https://www.mothur.org/wiki/Kulczynski
  expect_equal(kulczynski_mothur(forest, pasture), 1 - 9 / (33 + 31 - 2 * 9))
  # Mothur calls the vegan's kulczynski method Kulczynski-Cody
  # https://www.mothur.org/wiki/Kulczynskicody
  expect_equal(kulczynski_cody(forest, pasture), 1 - 0.5 * (9 / 33 + 9 / 31))
})
