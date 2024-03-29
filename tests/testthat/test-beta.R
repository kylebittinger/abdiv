context("beta diversity")

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
  expect_equal(
    canberra(x1, x2), from_vegdist(x1, x2, "canberra") * sum(x1 | x2))
  expect_equal(bray_curtis(x1, x2), from_vegdist(x1, x2, "bray"))
  expect_equal(
    weighted_kulczynski_second(x1, x2), from_vegdist(x1, x2, "kulczynski"))
  expect_equal(
    kulczynski_second(x1, x2),
    from_vegdist(x1, x2, "kulczynski", binary = TRUE))
  expect_equal(ruzicka(x1, x2), from_vegdist(x1, x2, "jaccard"))
  expect_equal(jaccard(x1, x2), from_vegdist(x1, x2, "jaccard", binary = TRUE))
  # gower not implemented
  expect_equal(
    modified_mean_character_difference(x1, x2),
    from_vegdist(x1, x2, "altGower"))
  expect_equal(morisita(x1, x2), from_vegdist(x1, x2, "morisita"))
  expect_equal(horn_morisita(x1, x2), from_vegdist(x1, x2, "horn"))
  expect_equal(binomial_deviance(x1, x2), from_vegdist(x1, x2, "binomial"))
  expect_equal(
    cy_dissimilarity(x1, x2, base=exp(1)), from_vegdist(x1, x2, "cao"))
  # We don't implement the "mahalanobis" method
  # abundance_jaccard is not equal to vegdist(x, method = "chao")
})

test_that("Distance functions are consistent with Legendre & Legendre", {
  # Example under 7.19
  q1 <- c(9, 3, 7, 3, 4, 9, 5, 4, 1, 6)
  q2 <- c(2, 3, 2, 1, 2, 9, 3, 2, 1, 6)
  expect_equal(1 - hamming(q1, q2) / length(q1), 0.4)
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
  # Under 7.37
  # x1 and x3, x2 and x3 at an angle of 90 degrees
  expect_equal(geodesic_metric(x1, x3), pi / 2)
  expect_equal(geodesic_metric(x2, x3), pi / 2)
  # x1 and x2 at an angle of 18.4 degrees
  expect_equal(geodesic_metric(x1, x2), 18.4349488229 * pi / 180)
  # Hellinger is chord on sqrt transformed proportions
  expect_equal(
    hellinger(x1, x2),
    chord(sqrt(x1 / sum(x1)), sqrt(x2 / sum(x2))))
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
    Inf)
})

test_that("Distance functions are consistent with scipy.spatial.distance", {
  # Using the iris dataset, results taken from scipy
  iris1 <- c(5.1, 3.5, 1.4, 0.2)
  iris2 <- c(4.9, 3.0, 1.4, 0.2)
  expect_equal(euclidean(iris1, iris2), 0.53851648)
  expect_equal(chebyshev(iris1, iris2), 0.5)
  expect_equal(manhattan(iris1, iris2), 0.7)
  expect_equal(correlation_distance(iris1, iris2), 0.0040013388)
  expect_equal(cosine_distance(iris1, iris2), 0.0014208365)
  expect_equal(minkowski(iris1, iris2, 3.2), 0.50817745)
  expect_equal(minkowski(iris1, iris2, 5.8), 0.50042326)

  # Using scipy's double-inp file
  d1 <- c(
    0.827893804941075, 0.903529398447625, 0.186218899467949, 0.892115131231046,
    0.0206185911937958, 0.344063672738573, 0.153377991283033, 0.57013723000098,
    0.551002073021156, 0.1792362258426, 0.808617512087658, 0.611548718431718,
    0.0123347178716485, 0.00144164353187104, 0.404430920904569,
    0.356139895949991, 0.128198571292975, 0.866330083384748, 0.869602778629158,
    0.361172737036377, 0.528353765877262, 0.144024108809012, 0.311245722713895,
    0.603128079689789, 0.923032479274252, 0.233212188113687,
    0.0319265226740344, 0.346620629499556, 0.298868772804637,
    0.0511674954204809, 0.258497583091449, 0.430202347804223,
    0.800397275171352, 0.93649319113681, 0.973709864996467, 0.471803845397223,
    0.452659168660786, 0.10564856785208, 0.588301971428541, 0.384609223767698,
    0.646150005343547, 0.101323972984882, 0.121615156165119, 0.515966892948466,
    0.845207447351023, 0.988517096224797, 0.762388307349013,
    0.0229116324361543, 0.577553098080238, 0.782069989682809,
    0.823918634584297, 0.339180010526023, 0.954631845161454, 0.37896779178677,
    0.0452653339964929, 0.836678647323859, 0.308263681104986,
    0.117393682079345, 0.0763199496916944, 0.299741665072218,
    0.579520865516023, 0.394235089254201, 0.117512638329726, 0.492823251395003,
    0.942129399622595, 0.0836539105384134, 0.686805969357184,
    0.358952796242944, 0.759293942716606, 0.562384946613145, 0.211074682803205,
    0.98246837046686, 0.266123014224624, 0.616227231500712, 0.50232545366075,
    0.0520285447666978, 0.58350906688421, 0.786464211888914, 0.250401238686751,
    0.672830864113599, 0.46107935345761, 0.482050877051591, 0.972040325102227,
    0.31000692852635, 0.768101712646175, 0.0795653930600708, 0.259338963788774,
    0.113785259040305, 0.388530307328445, 0.859909466007596,
    0.0521516787591828, 0.162090824857229, 0.185923609045766,
    0.624771651261048, 0.341512849552078, 0.703490336837803, 0.603756464001957,
    0.233896943442331, 0.010021048856099, 0.786605840396904)
  d2 <- c(
    0.803369411603336, 0.865326454554403, 0.746834041075404, 0.63624309199106,
    0.0512000630662547, 0.950334837263359, 0.469773260962682, 0.422130528845943,
    0.315345211983839, 0.299101484344266, 0.119066796728026, 0.348656771450934,
    0.828949364988505, 0.845481105080001, 0.91496730182119, 0.77087078371939,
    0.264015773212255, 0.210789702218961, 0.420763305505444, 0.67195002846547,
    0.145803168489306, 0.0180041273588613, 0.0840273343522001,
    0.0420676015688316, 0.137693351504131, 0.171671734102213,
    0.178822072765216, 0.822431043340212, 0.772909366686748, 0.206422362102598,
    0.959209203622721, 0.8312490243755, 0.66732893603699, 0.0463284790369077,
    0.764395409835898, 0.93593415256151, 0.191496631916303, 0.453659046940287,
    0.864083601653801, 0.0394152917817546, 0.560210199520548,
    0.926380616194166, 0.155599532594482, 0.617220810295012, 0.63355767528121,
    0.0976697546036804, 0.0447579568953987, 0.3248842796105, 0.57003771221495,
    0.906696296725681, 0.545846062150568, 0.683340128558149, 0.288724440954404,
    0.131633864701683, 0.232567330524599, 0.415012196318841, 0.383484546636606,
    0.814936577396873, 0.18670038494502, 0.317032217354302, 0.683209366268268,
    0.172972851892911, 0.923655735970264, 0.915294125215009, 0.722487998309662,
    0.855792062659806, 0.534488305925164, 0.487687327444911, 0.830827780450642,
    0.391662448932221, 0.345969512227397, 0.403351249902741, 0.655572644491301,
    0.713845240938024, 0.168393731459997, 0.176938214348644, 0.758868365517814,
    0.375058989288082, 0.752517624512621, 0.60839611525383, 0.114597230990799,
    0.623961448580955, 0.13076554820659, 0.853045875067092, 0.480160207012477,
    0.0816812218986355, 0.379313962274464, 0.149698699777684, 0.71290238783029,
    0.683097923743805, 0.763537594387651, 0.182400496325123, 0.576469584899234,
    0.88651132487316, 0.5784337085544, 0.970002662875512, 0.731820734790506,
    0.385140139393671, 0.17742918511934, 0.97634232292423)
  expect_equal(chebyshev(d1, d2), 0.89084734)
  expect_equal(manhattan(d1, d2), 32.420590)
  expect_equal(correlation_distance(d1, d2), 0.92507465)
  expect_equal(cosine_distance(d1, d2), 0.25695885)
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
  expect_equal(hamming(b1, b2), 0.46 * 100)

  # Other tests
  # mahalanobis() not tested
  # matching() is deprecated
  expect_equal(jaccard(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 0.6)
  expect_equal(jaccard(c(1, 0, 1), c(1, 1, 0)), 2 / 3)
  expect_equal(yule_dissimilarity(c(1, 0, 1, 1, 0), c(1, 1, 0, 1, 1)), 2)
  expect_equal(yule_dissimilarity(c(1, 0, 1), c(1, 1, 0)), 2)
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
  expect_equal(
    minkowski(c(1, 2, 3), c(1, 1, 5), 1.5), (2 ** 1.5 + 1) ** (2 / 3))
  expect_equal(minkowski(c(1, 2, 3), c(1, 1, 5), 2), sqrt(5))
  expect_equal(euclidean(c(1, 2, 3), c(1, 1, 5)), sqrt(5))
  expect_equal(
    cosine_distance(c(1, 2, 3), c(1, 1, 5)), 1 - 18 / (sqrt(14) * sqrt(27)))
  expect_equal(correlation_distance(c(1, 2, 3), c(1, 1, 5)), 0.1339746)
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
  expect_equal(cosine_distance(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(cosine_distance(c(100, 0, 0), c(0, 1, 0)), 1)
  expect_equal(
    cosine_distance(c(1, 0, 0), c(1, 1, 0)),  1 - 1 / (3 * sqrt(2 / 9)))
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
  expect_equal(hamming(c(1, 0, 0), c(0, 1, 0)), 2) # 2 / 3 in scipy
  expect_equal(hamming(c(1, 0, 0), c(1, 1, 0)), 1) # 1 / 3 in scipy
  expect_equal(hamming(c(1, 0, 0), c(2, 0, 0)), 1) # 1 / 3 in scipy
  expect_equal(hamming(c(1, 0, 0), c(3, 0, 0)), 1) # 1 / 3 in scipy
  expect_equal(jaccard(c(1, 0, 0), c(0, 1, 0)), 1)
  expect_equal(jaccard(c(1, 0, 0), c(1, 1, 0)), 0.5)
  expect_equal(jaccard(c(1, 0, 0), c(1, 2, 0)), 0.5)
  expect_equal(jaccard(c(1, 0, 0), c(1, 1, 1)), 2 / 3)
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
  expect_equal(yule_dissimilarity(c(1, 0, 0), c(0, 1, 0)), 2)
  expect_equal(yule_dissimilarity(c(1, 1, 0), c(0, 1, 0)), 0)
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
  # https://www.mothur.org/wiki/Anderberg
  expect_equal(sokal_sneath(forest, pasture), 1 - 9 / (2 * 33 + 2 * 31 - 3 * 9))
  # https://www.mothur.org/wiki/Hamming
  expect_equal(hamming(forest > 0, pasture > 0), 33 + 31 - 2 * 9)
  # https://www.mothur.org/wiki/Jclass
  expect_equal(jaccard(forest, pasture), 1 - 9 / (33 + 31 - 9))
  # jest not implemented
  # https://www.mothur.org/wiki/Kulczynski
  expect_equal(kulczynski_first(forest, pasture), 1 - 9 / (33 + 31 - 2 * 9))
  # https://www.mothur.org/wiki/Kulczynskicody
  expect_equal(kulczynski_second(forest, pasture), 1 - 0.5 * (9 / 33 + 9 / 31))
  # lennon not implemented
  # memchi2 not implemented
  # memchord has no example
  # https://www.mothur.org/wiki/Memeuclidean
  expect_equal(euclidean(forest > 0, pasture > 0), 6.7823299831)
  # https://www.mothur.org/wiki/Mempearson
  expect_equal(correlation_distance(forest > 0, pasture > 0), 1.7184212081)
  # ochiai not implemented
  # https://www.mothur.org/wiki/Sorclass
  expect_equal(sorenson(forest, pasture), 1 - 2 * 9 / (33 + 31))
  # sorest not implemented
  # https://www.mothur.org/wiki/Whittaker
  expect_equal(sorenson(forest, pasture), 1 - (2 - 2 * 55 / (33 + 31)))
  # https://www.mothur.org/wiki/Braycurtis
  expect_equal(bray_curtis(forest, pasture), 1 - 2 * 11 / (49 + 49))
  # canberra has no example
  # gower not implemented
  # hellinger has no example
  # https://www.mothur.org/wiki/jabund
  # Mothur implements the estimators for the distance, not the observed
  # distance. For the observed distance, u = 16/49 and v = 18/49
  #
  expect_equal(
    abundance_jaccard(forest, pasture),
    1 - (16 / 49) * (18 / 49) / ((16 / 49) + (18 / 49) - (16 / 49) * (18 / 49)))
  # manhattan has no example
  # https://www.mothur.org/wiki/Morisitahorn
  expect_equal(
    horn_morisita(forest, pasture),
    1 - 2 * (33 / (49 * 49)) / (99 / (49 ^ 2) + 131 / (49 ^ 2)))
  # odum has no example
  # soergel has no example
  # https://www.mothur.org/wiki/sorabund
  # See note for jabund
  expect_equal(
    abundance_sorenson(forest, pasture),
    1 - 2 * (16 / 49) * (18 / 49) / ((16 / 49) + (18 / 49)))
  # spearman has no example
  # speciesprofile has no example
  # structchi2 has no example
  # structchord has no example
  # structeuclidean has no example
  # structkulczynski has no example
  # structpearson has no example
  # thetan not implemented
  # thetayc not implemented
})

test_that("Cao index matches Cao (1997) paper", {
  # Columns of Table 7
  expect_equal(cy_dissimilarity(1, 0), 0.6494535986)
  expect_equal(cy_dissimilarity(10, 0), 1.6834893979)
  expect_equal(cy_dissimilarity(100, 90), 0.0018064951)
  expect_equal(cy_dissimilarity(10, 2), 0.3606262540)
})

test_that("Clark matches Clark 1952 paper", {
  Ai <- c(8.000, 10.059, 4.000, 4.206, 3.941)
  Bi <- c(7.969, 10.078, 1.016, 3.078, 2.250)
  Ci <- c(8.005,  9.968, 1.032, 2.853, 2.739)
  # 0.301 in paper
  expect_equal(clark_coefficient_of_divergence(Ai, Bi), 0.3008296906)
  # 0.047 in paper
  expect_equal(clark_coefficient_of_divergence(Bi, Ci), 0.0472068860)
  # 0.289 in paper
  expect_equal(clark_coefficient_of_divergence(Ai, Ci), 0.2888010841)
})

test_that("Mean character difference matches Cain 1958 paper", {
  # Table 1
  A <- c(100, 100,  75, 5,  90, 100, 90)
  B <- c( 90,  85, 100, 5, 100,  50, 95)
  # Table 2
  expect_equal(mean_character_difference(A, B), 115 / 7)
})

test_that("Empty vectors give correct results", {
  mt <- c(0, 0, 0, 0)
  expect_equal(euclidean(mt, mt), 0)
  expect_equal(rms_distance(mt, mt), 0)
  expect_identical(chord(mt, mt), NaN)
  expect_identical(hellinger(mt, mt), NaN)
  expect_identical(geodesic_metric(mt, mt), NaN)
  expect_identical(kullback_leibler_divergence(mt, mt), NaN)
  expect_equal(manhattan(mt, mt), 0)
  expect_equal(mean_character_difference(mt, mt), 0)
  expect_identical(modified_mean_character_difference(mt, mt), NaN)
  expect_equal(canberra(mt, mt), 0)
  expect_identical(clark_coefficient_of_divergence(mt, mt), NaN)
  expect_equal(chebyshev(mt, mt), 0)
  expect_identical(correlation_distance(mt, mt), NaN)
  expect_identical(cosine_distance(mt, mt), NaN)
  expect_identical(bray_curtis(mt, mt), NaN)
  expect_identical(weighted_kulczynski_second(mt, mt), NaN)
  expect_equal(minkowski(mt, mt), 0)
  expect_identical(morisita(mt, mt), NaN)
  expect_identical(horn_morisita(mt, mt), NaN)
  expect_equal(binomial_deviance(mt, mt), 0)
  expect_identical(cy_dissimilarity(mt, mt), NaN)
  expect_identical(ruzicka(mt, mt), NaN)
  expect_identical(jaccard(mt, mt), NaN)
  expect_identical(sorenson(mt, mt), NaN)
  expect_identical(kulczynski_first(mt, mt), NaN)
  expect_identical(kulczynski_second(mt, mt), NaN)
  expect_equal(rogers_tanimoto(mt, mt), 0)
  expect_equal(russel_rao(mt, mt), 1)
  expect_equal(sokal_michener(mt, mt), 0)
  expect_identical(sokal_sneath(mt, mt), NaN)
  expect_identical(yule_dissimilarity(mt, mt), NaN)
  expect_equal(hamming(mt, mt), 0)
  expect_identical(abundance_jaccard(mt, mt), NaN)
  expect_identical(abundance_sorenson(mt, mt), NaN)
})

test_that("Single observation gives correct results", {
  x <- c(0, 0, 0, 0)
  y <- c(0, 0, 0, 1)
  expect_equal(euclidean(x, y), 1)
  expect_equal(rms_distance(x, y), 0.5) # sqrt((0 + 0 + 0 + 1) / 4)
  expect_identical(chord(x, y), NaN)
  expect_identical(hellinger(x, y), NaN)
  expect_identical(geodesic_metric(x, y), NaN)
  expect_identical(kullback_leibler_divergence(x, y), NaN)
  expect_equal(manhattan(x, y), 1)
  expect_equal(mean_character_difference(x, y), 0.25)
  expect_identical(modified_mean_character_difference(x, y), 1)
  expect_equal(canberra(x, y), 1)
  expect_equal(clark_coefficient_of_divergence(x, y), 1)
  expect_equal(chebyshev(x, y), 1)
  expect_identical(correlation_distance(x, y), NaN)
  expect_identical(cosine_distance(x, y), NaN)
  expect_equal(bray_curtis(x, y), 1)
  expect_identical(weighted_kulczynski_second(x, y), NaN)
  expect_equal(minkowski(x, y), 1)
  expect_identical(morisita(x, y), NaN)
  expect_identical(horn_morisita(x, y), NaN)
  expect_equal(binomial_deviance(x, y), log(2)) # (0 + 0 + 1 * log(2)) / 1
  expect_equal(cy_dissimilarity(x, y), 0.6494535986)
  expect_equal(ruzicka(x, y), 1)
  expect_equal(jaccard(x, y), 1)
  expect_equal(sorenson(x, y), 1)
  expect_equal(kulczynski_first(x, y), 1)
  expect_identical(kulczynski_second(x, y), NaN)
  expect_equal(rogers_tanimoto(x, y), 2 / 5)
  expect_equal(russel_rao(x, y), 1)
  expect_equal(sokal_michener(x, y), 2 / 5)
  expect_equal(sokal_sneath(x, y), 1)
  expect_identical(yule_dissimilarity(x, y), NaN)
  expect_equal(hamming(x, y), 1)
  expect_equal(abundance_jaccard(x, y), NaN) # should it be 1, like Jaccard?
  expect_equal(abundance_sorenson(x, y), NaN) # should it be 1, like Sorenson?
})

test_that("Minkowski distance does not allow invalid p", {
  expect_error(minkowski(c(1,2,3), c(4,5,6), p = 0))
  expect_error(minkowski(c(1,2,3), c(4,5,6), p = c(1,2,3)))
})
