
<!-- README.md is generated from README.Rmd. Please edit that file -->

# abdiv

This package re-implements measures of ecological diversity from several
other software packages, including `vegan`, `scikit-bio`, and
`GUniFrac`.

## Installation

You can install abdiv from github with:

``` r
# install.packages("devtools")
devtools::install_github("kylebittinger/abdiv")
```

## Alpha diversity

Let’s say we’ve surveyed a field and counted the number of plants in
each of two sites. We’ve found five species in total, and we’d like to
summarize the diversity of the two sampling sites. The diversity within
each site is called α-diversity.

Here are the number of plants for each species at site 1 and site 2,
represented as vectors.

``` r
site1 <- c(2, 5, 16, 0, 1)
site2 <- c(0, 0, 8, 8, 8)
```

The two sites have about the same number of total plants, but the
ditsribution of species is much different. More than half the plants in
site 1 belong to a single species, whereas the plants in site 2 are
almost evenly distributed across three different species. To get
started, let’s look at a few ways to quantify the α-diversity for each
sample.

Richness measures the total number of species in each sample.

``` r
richness(site1)
## [1] 4
richness(site2)
## [1] 3
```

The Shannon index measures both the number of species and the evenness
of the relative abundance values. Site 2 has fewer species, but each
species has the same relative abundance, so the Shannon index is higher.

``` r
shannon(site1)
## [1] 0.9365995
shannon(site2)
## [1] 1.098612
```

Let’s summarize the diversity of site 1 and site 2 using all the
functions available in this library. The full set of within-sample
diversity functions is available as a character vector in
`alpha_diversities`. The term “α-diversity” means the diversity within a
single sample.

``` r
library(tidyverse)
tibble(Measure = alpha_diversities) %>%
  group_by(Measure) %>%
  summarize(Site1 = get(Measure)(site1), Site2 = get(Measure)(site2)) %>%
  pivot_longer(-Measure, "SampleID") %>%
  ggplot(aes(x=Measure, y=value, color=SampleID)) +
  geom_point() +
  scale_color_manual(values=c("#E64B35", "#4DBBD5")) +
  coord_flip() +
  theme_bw()
```

![](tools/readme/alpha-diversity-1.png)<!-- -->

We can see that site 1 is regarded as more diverse by some measures; it
has the most species. For other measures, site 2 is regarded as more
diverse; it has the most even distribution of species.

In our documentation, you can find more info on each α-diversity
function.

## Beta diversity

Having assessed the diversity within each sample, we can next ask about
the number of shared species between sites. If species are shared, how
similar is the distribution across species? There are many ways to
quantify the between-sample diversity or β-diversity.

You can think about β-diversity as either the similarity or
dissimilarity between sites. The functions in `abdiv` are written in
terms of dissimilarity: similar sites will have values close to zero,
and highly dissimilar sites will have values close to the maximum. This
way of thinking goes along with our intuition about diversity: sites
with greater dissimilarity will exhibit increased diversity if we
consider both sites together.

Let’s look at some examples. The Jaccard distance counts the fraction of
species present in only one site. The answer is 3 out of 5, or 0.6.

``` r
jaccard(site1, site2)
## [1] 0.6
```

The Bray-Curtis dissimilarity adds up the absolute differences between
species counts, then divides by the total counts. For our two sites,
that’s (2 + 5 + 8 + 8 + 7) / 48, or 0.625.

``` r
bray_curtis(site1, site2)
## [1] 0.625
```

Again, we’ll use a vector called `beta_diversities` to compute every
dissimilarity measure in the library.

``` r
tibble(Measure = beta_diversities) %>%
  group_by(Measure) %>%
  mutate(value = get(Measure)(site1, site2)) %>%
  ggplot(aes(x=Measure, y=value)) +
  geom_point(color="#4DBBD5") +
  scale_y_log10() +
  coord_flip() +
  theme_bw()
```

![](tools/readme/beta-diversity-1.png)<!-- -->

The dissimilarities are generally positive, and they have a range of
scales. Some dissimilarity measures range from 0 to 1, while others can
go up indefinitely.

As before, you can find more info on each β-diversity function in our
documentation.

## Phylogenetic diversity and UniFrac

A phylogenetic tree can be incorporated into diversity measures to add
information about the evolutionary history of species. One of the first
measures in this area was Faith’s phylogenetic diversity. Here, we take
an example from Faith and Richards (PMID 24832524) to show how
phylogenetic diversity works.

Here is the tree from the paper:

``` r
library(ggtree)
faith_tree <- ape::read.tree(text="(((a:5,(b:2,c:1):4):4,(d:5,e:3):1):20);")
faith_tree <- ape::rotateConstr(faith_tree, letters[5:1])
ggtree(faith_tree, ladderize = F) +
  geom_tiplab()
```

![](tools/readme/unnamed-chunk-7-1.png)<!-- -->

If all the species are present, the value of Faith’s phylogenetic
diversity (PD) is the sum of the branch lengths. Here, we expect the
total branch length to be 5 + 4 + 2 + 4 + 1 + 20 + 5 + 3 + 1 = 45,
adding from top to bottom.\[1\]

Here is how you would calculate Faith’s PD:

``` r
faith_pd(c(1, 1, 1, 1, 1), faith_tree)
## [1] 45
```

If species “a” is missing, we expect Faith’s PD to be reduced by 5, the
length of the branch leading to species “a”.

``` r
faith_pd(c(0, 1, 1, 1, 1), faith_tree)
## [1] 40
```

The practice of using phylogenetic information in diversity has been
especially popular in microbial ecology, where bacteria are suveyed
using the 16S rRNA marker gene. In addition to serving as a fingerprint
for bacteria, the gene sequence can be used to build a phylogenetic
tree.

In this area, the UniFrac distance is widely used to measure β-diversity
between bacterial communities. We’ll reproduce an example from the
UniFrac paper by Lozupone and Knight (PMID 16332807). Here is panel A of
figure 1:

``` r
unifrac_tree <- ape::read.tree(text = paste0(
  "((",
  "((A:0.8,(B:0.6,(C:0.4,(D:0.3,E:0.3):0.1):0.2):0.2):0.3,(F:0.4,G:0.4):0.7)",
  ":1.1,",
  "((H:0.8,(I:0.6,(J:0.4,(K:0.3,L:0.3):0.1):0.2):0.2):0.3,(M:0.4,N:0.4):0.7)",
  ":1.1):0.3)root;"))
unifrac_tree <- ape::rotateConstr(unifrac_tree, LETTERS[14:1])
unifrac_paper_data <- tibble(
  Species = LETTERS[1:14],
  PanelA = factor(c(1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2, 1)),
  PanelB = factor(c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2)))
ggtree(unifrac_tree, ladderize = F) %<+%
  unifrac_paper_data +
  geom_tippoint(aes(shape=PanelA), x=2.6, size=3) +
  scale_shape_manual(values = c(1, 15))
```

![](tools/readme/unnamed-chunk-10-1.png)<!-- -->

In this example, we’ve measured two bacterial communities, where the
species detected in each are labeled with circles and squares. They have
no species in common, but the species from each community spans the
entire tree. The communities are not very phylogenetically distinct. In
fact, if we measure the branch length leading uniquely to either a
circle or a square, we sum up the branch lengths on the far right-hand
side of the plot. This is about half the total branch length for the
tree.

``` r
a_circle <- c(1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1)
a_square <- c(0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0)
unweighted_unifrac(a_circle, a_square, unifrac_tree)
## [1] 0.5378151
```

We can increase the UniFrac distance if we re-arrange the species in
each community so that all the circles are in the upper part of the tree
and all the squares are in the lower part.

``` r
ggtree(unifrac_tree, ladderize = F) %<+%
  unifrac_paper_data +
  geom_tippoint(aes(shape=PanelB), x=2.6, size=3) +
  scale_shape_manual(values = c(1, 15))
```

![](tools/readme/unnamed-chunk-12-1.png)<!-- -->

With this arrangement, the species in each community are
phylogenetically very different. Except for the root, each branch of the
tree leads uniquely to a species present in square or circle, and never
to a species present in both. Therefore, the UniFrac distance is close
to 1.

``` r
b_circle <- c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
b_square <- c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
unweighted_unifrac(b_circle, b_square, unifrac_tree)
## [1] 0.9747899
```

1.  In the paper, they give the value as 41, but they don’t assign a
    length to the edge connecting species “b” and “c”. Looking at the
    figure, I’ve estimated that the length should be 4, so we get 45 for
    our example, rather than 41.
