# To avoid dependence on ape library,
# generate code to make trees for testing
library(ape)
library(tidyverse)
library(ggtree)

code_for_tree <- function (varname, ...) {
  cat(paste0(varname, " <- "))
  tree_newick <- paste0(...)
  tree <- ape::read.tree(text=tree_newick)
  dput(tree)
}

faith_richards_newick <- "(((a:5,(b:2,c:1):3):4,(d:5,e:3):1):20);"
faith_richards_newick %>%
  ape::read.tree(text=.) %>%
  ape::rotateConstr(letters[5:1]) %>%
  ggtree::ggtree(ladderize = F) +
  geom_tiplab()

leprieur_newick <- "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"
leprieur_newick %>%
  ape::read.tree(text=.) %>%
  #ape::rotateConstr(letters[4:1]) %>%
  ggtree::ggtree(ladderize = F) +
  geom_tiplab()

code_for_tree(
  "skbio_t1",
  "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,",
  "(OTU4:0.75,OTU5:0.75):1.25):0.0)root;")

code_for_tree(
  "skbio_t1_w_extra_tips",
  "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,(OTU5:0.25,",
  "(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0)root;")

code_for_tree(
  "skbio_minimal",
  "(OTU1:0.25, OTU2:0.25)root;")

code_for_tree(
  "skbio_unobs_root",
  "((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)root;")

code_for_tree(
  "no_branch_lengths",
  "((((OTU1,OTU2),OTU3)),(OTU4,OTU5));")

.ag <- "((A:0.8,(B:0.6,(C:0.4,(D:0.3,E:0.3):0.1):0.2):0.2):0.3,(F:0.4,G:0.4):0.7)"
.hn <- "((H:0.8,(I:0.6,(J:0.4,(K:0.3,L:0.3):0.1):0.2):0.2):0.3,(M:0.4,N:0.4):0.7)"
unif_paper_newick <- paste0("((", .ag, ":1.1,", .hn, ":1.1):0.3)root;")
unif_paper_tree <- unif_paper_newick %>%
  ape::read.tree(text=.) %>%
  ape::rotateConstr(LETTERS[14:1])

unif_paper_cts_wide <- tibble(
  Tip = LETTERS[1:14],
  PanelA = c(1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2, 1),
  PanelB = c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2),
  PanelD = c(1, 1, 1, 3, 3, 1, 3, 1, 2, 2, 2, 2, 2, 3))

unif_paper_cts <- unif_paper_cts_wide %>%
  pivot_longer(
    starts_with("Panel"), names_to = "Panel",
    names_pattern = "Panel([ABD])", values_to = "Sample") %>%
  select(Tip, Panel, Sample) %>%
  mutate(Sample = factor(Sample)) %>%
  arrange(Panel, Tip) %>%
  mutate(Position = 2.6)

a_circle <- as.integer(unif_paper_cts_wide$PanelA == 1)
a_square <- as.integer(unif_paper_cts_wide$PanelA == 2)

b_circle <- as.integer(unif_paper_cts_wide$PanelB == 1)
b_square <- as.integer(unif_paper_cts_wide$PanelB == 2)

d_circle <- as.integer(unif_paper_cts_wide$PanelD == 1)
d_square <- as.integer(unif_paper_cts_wide$PanelD == 2)
d_triangle <- as.integer(unif_paper_cts_wide$PanelD == 3)

panelA <- unif_paper_cts %>%
  filter(Panel %in% "A") %>%
  mutate(
    Circle = as.integer(Sample %in% 1),
    Square = as.integer(Sample %in% 2)) %>%
  select(Tip, Circle, Square)

panelB <- unif_paper_cts %>%
  filter(Panel %in% "B") %>%
  mutate(
    Circle = as.integer(Sample %in% 1),
    Square = as.integer(Sample %in% 2)) %>%
  select(Tip, Circle, Square)

panelD <- unif_paper_cts %>%
  filter(Panel %in% "D") %>%
  mutate(
    Circle = as.integer(Sample %in% 1),
    Square = as.integer(Sample %in% 2),
    Triangle = as.integer(Sample %in% 3)) %>%
  select(Tip, Circle, Square, Triangle)

unif_paper_tree %>%
  ggtree::ggtree(ladderize = F) %<+%
  unif_paper_cts_wide +
  geom_tippoint(aes(shape=PanelA, x=Position), size=3) +
  scale_shape_manual(values = c(1, 15, 17))

unif_paper_tree %>%
  ggtree::ggtree(ladderize = F) %<+%
  filter(unif_paper_cts, Panel %in% "B") +
  geom_tippoint(aes(shape=Sample, x=Position), size=3) +
  scale_shape_manual(values = c(1, 15, 17))

unif_paper_tree %>%
  ggtree::ggtree(ladderize = F) %<+%
  filter(unif_paper_cts, Panel %in% "D") +
  geom_tippoint(aes(shape=Sample, x=Position), size=3) +
  scale_shape_manual(values = c(1, 15, 17))

