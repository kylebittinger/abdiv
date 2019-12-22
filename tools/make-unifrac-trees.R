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

code_for_tree(
  "no_branch_lengths",
  "((((OTU1,OTU2),OTU3)),(OTU4,OTU5));")

code_for_tree(
  "skbio_t1",
  "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,",
  "(OTU4:0.75,OTU5:0.75):1.25):0.0)root;")

code_for_tree(
  "skbio_minimal",
  "(OTU1:0.25, OTU2:0.25)root;")

code_for_tree(
  "skbio_unobs_root",
  "((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)root;")

code_for_tree(
  "skbio_t1_w_extra_tips",
  "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,(OTU5:0.25,",
  "(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0)root;")




