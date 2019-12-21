# Phylogenetic tree from Figure 1 of Lozupone and Knight (2005)

.ag <- "((A:0.8,(B:0.6,(C:0.4,(D:0.3,E:0.3):0.1):0.2):0.2):0.3,(F:0.4,G:0.4):0.7)"
.hn <- "((H:0.8,(I:0.6,(J:0.4,(K:0.3,L:0.3):0.1):0.2):0.2):0.3,(M:0.4,N:0.4):0.7)"
lozupone_newick <- paste0("((", .ag, ":1.1,", .hn, ":1.1):0.3)root;")
lozupone_tree <- ape::read.tree(text=lozupone_newick)
lozupone_tree <- ape::rotateConstr(lozupone_tree, LETTERS[14:1])
plot(lozupone_tree)
usethis::use_data(lozupone_tree)

lozupone_panel_a <- data.frame(
  Species = LETTERS[1:14],
  SampleID = c("Circle", "Square")[c(1, 2, 1, 1, 2, 1, 2, 2, 1, 2, 2, 1, 2, 1)],
  Counts = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
usethis::use_data(lozupone_panel_a)

lozupone_panel_b <- data.frame(
  Species = LETTERS[1:14],
  SampleID = c("Circle", "Square")[c(1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2)],
  Counts = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
usethis::use_data(lozupone_panel_b)
