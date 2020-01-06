leprieur_newick <- "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"
leprieur_tree <- ape::read.tree(text=leprieur_newick)
plot(leprieur_tree)
usethis::use_data(leprieur_tree, compress = "gzip")
