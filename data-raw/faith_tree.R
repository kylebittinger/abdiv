faith_newick <- "(((a:5,(b:2,c:1):4):4,(d:5,e:3):1):20);"
faith_tree <- ape::read.tree(text=faith_newick)
faith_tree <- ape::rotateConstr(faith_tree, letters[5:1])
plot(faith_tree)
usethis::use_data(faith_tree, compress = "gzip")
