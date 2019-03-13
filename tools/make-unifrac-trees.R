library(ape)

# To avoid dependence on ape library,
# generate code to make trees for testing


t1 <- "(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;"
cat("t1 <- ")
dput(read.tree(text=t1))

t_nbl <- "((((OTU1,OTU2),OTU3)),(OTU4,OTU5));"
cat("t_no_branch_lengths <-")
dput(read.tree(text=t_nbl))

