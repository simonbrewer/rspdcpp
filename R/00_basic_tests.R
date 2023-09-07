library(Rcpp)
sourceCpp("./src/list2mat.cpp")

mylist <- list(c(3,4,99,1,222), c(1,2,3,4,5))

make_mat(mylist)
print_list(mylist)

myorder <- list(c(5,4,3,2,1) - 1, c(5,4,3,2,1) - 1)

make_mat_order(mylist, myorder, 0, 10)
make_spd(mylist, myorder, 0, 10)
