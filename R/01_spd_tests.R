library(Rcpp)
sourceCpp("./src/makespd.cpp")

load("./data/uncalDates.RData")

test_dates = list(grid = list(uncalsample.norm$grids[[1]],
                              uncalsample.norm$grids[[2]])
)

make_spd(test_dates$grid, 0, 4000)
make_spd(uncalsample.norm$grids, 0, 4000) -> x

