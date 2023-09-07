library(Rcpp)
sourceCpp("./src/makespd.cpp")

load("./data/uncalDates.RData")
load("./data/spdExample.RData")

test_dates = list(grid = list(uncalsample.norm$grids[[1]],
                              uncalsample.norm$grids[[2]])
)

make_spd(test_dates$grid, 0, 4000)

## Results comparison
spd_test <- make_spd(uncalsample.norm$grids, 0, 4000)

plot(spd.uncalsample.norm$grid$calBP, spd.uncalsample.norm$grid$PrDens, type = 'l')
lines(spd_test$calbp, spd_test$d, col = 2)

## RCARBON comparison below here
library(rcarbon)

