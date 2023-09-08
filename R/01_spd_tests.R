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

plot(spd.uncalsample.norm$grid$calBP, spd.uncalsample.norm$grid$PrDens, 
     type = 'p', xlab = "Age BP", ylab = "D")
lines(spd_test$calbp, spd_test$d, col = 2)

df = data.frame(ages = rev(spd.uncalsample.norm$grid$calBP), 
                rcarbon = rev(spd.uncalsample.norm$grid$PrDens), rcpp = 
                  spd_test$d)
head(df, 10)

## RCARBON speed comparison below here
# library(rcarbon)

# system.time(spd_rcarbon <- spd(uncalsample.norm, timeRange = c(4000,0), spdnormalised = TRUE, verbose = FALSE))

# system.time(spd_rcpp <- make_spd(uncalsample.norm$grids, 0, 4000))


