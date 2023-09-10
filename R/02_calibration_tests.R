## Comparison of rcarbon and rcpp calibration
library(Rcpp)
sourceCpp("./src/calibrate.cpp")

load("data/dates.RData")
system.time(dates_rcarbon <- 
              rcarbon::calibrate(uncal.samples, errors, calCurves = 'intcal20', 
                                 ids = ids, normalised = TRUE, verbose = FALSE, timeRange = c(10000, 0)))

## RCPP - need to load CalCurve
cc <- read.table("./data/intcal20.14c", header = FALSE, sep = ',')[, 1:3]
colnames(cc) <- c("CALBP","C14BP","Error")
cc <- cc[order(cc$CALBP), ]
## Full date range
system.time(dates_rcpp <- calibrate(uncal.samples, errors, cc, 0, 55000,
                                    normalize = TRUE)) 
plot(dates_rcarbon$grids[[1]])
lines(dates_rcpp[[1]], col = 2)
## Differeces
rev(dates_rcpp[[1]]$prdens) - dates_rcarbon$grids[[1]]$PrDens

## Subset range 0-10000
system.time(dates_rcpp <- calibrate(uncal.samples, errors, cc, 0, 10000,
                                    normalize = TRUE)) 
plot(dates_rcarbon$grids[[1]])
lines(dates_rcpp[[1]], col = 2)
## Differeces
rev(dates_rcpp[[1]]$prdens) - dates_rcarbon$grids[[1]]$PrDens



