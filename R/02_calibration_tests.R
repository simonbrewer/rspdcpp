## Comparison of rcarbon and rcpp calibration

sourceCpp("./src/calibrate.cpp")

load("data/dates.RData")


## RCPP - need to load CalCurve
cc <- read.table("./data/intcal20.14c", header = FALSE, sep = ',')[, 1:3]
colnames(cc) <- c("CALBP","C14BP","Error")
cc <- cc[order(cc$CALBP), ]
dates_rcpp <- calibrate(age, error, cc, 0, 10000)
