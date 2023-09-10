## Compare uncalibrate methods

load("./data/model.RData")

## RCarbon
system.time(
  uncal_rcarbon <- rcarbon::uncalibrate(sim.growthModel.spd, 
                                        calCurves = 'intcal20', 
                                        verbose = F)
)

## RCPP
sourceCpp("./src/uncalibrate.cpp")
cc <- read.table("./data/intcal20.14c", header = FALSE, sep = ',')[, 1:3]
colnames(cc) <- c("CALBP","C14BP","Error")
cc <- cc[order(cc$CALBP), ]
system.time(
  uncal_rcpp <- uncalibrate(sim.growthModel.spd$CalBP,
                            sim.growthModel.spd$PrDens,
                            cc, 
                            min(cc$C14BP), 
                            max(cc$C14BP))
)

plot(uncal_rcarbon$CRA, rev(uncal_rcpp$cra))
abline(0,1)
plot(uncal_rcarbon$CRA, uncal_rcarbon$PrDens, type = 'l')
lines(uncal_rcpp$cra, uncal_rcpp$prdens, col = 2)

hist(uncal_rcarbon$PrDens - rev(uncal_rcpp$prdens))
