#load packages and other files
library(mvtnorm)
library(bbmle)
source("generateDistanceMatrix.R")
source("simulateActivitySpace.R")
normalize = function(x){x / sum(x)}


# load data and estimation results
Locations = read.table(
  'locations.txt',
  col.names=c('location_code','location_type','xcor','ycor'))
load("aspaceSize.RData")
load("timeHome.RData")
load("where.RData")
load("propLoc.RData")
load("freqDurn.RData")


# pick an obviously fake part_code for this hypothetical individual
part_code = c("EXAMPLE1","EXAMPLE2")
home_code = c("lot_001","lot_002")
participants = data.frame(part_code, home_code, stringsAsFactors=FALSE)
home = Locations[Locations$location_code == home_code,]
home$home_code = home$location_code
dist = generateDistanceMatrix(home, Locations)


# simulate individual's activity space
sim_data = simulate_aspace(
  part_codes = 'EXAMPLE',
  Locations,
  dist,
  participants,
  aspaceSize,
  timeHome,
  where,
  propLoc,
  freqDurn,
  aspaceSizeModel = 'Negative binomial',
  propLocModel = "L9",
  whereModel = "L8",
  fdModel = "L8",
  dist_bins = seq(0, max(dist) + 100, 100))
