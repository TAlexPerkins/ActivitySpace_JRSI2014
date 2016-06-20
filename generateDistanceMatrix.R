# function to calculate a distance matrix

generateDistanceMatrix <- function(homes, Locations)
{
  # calculate distances between each home and every other location in the city
  dist <- t(sapply(1:nrow(homes), function(i) {
    sqrt((Locations$xcor - homes$xcor[i])^2 + (Locations$ycor - homes$ycor[i])^2)
  }))
  rownames(dist) = homes$home_code
  colnames(dist) = Locations$location_code

  return(dist)
}
