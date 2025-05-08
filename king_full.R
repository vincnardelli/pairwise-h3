library(sf)
library(dplyr)
library(spatialreg)
library(spdep)
library(h3jsr)


h3_res <- 7

data <- read_sf("kingcounty/kc_house.shp")

data$price <- data$price - mean(data$price)
data$sqft_liv <- data$sqft_liv - mean(data$sqft_liv)

coords <- st_coordinates(data)
nb<-knn2nb(knearneigh(coords, k=5, longlat = T))
listw <- nb2listw(nb)

model <- errorsarlm(price ~ sqft_liv, data=data, listw=listw)
model <- GMerrorsar(price ~ sqft_liv, data=data, listw=listw)

summary(model)


#summary(model)
