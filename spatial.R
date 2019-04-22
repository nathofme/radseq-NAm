### Spatial & Environmental Data on North American Starlings

library(raster) # to download data from WorldClim
library(sp)
library(rgdal)

setwd("~/Documents/starlingRAD/analysis/spatial")
BIO1 <- raster("wc2.0_bio_5m_01.tif")
BIO2 <- raster("wc2.0_bio_5m_02.tif")
BIO3 <- raster("wc2.0_bio_5m_03.tif")
BIO4 <- raster("wc2.0_bio_5m_04.tif")
BIO5 <- raster("wc2.0_bio_5m_05.tif")
BIO6 <- raster("wc2.0_bio_5m_06.tif")
BIO7 <- raster("wc2.0_bio_5m_07.tif")
BIO8 <- raster("wc2.0_bio_5m_08.tif")
BIO9 <- raster("wc2.0_bio_5m_09.tif")
BIO10 <- raster("wc2.0_bio_5m_10.tif")
BIO11 <- raster("wc2.0_bio_5m_11.tif")
BIO12 <- raster("wc2.0_bio_5m_12.tif")
BIO13 <- raster("wc2.0_bio_5m_13.tif")
BIO14 <- raster("wc2.0_bio_5m_14.tif")
BIO15 <- raster("wc2.0_bio_5m_15.tif")
BIO16 <- raster("wc2.0_bio_5m_16.tif")
BIO17 <- raster("wc2.0_bio_5m_17.tif")
BIO18 <- raster("wc2.0_bio_5m_18.tif")
BIO19 <- raster("wc2.0_bio_5m_19.tif")
# definitions of WorldClim2.0 variables
#BIO1 = Annual Mean Temperature
#BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
#BIO3 = Isothermality (BIO2/BIO7) (* 100)
#BIO4 = Temperature Seasonality (standard deviation *100)
#BIO5 = Max Temperature of Warmest Month
#BIO6 = Min Temperature of Coldest Month
#BIO7 = Temperature Annual Range (BIO5-BIO6)
#BIO8 = Mean Temperature of Wettest Quarter
#BIO9 = Mean Temperature of Driest Quarter
#BIO10 = Mean Temperature of Warmest Quarter
#BIO11 = Mean Temperature of Coldest Quarter
#BIO12 = Annual Precipitation
#BIO13 = Precipitation of Wettest Month
#BIO14 = Precipitation of Driest Month
#BIO15 = Precipitation Seasonality (Coefficient of Variation)
#BIO16 = Precipitation of Wettest Quarter
#BIO17 = Precipitation of Driest Quarter
#BIO18 = Precipitation of Warmest Quarter
#BIO19 = Precipitation of Coldest Quarter