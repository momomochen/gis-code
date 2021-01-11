
library(sf) 
library(tmap)
library(spdep) 
library(tidyverse) 
library(viridis) 
library(corrplot) 
library(Hmisc) 
library(PerformanceAnalytics)
library(spgwr) 
library(car) 
shp <- st_read('./data/shp/province_all_data_simple.shp')
shp.point <- st_centroid(shp)


p.aqi<- tm_shape(shp = shp) +
  tm_fill("aqi",
          breaks=c(0,50,100,150,200),
          title = 'AQI concentration') +
  tm_borders(alpha = 0.5)+
  tm_layout(title = 'AQI',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))
tmap_save(p.aqi, filename = "chinamainland_aqi.png")



p.pm25 <- tm_shape(shp = shp) +
  tm_fill("pm2_5",
          breaks=c(0,35,75,115,150,250),
          title = 'PM2.5 concentration') +
  tm_borders(alpha = 0.5)+
  tm_layout(title = 'PM2.5',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))
tmap_save(p.pm25, filename = "chinamainland_pm2.5.png")



p.pm10 <- tm_shape(shp = shp) +
  tm_fill("pm10",
          breaks=c(0,50,150,250,350,420),
          title = 'PM10 concentration') +
  tm_borders(alpha = 0.5)+
  tm_layout(title = 'PM10',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))
tmap_save(p.pm10, filename = "chinamainland_pm10.png")



set.ZeroPolicyOption(TRUE)
coords <- coordinates(as(shp.point,'Spatial'))
knn_wards <- knearneigh(shp.point, k=3)# nearest 
LWard_knn <- knn2nb(knn_wards)
plot(LWard_knn, coords, col="red")

Lward.knn_3_weight <- nb2listw(LWard_knn, style="C")

moran.test(shp.point$pm2_5, Lward.knn_3_weight)


moran.test(shp.point$pm10, Lward.knn_3_weight)


moran.test(shp.point$aqi, Lward.knn_3_weight)


local.aqi <- localmoran(x = shp.point$aqi, listw = Lward.knn_3_weight)
moran.map.aqi <- cbind(shp, local.aqi)
tm_shape(moran.map.aqi) +
  tm_polygons(col = 'Ii',
              style = 'quantile',
              title = 'Local Moran Values for AQI')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))


local.pm25 <- localmoran(x = shp.point$pm2_5, listw = Lward.knn_3_weight)
moran.map.pm25 <- cbind(shp, local.pm25)
tm_shape(moran.map.pm25) +
  tm_polygons(col = 'Ii',
              style = 'quantile',
              title = 'Local Moran Values for PM2.5')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))


local.pm10 <- localmoran(x = shp.point$pm10, listw = Lward.knn_3_weight)
moran.map.pm10 <- cbind(shp, local.pm10)
tm_shape(moran.map.pm10) +
  tm_polygons(col = 'Ii',
              style = 'quantile',
              title = 'Local Moran Values for PM10')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))






library(spdep)
lg1 <- localG(shp.point$aqi, listw=Lward.knn_3_weight, zero.policy=T)
shp$lg1 <- lg1[]
tm_shape(shp) +
  tm_polygons(col = 'lg1',
              title = "Getis-Ord Gi*
(z scores)") + 
  tm_layout(title = 'Hotspot Analysis for AQI',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))

lg1 <- localG(shp.point$pm2_5, listw=Lward.knn_3_weight, zero.policy=T)
shp$lg1 <- lg1[]
tm_shape(shp) +
  tm_polygons(col = 'lg1',
              title = "Getis-Ord Gi*
(z scores)") + 
  tm_layout(title = 'Hotspot Analysis for PM2.5',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))

lg1 <- localG(shp.point$pm10, listw=Lward.knn_3_weight, zero.policy=T)
shp$lg1 <- lg1[]
tm_shape(shp) +
  tm_polygons(col = 'lg1',
              title = "Getis-Ord Gi*
(z scores)") + 
  tm_layout(title = 'Hotspot Analysis for PM10',title.position = c("center","top"))+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))





df <- shp %>% select(aqi,气温,降水,湿度,人均生) %>% st_drop_geometry()
chart.Correlation(df, histogram=TRUE, pch=19, method = 'pearson') # see the correlation
names(df) <- c('aqi','tem','pre','rhu','gdp')

model <- lm(aqi ~ tem + pre + rhu + gdp, data = df) # OLS model
summary(model)
AIC(model)
plot(model) 
durbinWatsonTest(model) 
shp$residuals <- residuals(model)
p.res <- tm_shape(shp = shp) +
  tm_polygons("residuals")
tmap_save(p.res, filename = "ols_residual.png")





shp.sp <- as(shp.point, 'Spatial')
GWRbandwidth <- gwr.sel(aqi ~ 气温+降水+湿度+人均生, data=shp.sp, adapt=T, method="cv") #calculate kernel bandwidth

GWRModel <- gwr(aqi ~ 气温+降水+湿度+人均生, data=shp.sp, adapt=GWRbandwidth,
                hatmatrix = TRUE,se.fit = TRUE)
GWRModel 



results<-as.data.frame(GWRModel$SDF)
names(results)

shp$tem_coef <- results$气温 
shp$pre_coef <- results$降水  
shp$rhu_coef <- results$湿度 
shp$gdp_coef <- results$人均生  

tm_shape(shp) +
  tm_polygons(col = 'tem_coef',
              title ='Coefficient of 
Temperature')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))

tm_shape(shp) +
  tm_polygons(col = 'pre_coef',
              title ='Coefficient of 
Precipitation')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))


tm_shape(shp) +
  tm_polygons(col = 'rhu_coef',
              title ='Coefficient of 
Humidity')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))


tm_shape(shp) +
  tm_polygons(col = 'gdp_coef',
              title ='Coefficient of 
GDP per capita')+
  tm_scale_bar(position=c("right","bottom"),text.size=0.4)+
  tm_compass(type="4star",position=c("left","top"))


