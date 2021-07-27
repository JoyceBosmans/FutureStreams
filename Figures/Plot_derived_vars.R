library(raster)
library(sf)
library(ggplot2)
#library(rgeos)
#library(ncdf4)
library(inlmisc)

inputdir = '/vol/milkundata/FutureStreams/'	# Q-max, Q-wm, WT-hq, WT-range

# bbox for 3 insets
insets_c <- data.frame(
  x = c(-57,20,88.5),
  y = c(-2.5,45,24.5)
)

#1 -57.257 -2.622	## -57.2083, -2.625 Amazone
#2 20.132 45.210	## 20.125, 45.2083 Danube
#3 88.380 24.376	## 88.375 24.375 Ganges

insets_side = 4
insets_ext <- foreach(i = 1:nrow(insets_c)) %do% extent(c(insets_c[i,'x']-insets_side, insets_c[i,'x']+insets_side,
                                                          insets_c[i,'y']-insets_side, insets_c[i,'y']+insets_side))
                                                          
####

climate_models <- c('gfdl','hadgem','ipsl','miroc','noresm')

stack_Q_max <- stack()
for (gcm in climate_models){
	Q_max_hist <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_max_fut  <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_max_diff <- Q_max_fut - Q_max_hist	#TO DO: log10
	stack_Q_max <- addLayer(stack_Q_max,Q_max_diff)
}
mean_Q_max <- calc(stack_Q_max, fun=mean)

mean_Q_max_1 <- crop(mean_Q_max,insets_ext[[1]])
plot(mean_Q_max_1)

par(mfrow=c(1,3))		#To DO: set to 4*3 for all 4 variables

df_Q_max_1 <- rasterToPoints(crop(mean_Q_max,insets_ext[[1]]), spatial=TRUE) %>% as.data.frame()
ggplot() + geom_raster(data = df_Q_max_1 , aes(x = x, y = y, fill=layer)) 

df_Q_max_2 <- rasterToPoints(crop(mean_Q_max,insets_ext[[2]]), spatial=TRUE) %>% as.data.frame()
ggplot() + geom_raster(data = df_Q_max_2 , aes(x = x, y = y, fill=layer)) 

df_Q_max_3 <- rasterToPoints(crop(mean_Q_max,insets_ext[[3]]), spatial=TRUE) %>% as.data.frame()
ggplot() + geom_raster(data = df_Q_max_3 , aes(x = x, y = y, fill=layer)) 

# TO DO
# apply Q 10 m3/s as cutoff


### TRYING STUFF:

Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')

df_Q_max_1 <- rasterToPoints(mean_Q_max_1, spatial=TRUE) %>% as.data.frame()
ggplot() + geom_raster(data = df_Q_max_1 , aes(x = x, y = y, fill=layer))

