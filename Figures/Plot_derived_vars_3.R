library(raster)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(inlmisc)
library(foreach)
library(ggpubr)
library(plyr)

inputdir = '/vol/milkundata/FutureStreams/'	# Q-max, Q-wm, WT-hq, WT-range

crs_custom <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"

crop_main <- c(-180, 180, -60, 90)
cropping_poly <- st_bbox(extent(crop_main), crs = st_crs(4326)) %>% st_as_sfc(.)

world <- rnaturalearth::ne_countries(returnclass = "sf")[,1] %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
bb <- rnaturalearth::ne_download(type = "wgs84_bounding_box", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)
graticules <- rnaturalearth::ne_download(type = "graticules_30", category = "physical",returnclass = "sf") %>%
  st_crop(cropping_poly) %>%
  st_transform(crs_custom)

### Q-mean from hist as mask
climate_models <- c('gfdl','hadgem','ipsl','miroc','noresm')

stack_Q_mean <- stack()
for (gcm in climate_models){
	Q_mean_hist  <- raster(paste0(inputdir,'Q-mean/Q-mean_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	stack_Q_mean <- addLayer(stack_Q_mean,Q_mean_hist)
}
mean_Q_mean <- calc(stack_Q_mean, fun=mean)

cutoff <- 50
qm               <- mean_Q_mean
qm[qm < cutoff]  <- NA
qm[qm >= cutoff] <- 1
#crs(qm) <- 4326

### Q-max ###
stack_Q_max <- stack()
for (gcm in climate_models){
	Q_max_hist <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_max_fut  <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_max_diff <- log10(Q_max_fut) - log10(Q_max_hist)	#TO DO: log10
	stack_Q_max <- addLayer(stack_Q_max,Q_max_diff)
}
mean_Q_max <- calc(stack_Q_max, fun=mean)
mean_Q_max <- mask(mean_Q_max,qm)

### Q-wm ###
stack_Q_wm <- stack()
for (gcm in climate_models){
	Q_wm_hist <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_wm_fut  <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_wm_diff <- log10(Q_wm_fut) - log10(Q_wm_hist)	#TO DO: log10
	stack_Q_wm <- addLayer(stack_Q_wm,Q_wm_diff)
}
mean_Q_wm <- calc(stack_Q_wm, fun=mean)
mean_Q_wm <- mask(mean_Q_wm,qm)

plotvars <- stack(mean_Q_max,mean_Q_wm)
vars     <- c('Q-max','Q-wm')
names(plotvars) <- vars

df <- foreach(i = seq_along(vars),.combine='rbind') %do% {
    as(plotvars[[i]] %>% projectRaster(.,crs=crs_custom) %>% mask(.,bb), "SpatialPixelsDataFrame") %>%
      as.data.frame(.) %>%
      mutate(dervar = vars[i]) %>%
      dplyr::select(x,y,value = colnames(.)[1],dervar)
  }
  
### manually make dfs:

mean_Q_max <- projectRaster(mean_Q_max,crs=crs_custom)
df_q_max <- as.data.frame(mean_Q_max,xy=T,na.rm=T)
colnames(df_q_max) <- c('x','y','value')
df_q_max$dervar <- 'Q-max'

mean_Q_wm <- mask(mean_Q_wm,qm)
mean_Q_wm <- projectRaster(mean_Q_wm,crs=crs_custom)
df_q_wm <- as.data.frame(mean_Q_wm,xy=T,na.rm=T)
colnames(df_q_wm) <- c('x','y','value')
df_q_wm$dervar <- 'Q-wm'

df <- rbind(df_q_max,df_q_wm)

Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')
lim <- c(Tlim_lo-Tbreaks,Tlim_up+Tbreaks)

p <- ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = "black", lwd = NA) +
  geom_raster(data = df_q_max, aes(x=x, y=y, fill=value), alpha=0.8) +
  scale_fill_gradientn(colors = color_scheme,breaks = br,labels = lab,limits = lim,na.value = 'transparent') +
  facet_wrap(~dervar) +
  theme_minimal() +
  theme(text = element_text(size = 12),
        panel.grid.major = element_line(color=NA),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'bottom',
        legend.key.width = unit(6,'line'),
        strip.background = element_rect('white'),
        strip.background.x = element_blank(),
        strip.background.y = element_blank(),
        strip.text = element_text(angle = 0, vjust = -1, size = 12),
        legend.title = element_blank()
  )
######

main <- rasterToPoints(mean_Q_mean %>% mask(.,qm) %>% crop(crop_main) %>% mask(.,bb)) %>%
  as.data.frame() %>%
  rename(qav = layer) %>%
  st_as_sf(coords = c('x','y'))
  
ggplot() +
  geom_sf(data = bb, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = graticules, fill = NA, color = "black", lwd = 0.1) +
  geom_sf(data = world, fill = 'black', lwd = NA) + # '#e4cead' yellow ochre
  geom_point(data = df, aes(x=x, y=y, color=value), shape=15, stroke = 0) +
  facet_wrap(~dervar)+
  scale_fill_gradientn(colors = color_scheme,breaks = br,labels = lab,limits = lim,na.value = 'transparent')
