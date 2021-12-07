library(raster)
library(sf)
library(ggplot2)
library(RColorBrewer)
library(inlmisc)
library(foreach)
library(ggpubr)
library(plyr) 

inputdir = '/vol/milkundata/FutureStreams/'	# Q-max, Q-wm, WT-hq, WT-range

### inset function
inset_function <- function(mean_of_var){
	foreach(j = 1:length(insets_ext),.combine = 'rbind') %do% {
  as(mean_of_var %>% mask(.,qm) %>% crop(insets_ext[[j]]), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(inset_no = j) #%>%
    #dplyr::select(x,y,value = V1,inset_no) %>%
    #filter(!is.na(value))
}}

### plot function
plot_function <- function(inset,color_scheme,br,lab,lim,ylab_text){
  ggplot() +
  geom_tile(data = inset, aes(x=x, y=y, fill=V1), alpha=0.8) +
  scale_fill_gradientn(colors = color_scheme,breaks = br,labels = lab,limits = lim,na.value = 'transparent') +
  facet_wrap('inset_no',nrow = 1, scales = 'free') + 
  coord_cartesian(expand = F) +
  ylab(ylab_text) + 
  theme_minimal(base_size = 50, base_rect_size = 50) +
  theme(
    #text = element_text(size = 50),
    panel.background = element_rect(fill = 'black'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.5, "lines"),
    # axis.text = element_blank(),		# theme legend
    axis.ticks = element_line(color = 'black',size=1),
    axis.title.x = element_blank(),
    axis.text=element_text(size=12),
    legend.position = 'bottom',
    legend.key.width = unit(9,'line'),
    legend.key.height = unit(2,'line'),
    legend.margin=margin(-25,-15,-15,-15),
    legend.text=element_text(size=24),
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1,
    plot.margin = margin(0, 1, 0, 0, "cm")
  ) }

# bbox for 3 insets
insets_c <- data.frame(
  x = c(-57,20,88.5),
  y = c(-2.5,45,24.5)
)

insets_side = 4
insets_ext <- foreach(i = 1:nrow(insets_c)) %do% extent(c(insets_c[i,'x']-insets_side, insets_c[i,'x']+insets_side,
                                                          insets_c[i,'y']-insets_side, insets_c[i,'y']+insets_side))
                                                          
### Q-mean from hist as mask
climate_models <- c('gfdl','hadgem','ipsl','miroc','noresm')

stack_Q_mean <- stack()
for (gcm in climate_models){
	Q_mean_hist  <- raster(paste0(inputdir,'Q-mean/Q-mean_',gcm,'_hist_1976-01-31_to_2005-12-31.nc4'))
	stack_Q_mean <- addLayer(stack_Q_mean,Q_mean_hist)
}
mean_Q_mean <- calc(stack_Q_mean, fun=mean)

cutoff <- 10
qm               <- mean_Q_mean
qm[qm < cutoff]  <- NA
qm[qm >= cutoff] <- 1

### Q-max ###
stack_Q_hist <- stack()
stack_Q_fut  <- stack()
for (gcm in climate_models){
	Q_max_hist <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_max_fut  <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	stack_Q_hist <- addLayer(stack_Q_hist,Q_max_hist)
	stack_Q_fut  <- addLayer(stack_Q_fut,Q_max_fut)
}
mean_of_var <- log10(calc(stack_Q_fut, fun=mean)) - log10(calc(stack_Q_hist, fun=mean))

#boundaries / ranges for plotting Q-max
Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5

mean_of_var[mean_of_var > Tlim_up] <- Tlim_up+Tbreaks
mean_of_var[mean_of_var < Tlim_lo] <- Tlim_lo-Tbreaks

inset <- inset_function(mean_of_var)

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')
lim <- c(Tlim_lo-Tbreaks,Tlim_up+Tbreaks)

Q_max_plot <- plot_function(inset,color_scheme,br,lab,lim,'Q-max') 

stack_Q_hist <- stack()
stack_Q_fut  <- stack()
for (gcm in climate_models){
	Q_mean_hist  <- raster(paste0(inputdir,'Q-mean/Q-mean_',gcm,'_hist_1976-01-31_to_2005-12-31.nc4'))
	Q_mean_fut   <- raster(paste0(inputdir,'Q-mean/Q-mean_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc4'))
	stack_Q_hist <- addLayer(stack_Q_hist,Q_mean_hist)
	stack_Q_fut  <- addLayer(stack_Q_fut,Q_mean_fut)
}
mean_of_var <- log10(calc(stack_Q_fut, fun=mean)) - log10(calc(stack_Q_hist, fun=mean))

inset <- inset_function(mean_of_var)

# plot using same color scheme as Q-max
Q_mean_plot <- plot_function(inset,color_scheme,br,lab,lim,'Q-mean')

### Q-wm ###
stack_Q_hist <- stack()
stack_Q_fut  <- stack()
for (gcm in climate_models){
	Q_wm_hist <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_wm_fut  <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	stack_Q_hist <- addLayer(stack_Q_hist,Q_wm_hist)
	stack_Q_fut  <- addLayer(stack_Q_fut,Q_wm_fut)
}
mean_of_var <- log10(calc(stack_Q_fut, fun=mean)) - log10(calc(stack_Q_hist, fun=mean))

inset <- inset_function(mean_of_var)

# plot using same color scheme as Q-max
Q_wm_plot <- plot_function(inset,color_scheme,br,lab,lim,'Q-wm')

### Q-dm ###
stack_Q_hist <- stack()
stack_Q_fut  <- stack()
for (gcm in climate_models){
	Q_dm_hist <- raster(paste0(inputdir,'Q-dm/Q-dm_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_dm_fut  <- raster(paste0(inputdir,'Q-dm/Q-dm_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	stack_Q_hist <- addLayer(stack_Q_hist,Q_dm_hist)
	stack_Q_fut  <- addLayer(stack_Q_fut,Q_dm_fut)
}
mean_of_var <- log10(calc(stack_Q_fut, fun=mean)) - log10(calc(stack_Q_hist, fun=mean))

inset <- inset_function(mean_of_var)

# plot using same color scheme as Q-max
Q_dm_plot <- plot_function(inset,color_scheme,br,lab,lim,'Q-dm')

### WT-hq ###
stack_WT_hist <- stack()
stack_WT_fut  <- stack()
for (gcm in climate_models){
	WT_hq_hist <- raster(paste0(inputdir,'WT-hq/WT-hq_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	WT_hq_fut  <- raster(paste0(inputdir,'WT-hq/WT-hq_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	stack_WT_hist <- addLayer(stack_WT_hist,WT_hq_hist)
	stack_WT_fut  <- addLayer(stack_WT_fut,WT_hq_fut)
}
mean_of_var <- calc(stack_WT_fut, fun=mean) - calc(stack_WT_hist, fun=mean)

Tbreaks = 0.25
Tlim_up = 5
Tlim_lo = -5

mean_of_var[mean_of_var > Tlim_up] <- Tlim_up+Tbreaks
mean_of_var[mean_of_var < Tlim_lo] <- Tlim_lo-Tbreaks

inset <- inset_function(mean_of_var)

#plot
color_scheme <- inlmisc::GetColors(length(Tlim_lo:Tlim_up),scheme='sunset')
br <- seq((Tlim_lo-1),(Tlim_up+1),1)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,1),paste0('+',seq(1,Tlim_up,1)),paste0('>',Tlim_up)),'°C')
lim <- c(Tlim_lo-1,Tlim_up+1)

WT_hq_plot <- plot_function(inset,color_scheme,br,lab,lim,'WT-hq')

### WT-cq ###
stack_WT_hist <- stack()
stack_WT_fut  <- stack()
for (gcm in climate_models){
	WT_cq_hist <- raster(paste0(inputdir,'WT-cq/WT-cq_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	WT_cq_fut  <- raster(paste0(inputdir,'WT-cq/WT-cq_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	stack_WT_hist <- addLayer(stack_WT_hist,WT_cq_hist)
	stack_WT_fut  <- addLayer(stack_WT_fut,WT_cq_fut)
}
mean_of_var <- calc(stack_WT_fut, fun=mean) - calc(stack_WT_hist, fun=mean)

mean_of_var[mean_of_var > Tlim_up] <- Tlim_up+Tbreaks
mean_of_var[mean_of_var < Tlim_lo] <- Tlim_lo-Tbreaks

inset <- inset_function(mean_of_var)

WT_cq_plot <- plot_function(inset,color_scheme,br,lab,lim,'WT-cq')

### WT-range ###
stack_WT_range <- stack()
for (gcm in climate_models){
	WT_range_hist <- raster(paste0(inputdir,'WT-range/WT-range_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	WT_range_fut  <- raster(paste0(inputdir,'WT-range/WT-range_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	WT_range_diff <- WT_range_fut - WT_range_hist	
	stack_WT_range <- addLayer(stack_WT_range,WT_range_diff)
}
mean_of_var <- calc(stack_WT_range, fun=mean)

mean_of_var[mean_of_var > Tlim_up] <- Tlim_up+Tbreaks
mean_of_var[mean_of_var < Tlim_lo] <- Tlim_lo-Tbreaks

inset <- inset_function(mean_of_var)

# plot using same color settings as WT-hq
WT_range_plot <- plot_function(inset,color_scheme,br,lab,lim,'WT-range')

### WT-mean ###
stack_WT_mean <- stack()
for (gcm in climate_models){
	WT_mean_hist <- raster(paste0(inputdir,'WT-mean/WT-mean_',gcm,'_hist_1976-01-31_to_2005-12-31.nc4'))
	WT_mean_fut  <- raster(paste0(inputdir,'WT-mean/WT-mean_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc4'))
	WT_mean_diff <- WT_mean_fut - WT_mean_hist	
	stack_WT_mean <- addLayer(stack_WT_mean,WT_mean_diff)
}
mean_of_var <- calc(stack_WT_mean, fun=mean)

mean_of_var[mean_of_var > Tlim_up] <- Tlim_up+Tbreaks
mean_of_var[mean_of_var < Tlim_lo] <- Tlim_lo-Tbreaks

inset <- inset_function(mean_of_var)

WT_mean_plot <- plot_function(inset,color_scheme,br,lab,lim,'WT-mean')

#~ all_plot <- ggarrange(Q_max_plot, Q_wm_plot, WT_hq_plot, WT_range_plot, 
#~           ncol = 1, nrow = 4
#~           )
          
#~ ggsave('figure4.pdf',all_plot,width = 200,height = 263,dpi = 300,units = 'mm', scale = 5)#, limitsize = F)

ggsave('figure4-Q-max.png',Q_max_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-Q-wm.png',Q_wm_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-Q-dm.png',Q_dm_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-Q-mean.png',Q_mean_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)

ggsave('figure4-WT-hq.png',WT_hq_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-WT-cq.png',WT_cq_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-WTrange.png',WT_range_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
ggsave('figure4-WT-mean.png',WT_mean_plot,width = 75,height = 35,dpi = 300,units = 'mm', scale = 5, limitsize=F)
