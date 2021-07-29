library(raster)
library(sf)
library(ggplot2)
library(inlmisc)
library(foreach)
library(ggpubr)

inputdir = '/vol/milkundata/FutureStreams/'	# Q-max, Q-wm, WT-hq, WT-range

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
	Q_mean_hist  <- raster(paste0(inputdir,'Q-mean/Q-mean_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	stack_Q_mean <- addLayer(stack_Q_mean,Q_mean_hist)
}
mean_Q_mean <- calc(stack_Q_mean, fun=mean)

cutoff <- 50
qm               <- mean_Q_mean
qm[qm < cutoff]  <- NA
qm[qm >= cutoff] <- 1

### Q-max
stack_Q_max <- stack()
for (gcm in climate_models){
	Q_max_hist <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_max_fut  <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_max_diff <- log10(Q_max_fut) - log10(Q_max_hist)	#TO DO: log10
	stack_Q_max <- addLayer(stack_Q_max,Q_max_diff)
}
mean_of_var <- calc(stack_Q_max, fun=mean)

inset <- foreach(j = 1:length(insets_ext),.combine = 'rbind') %do% {
  as(mean_of_var %>% mask(.,qm) %>% crop(insets_ext[[j]]), "SpatialPixelsDataFrame") %>%
    as.data.frame(.) %>%
    mutate(inset_no = j) #%>%
    #dplyr::select(x,y,value = V1,inset_no) %>%
    #filter(!is.na(value))
}

#boundaries / ranges for plotting Q-max
Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')

Q_max_plot <- ggplot() +
  geom_tile(data = inset, aes(x=x, y=y, fill=V1), alpha=0.8) +
  scale_fill_gradientn(
    colors = color_scheme,breaks = br,labels = lab,na.value = 'transparent') +
  facet_wrap('inset_no',nrow = 1, scales = 'free') +
  coord_cartesian(expand = F) +
  ylab(' ') +
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
    strip.background = element_rect('white'),
    strip.background.x = element_blank(),
    strip.background.y = element_blank(),
    strip.text = element_blank(),
    legend.title = element_blank(),
    aspect.ratio = 1,
    plot.margin = margin(0, 0, -3, 0, "cm")
  )
