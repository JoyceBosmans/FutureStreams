library(raster)
library(sf)
library(ggplot2)
#library(rgeos)
#library(ncdf4)
library(inlmisc)
library(foreach)
library(ggpubr)

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
plot(qm)

### Q-max
stack_Q_max <- stack()
for (gcm in climate_models){
	Q_max_hist <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_max_fut  <- raster(paste0(inputdir,'Q-max/Q-max_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_max_diff <- log10(Q_max_fut) - log10(Q_max_hist)	#TO DO: log10
	stack_Q_max <- addLayer(stack_Q_max,Q_max_diff)
}
mean_Q_max <- calc(stack_Q_max, fun=mean)

#mask grid cells with mean Q below cutoff
mean_Q_max <- mean_Q_max * qm

df_Q_max_1 <- rasterToPoints(crop(mean_Q_max,insets_ext[[1]]), spatial=TRUE) %>% as.data.frame()
df_Q_max_2 <- rasterToPoints(crop(mean_Q_max,insets_ext[[2]]), spatial=TRUE) %>% as.data.frame()
df_Q_max_3 <- rasterToPoints(crop(mean_Q_max,insets_ext[[3]]), spatial=TRUE) %>% as.data.frame()

df_Q_max_1$`Q-max` <- df_Q_max_1$Var.1
df_Q_max_2$`Q-max` <- df_Q_max_2$Var.1
df_Q_max_3$`Q-max` <- df_Q_max_3$Var.1

#boundaries / ranges for plotting Q-max
Tbreaks = 0.1
Tlim_up = 0.5
Tlim_lo = -0.5

color_scheme <- inlmisc::GetColors(length(seq(Tlim_lo,Tlim_up,Tbreaks)),scheme='sunset',reverse = T)
br <- seq((Tlim_lo-Tbreaks),(Tlim_up+Tbreaks),Tbreaks)
lab <- paste0(c(paste0('<',Tlim_lo),seq(Tlim_lo,0,Tbreaks),paste0('+',seq(Tbreaks,Tlim_up,Tbreaks)),paste0('>',Tlim_up)),'')

df_Q_max_1$t <- df_Q_max_1$`Q-max`
df_Q_max_1$`Q-max`[df_Q_max_1$`Q-max` > Tlim_up] <- Tlim_up+Tbreaks
df_Q_max_1$`Q-max`[df_Q_max_1$`Q-max` < Tlim_lo] <- Tlim_lo-Tbreaks

Qmax1 <- ggplot() + geom_raster(data = df_Q_max_1 , aes(x = x, y = y, fill=`Q-max`)) + ylab('Q-max') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.3, -.3, 0, "cm")) + scale_color_gradientn(colours = 'rainbow' ,breaks = br,labels = lab, na.value = 'transparent')
Qmax2 <- ggplot() + geom_raster(data = df_Q_max_2 , aes(x = x, y = y, fill=`Q-max`)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm")) + scale_color_gradientn(colors = color_scheme,values=br)
Qmax3 <- ggplot() + geom_raster(data = df_Q_max_3 , aes(x = x, y = y, fill=`Q-max`)) + ylab('') + xlab('')  + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm")) + scale_color_gradientn(colors = color_scheme,values=br)
    
    
### Q-wm
stack_Q_wm <- stack()
for (gcm in climate_models){
	Q_wm_hist <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	Q_wm_fut  <- raster(paste0(inputdir,'Q-wm/Q-wm_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	Q_wm_diff <- log10(Q_wm_fut) - log10(Q_wm_hist)	#TO DO: log10
	stack_Q_wm <- addLayer(stack_Q_wm,Q_wm_diff)
}
mean_Q_wm <- calc(stack_Q_wm, fun=mean)

#mask grid cells with mean Q below cutoff
mean_Q_wm <- mean_Q_wm * qm

df_Q_wm_1 <- rasterToPoints(crop(mean_Q_wm,insets_ext[[1]]), spatial=TRUE) %>% as.data.frame()
df_Q_wm_2 <- rasterToPoints(crop(mean_Q_wm,insets_ext[[2]]), spatial=TRUE) %>% as.data.frame()
df_Q_wm_3 <- rasterToPoints(crop(mean_Q_wm,insets_ext[[3]]), spatial=TRUE) %>% as.data.frame()

df_Q_wm_1$`Q-wm` <- df_Q_wm_1$Var.1
df_Q_wm_2$`Q-wm` <- df_Q_wm_2$Var.1
df_Q_wm_3$`Q-wm` <- df_Q_wm_3$Var.1

Qwm1 <- ggplot() + geom_raster(data = df_Q_wm_1 , aes(x = x, y = y, fill=`Q-wm`)) + ylab('Q-wm') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.3, -.3, 0, "cm"))
Qwm2 <- ggplot() + geom_raster(data = df_Q_wm_2 , aes(x = x, y = y, fill=`Q-wm`)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"))
Qwm3 <- ggplot() + geom_raster(data = df_Q_wm_3 , aes(x = x, y = y, fill=`Q-wm`)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"))
    
### WT-hq
stack_WT_hq <- stack()
for (gcm in climate_models){
	WT_hq_hist <- raster(paste0(inputdir,'WT-hq/WT-hq_',gcm,'_hist_1976-01-31_to_2005-12-31.nc'))
	WT_hq_fut  <- raster(paste0(inputdir,'WT-hq/WT-hq_',gcm,'_rcp8p5_2081-01-31_to_2099-12-31.nc'))
	WT_hq_diff <- WT_hq_fut - WT_hq_hist	
	stack_WT_hq <- addLayer(stack_WT_hq,WT_hq_diff)
}
mean_WT_hq <- calc(stack_WT_hq, fun=mean)

#mask grid cells with mean Q below cutoff
mean_WT_hq <- mean_WT_hq * qm

df_WT_hq_1 <- rasterToPoints(crop(mean_WT_hq,insets_ext[[1]]), spatial=TRUE) %>% as.data.frame()
df_WT_hq_2 <- rasterToPoints(crop(mean_WT_hq,insets_ext[[2]]), spatial=TRUE) %>% as.data.frame()
df_WT_hq_3 <- rasterToPoints(crop(mean_WT_hq,insets_ext[[3]]), spatial=TRUE) %>% as.data.frame()

df_WT_hq_1$`WT-hq` <- df_WT_hq_1$Var.1
df_WT_hq_2$`WT-hq` <- df_WT_hq_2$Var.1
df_WT_hq_3$`WT-hq` <- df_WT_hq_3$Var.1

WThq1 <- ggplot() + geom_raster(data = df_WT_hq_1 , aes(x = x, y = y, fill=`WT-hq`)) + ylab('WT-hq') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.3, -.3, 0, "cm"))
WThq2 <- ggplot() + geom_raster(data = df_WT_hq_2 , aes(x = x, y = y, fill=`WT-hq`)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"))
WThq3 <- ggplot() + geom_raster(data = df_WT_hq_3 , aes(x = x, y = y, fill=`WT-hq`)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"))
    
ggarrange(
    ggarrange(Qmax1, Qmax2, Qmax3, 
          #labels = c("Amazone", "Danube", "Ganges"),
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1
          ),
    ggarrange(Qwm1, Qwm2, Qwm3, 
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1
          ),
    ggarrange(WThq1, WThq2, WThq3, 
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1
          ),
    ncol = 1, nrow = 4)
    

### changing legend: #######################
    
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)} 
  
Qmax1 <- ggplot() + geom_raster(data = df_Q_max_1 , aes(x = x, y = y, fill=Var.1)) + ylab('Q-max') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.3, -.3, 0, "cm"))
Qmax2 <- ggplot() + geom_raster(data = df_Q_max_2 , aes(x = x, y = y, fill=Var.1)) + ylab('') + xlab('') + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"))
Qmax3 <- ggplot() + geom_raster(data = df_Q_max_3 , aes(x = x, y = y, fill=Var.1)) + ylab('') + xlab('')  + theme(
    panel.background = element_blank(),plot.margin = margin(-.1, -.1, -.3, -.1, "cm"),legend.title = element_blank(),legend.text = element_text(size=18),legend.position="bottom") 

Qmax_legend <- get_legend(Qmax3)

ggarrange(
    ggarrange(Qmax1 + theme(legend.position="none"), Qmax2 + theme(legend.position="none"), Qmax3 + theme(legend.position="none"), 
          #labels = c("Amazone", "Danube", "Ganges"),
          common.legend = TRUE, legend = "bottom",
          ncol = 3, nrow = 1
          ),
    ncol = 1, nrow = 3)
    

