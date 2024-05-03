
library(stringr)
library(ggplot2)
library(dplyr)
library(rgdal)
library(rnaturalearth)
library(data.table)
library(ggpubr)
# library(grid)

as.num <- function(x) { x <- as.numeric(as.character(x))}

world <- ne_countries(scale = "medium", returnclass = "sf")
four_cols = c( '#F22300', '#3C9AB2','#5E0079', '#739D1C')

############# load Ecoregions shapefile ############# 

shp_ecoregions <- readOGR(dsn = "data/ICES_ecoregions/ICES_ecoregions_20171207_erase_ESRI.shp")
shp_ecoregions <- spTransform(shp_ecoregions,
                              CRS("+proj=longlat +ellps=WGS84 +datum=WGS84
                                     +no_defs"))
df_ecoregion <- fortify(shp_ecoregions)
df_ecoregion$OBJECTID <- as.num(df_ecoregion$id) + 1

df_ecoregion2 <- as.data.frame(shp_ecoregions)
df_ecoregion2$OBJECTID <- as.num(df_ecoregion2$OBJECTID)

df_ecoregion3 <- merge(data.table(df_ecoregion[,c("long","lat","OBJECTID", "order","piece")]),
                       data.table(df_ecoregion2[,c("OBJECTID", "Ecoregion")]), by = "OBJECTID")
df_ecoregion3 <- df_ecoregion3 %>% dplyr::arrange(Ecoregion,order)

############# load and process taxonomic data ############# 
df_env <- read.csv2(file = 'data/raw/df_env.csv')
df_env <- df_env %>% dplyr::select(-X)

df_PCA <- read.csv2(file = 'data/raw/df_PCA.csv')
df_PCA <- df_PCA %>% dplyr::select(-X)

df_PCA <- merge(df_PCA, df_env)

df_trend_PCA <- read.csv2(file = 'data/intermediate/df_temporal_trends_PCA.csv')
df_trend_PCA <- df_trend_PCA %>% dplyr::select(-X)

############# load and process functional  data ############# 

df_MFA <- read.csv2(file = 'data/raw/df_MFA.csv')
df_MFA <- df_MFA %>% dplyr::select(-X)

df_MFA <- merge(df_MFA, df_env)

df_trend_MFA <- read.csv2(file = 'data/intermediate/df_temporal_trends_MFA.csv')
df_trend_MFA <- df_trend_MFA %>% dplyr::select(-X)

################################################################################
################### temporal evolution by bioregion  ###########################
################################################################################

seuil_pval = 0.05
for(rep_region in unique(df_trend_PCA$Ecoregion )){
  
  pl1_celtic_sea_taxo <- df_pca_taxo_sm_melt[  df_pca_taxo_sm_melt$Ecoregion == rep_region &
                                                 df_pca_taxo_sm_melt$dimension %in% c("PC1", "PC2"), ]
  
  pl1_celtic_sea_taxo$facet <- "Taxonomic"
  
  pl1_celtic_sea_fonctio <- df_trend_MFA[  df_trend_MFA$Ecoregion == rep_region&
                                             df_trend_MFA$dimension %in% c("MF1", "MF2"), ]
  pl1_celtic_sea_fonctio$facet <- "Functional"
  
  pl1_celtic_sea_both <- rbind(pl1_celtic_sea_taxo, pl1_celtic_sea_fonctio)
  
  pl1_celtic_sea_both$dimension <- factor(pl1_celtic_sea_both$dimension, exclude = T)
  
  pl1_celtic_sea_both$facet <- as.factor(pl1_celtic_sea_both$facet)
  pl1_celtic_sea_both$facet  <- factor(pl1_celtic_sea_both$facet, levels = c( "Taxonomic", "Functional"))
  

  a <- ggplot(pl1_celtic_sea_both[pl1_celtic_sea_both$pval <= seuil_pval, ]) +
    scale_fill_manual(drop = FALSE, values = four_cols) + 
    scale_color_manual(drop = FALSE, values = four_cols)  + 
    
    geom_point(aes(x = Year, y = mean_value, color = dimension   ),
               size = 1, alpha = 0.7) +
    
    geom_smooth(aes(x = Year, y = mean_value, color = dimension, fill = dimension),
                method = 'lm', alpha = 0.2) +
    
    
    geom_point(data = pl1_celtic_sea_both[pl1_celtic_sea_both$pval >= seuil_pval, ],
               aes(x = Year, y = mean_value, color = dimension),
               size = 2, alpha = 0.1) +
    
    facet_grid(facet ~ Ecoregion, scales = 'free_y') + 
    theme_classic() + xlim(1993, 2022.3) + ylab("Mean values") + 
    
    theme(panel.background = element_rect(fill = NA, color = 'black'),
          legend.title = element_blank(),
          panel.grid.major = element_line(color = 'grey85'),
          panel.grid.minor = element_line(color = 'grey95'),
          legend.position = "none",
          axis.title = element_blank())
  a
  
  
  rep_region2 <- str_replace_all(rep_region, " ", "_")
  rep_region2 <- str_replace_all(rep_region2, "\n", "_")
  rep_region2 <- str_replace_all(rep_region2, "-", "_")
  
  assign(paste0("plt_", rep_region2), a, .GlobalEnv)
  print(paste0("plt_", rep_region2))
}

# plt_Adriatic_Sea
# plt_Aegean_Levantine_Sea
# plt_Baltic_Sea
# plt_Bay_of_Biscay_and_Iberian_Coast
# plt_Celtic_Seas
# plt_Greater_North_Sea
# plt_Ionian_Sea_and_Central_Med._Sea
# plt_Western_Med._Sea


df_legend <- pl1_celtic_sea_both[c(1,2,99,100), ]
df_legend$dimension <- factor(df_legend$dimension, exclude = T)
levels(df_legend$dimension)

plt_all <- ggplot(data = world) +
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = '#006FC4', alpha = 0.5,
               aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
  geom_point(data = df_legend, 
             aes(x = 1, y = 0, color = dimension)) + 
  scale_color_manual(drop = FALSE, values = four_cols)  + 
  guides(color=guide_legend(ncol=2)) +
  
  theme(legend.title =  element_blank(),
        legend.position = c(0.65, 0.55))  +
  
  annotation_custom( ggplotGrob(plt_Bay_of_Biscay_and_Iberian_Coast), 
                     xmin = -49, xmax = -19, ymin = 27, ymax = 46) +
  
  annotation_custom( ggplotGrob(plt_Western_Med._Sea), 
                     xmin = -22, xmax = 8, ymin = 14, ymax = 32) + 
  annotation_custom( ggplotGrob(plt_Ionian_Sea_and_Central_Med._Sea), 
                     xmin = 7, xmax = 35, ymin = 14, ymax = 32) + 
  
  
  annotation_custom( ggplotGrob(plt_Celtic_Seas), 
                     xmin = -49, xmax = -19, ymin = 46, ymax = 65) + 
  
  annotation_custom( ggplotGrob(plt_Aegean_Levantine_Sea), 
                     xmin = 32, xmax = 62, ymin = 27, ymax = 46) + 
  annotation_custom( ggplotGrob(plt_Adriatic_Sea), 
                     xmin = 32, xmax = 62, ymin = 46, ymax = 65) +
  
  
  annotation_custom( ggplotGrob(plt_Greater_North_Sea), 
                     xmin = -21, xmax = 9, ymin = 62, ymax = 82) +
  annotation_custom( ggplotGrob(plt_Baltic_Sea), 
                     xmin = 8, xmax = 36, ymin = 62, ymax = 82) +
  
  geom_segment(aes(x = 16, y = 43, xend = 32,  yend = 55), size = 1.5, color = '#00406F')  +
  geom_segment(aes(x = 25, y = 37, xend = 32,  yend = 42), size = 1.5, color = '#00406F') +
  
  geom_segment(aes(x = 18, y = 36, xend = 19,  yend = 30), size = 1.5, color = '#00406F')  +
  geom_segment(aes(x = 8, y = 39, xend = 0,  yend = 30), size = 1.5, color = '#00406F') + 
  
  geom_segment(aes(x = -10, y = 41, xend = -23,  yend = 35), size = 1.5, color = '#00406F')  +
  geom_segment(aes(x = -10, y = 55, xend = -23,  yend = 55), size = 1.5, color = '#00406F') + 
  
  geom_segment(aes(x = 4, y = 57, xend = -3,  yend = 68), size = 1.5, color = '#00406F')  +
  geom_segment(aes(x = 20, y = 57, xend = 22,  yend = 65), size = 1.5, color = '#00406F') +
  
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold'),
        plot.margin = margin(7.3,6,7,6, "cm")) 


ggsave(file = 'figures/figure3.jpg', plot = plt_all,
       width = 3.5,  height = 3.2, scale = 3)
