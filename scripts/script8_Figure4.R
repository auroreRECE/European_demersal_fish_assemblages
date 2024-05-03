library(scales)

# library(stringr)
# library(ggplot2)
# library(dplyr)
# library(rgdal)
# library(rnaturalearth)
# library(data.table)
# library(ggpubr)
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

df_all_taxo <- read.csv2(file = 'data/intermediate/df_all_taxo.csv')

############# load and process functional  data ############# 

df_MFA <- read.csv2(file = 'data/raw/df_MFA.csv')
df_MFA <- df_MFA %>% dplyr::select(-X)

df_MFA <- merge(df_MFA, df_env)

df_trend_MFA <- read.csv2(file = 'data/intermediate/df_temporal_trends_MFA.csv')
df_trend_MFA <- df_trend_MFA %>% dplyr::select(-X)

df_all_func <- read.csv2(file = 'data/intermediate/df_all_func.csv')

################################################################################
################### temporal evolution by bioregion  ###########################
################################################################################


param_lineheit <- 0.6

df_all2$dim_legend <-  as.factor(ifelse(df_all2$facet == "Taxonomic",
                                        substr(df_all2$dim, 1, 3), 
                                        substr(df_all2$dim, 5, 8) ))
df_all3 <-  df_all2[ df_all2$dim %in% c("PC1/MF1", "PC2/MF2"), ]

df_all3$dim_legend <- factor(df_all3$dim_legend ,  exclude = T,
                             levels = c(  "PC1",  "PC2", "MF1" , "MF2"))
levels(df_all3$dim_legend)


levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "chloro_std_scale"] <- "Var. surf. chloro." 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "mlotst_std_scale"] <- "Var. mixed layer depth" 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "oxy_bottom_std_scale"] <- "Var. bottom oxy." 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "sal_surf_scale"] <- "Surface salinity" 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "sal_surf_std_scale"] <- "Var. surf. sal" 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "temp_surf_std_scale"] <- "Var. surf. temp." 
levels(df_all3$variable_env)[levels(df_all3$variable_env) ==  "temp_bottom_std_scale"] <- "Var. bottom temp." 


for(rep_region in unique(df_all3$Ecoregion )){
  print(rep_region)
  dede <-  df_all3[df_all3$Ecoregion == rep_region, ]
  
  if(dim(dede)[1] != 0){

    dede$variable_env2 <- dede$variable_env
    dede$variable_env2 <- ifelse(dede$variable_env2 == "Mixed layer depth", 
                                 "Mixed layer\ndepth", 
                                 str_replace_all(dede$variable_env2, " ", "\n"))
    dede$variable_env2 <- as.factor( dede$variable_env2)
    levels(dede$variable_env2)
    levels(dede$variable_env2)[levels(dede$variable_env2) ==  "Var.\nsurf.\nchloro."] <- "Var. surf.\nchloro." 
    levels(dede$variable_env2)[levels(dede$variable_env2) ==  "Var.\nbottom\noxy."] <- "Var. bottom\noxy." 
    levels(dede$variable_env2)[levels(dede$variable_env2) ==  "Var.\nbottom\ntemp."] <- "Var. bottom\ntemp." 
    levels(dede$variable_env2)[levels(dede$variable_env2) ==  "Var.\nmixed\nlayer\ndepth"] <- "Var. mixed\nlayer depth"
  
    plt_region <- ggplot(dede, 
                         aes(x = Value, y = variable_env, color = variable_env)) + 
      geom_pointrange(aes(xmin =Value -  std_error, xmax = Value + std_error,
                          fill = variable_env), shape=21, fatten = 1.1, size = 1,
                      position = position_dodge2(width = 0.35)) + 
      geom_text(data = dede, aes(x = Value, y = variable_env, label = variable_env2),
                vjust = -0.2,  
                lineheight = param_lineheit) + 
      theme_classic() + 
      scale_color_manual(values = hue_pal()(12), drop = FALSE)  +
      scale_fill_manual(values = hue_pal()(12),drop = FALSE)  +
      facet_nested( facet + dim_legend ~  Ecoregion, scales = 'free', space='free_y') + 
      xlab("Slope") +
      theme(panel.background = element_rect(color = "black", fill = NA),
            panel.grid.major.y = element_line(color = 'grey90'),
            legend.title = element_blank(),
            axis.title = element_blank(),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    
    rep_region2 <- str_replace_all(rep_region, " ", "_")
    rep_region2 <- str_replace_all(rep_region2, "\n", "_")
    rep_region2 <- str_replace_all(rep_region2, "-", "_")
    
    assign(paste0("plt2_", rep_region2), plt_region, .GlobalEnv)
    print(paste0("plt2_", rep_region2))
  }
  
  
}

# plt2_Adriatic_Sea
# plt2_Aegean_Levantine_Sea
# plt2_Baltic_Sea
# plt2_Bay_of_Biscay_and_Iberian_Coast
# plt2_Celtic_Seas
# plt2_Ionian_Sea_and_Central_Med._Sea
# plt2_Spatial_distribution

plt_all_temp <- ggplot(data = world) +
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = '#006FC4', alpha = 0.5,
               aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
  
  theme(legend.title =  element_blank(),
        legend.position = c(0.65, 0.55))  +
  
  
  annotation_custom( ggplotGrob(plt2_Aegean_Levantine_Sea), 
                     xmin = 17.5, xmax = 49.5, ymin = 23, ymax = 34) + 
  
  annotation_custom( ggplotGrob(plt2_Adriatic_Sea), 
                     xmin = 31, xmax = 61, ymin = 33.5, ymax = 47) +
  
  
  annotation_custom( ggplotGrob(plt2_Ionian_Sea_and_Central_Med._Sea), 
                     xmin = -15, xmax = 17, ymin = 23, ymax = 34) + 
  
  annotation_custom( ggplotGrob(plt2_Bay_of_Biscay_and_Iberian_Coast), 
                     xmin = -46, xmax = -14, ymin = 32, ymax = 48) +
  annotation_custom( ggplotGrob(plt2_Celtic_Seas), 
                     xmin = -46, xmax = -14, ymin = 48, ymax = 65) + 
  
  annotation_custom( ggplotGrob(plt2_Baltic_Sea), 
                     xmin = 31, xmax = 61,  ymin = 47, ymax = 65) +
  
theme(axis.title = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.position = 'none',
      plot.title = element_text(face = 'bold'),
      plot.margin = margin(7.3,6,7,6, "cm")) 


ggsave(file = 'figures/taxo_fonctio/figure4_all_temporel.jpg', plot = plt_all_temp,
       width = 3.5,  height = 3.2, scale = 3)


########

df_spatial <-  df_all3[df_all3$Ecoregion == "Spatial distribution", ]
df_spatial$dim <- factor(df_spatial$dim, exclude = T)

df_spatial$variable_env2 <- df_spatial$variable_env

df_spatial$variable_env2 <- str_replace_all(df_spatial$variable_env2, " ", "\n")
df_spatial$variable_env2 <- as.factor(df_spatial$variable_env2)
levels(df_spatial$variable_env2)
levels(df_spatial$variable_env2)[levels(df_spatial$variable_env2) ==  "Mixed\nlayer\ndepth"] <- "Mixed\nlayer depth" 
levels(df_spatial$variable_env2)[levels(df_spatial$variable_env2) ==  "Var.\nbottom\noxy."] <- "Var. bottom\noxy." 
levels(df_spatial$variable_env2)[levels(df_spatial$variable_env2) ==  "Var.\nbottom\ntemp."] <- "Var. bottom\ntemp." 
levels(df_spatial$variable_env2)[levels(df_spatial$variable_env2) ==  "Var.\nmixed\nlayer\ndepth"] <- "Var. mixed\nlayer depth"
levels(df_spatial$variable_env2)[levels(df_spatial$variable_env2) ==  "Var.\nsurf.\nchloro."] <- "Var. surf.\nchloro."


plt2_spa <- ggplot(df_spatial, 
                  aes(x = Value, y = variable_env, color = variable_env)) + 
  geom_vline(xintercept = 0 , color = 'grey50', size = 0.1) +
  
  geom_pointrange(aes(xmin =Value -  std_error, xmax = Value + std_error,
                      fill = variable_env), shape=21, fatten = 1.1, size = 1,
                  position = position_dodge2(width = 0.35)) + 
  geom_text(data = df_spatial, aes(x = Value, y = variable_env, label = variable_env2),
            vjust = -0.3,  
            lineheight = param_lineheit) + 
  theme_classic() + 
  scale_color_manual(values = hue_pal()(12), drop = FALSE)  +
  scale_fill_manual(values = hue_pal()(12),drop = FALSE)  +
  # xlim(-1.15, 0.2) +
  facet_nested(~ Ecoregion + facet + dim_legend , scales = 'free', space='free_y') + 
  xlab("Slope") +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = 'grey90'),
        legend.title = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
plt2_spa

ggsave(file = 'figures/figure4_all_spatial.jpg', plot = plt2_spa,
       width = 3.5,  height = 1.2, scale = 3)

