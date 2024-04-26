################################################################################
######## Aurore Receveur
######## 26/10/2022 
################################################################################

two_cols = c( '#F22300', '#3C9AB2')
three_cols = c('#21908C','#5E0079', '#F92867')
four_cols = c( '#F22300', '#3C9AB2','#5E0079', '#739D1C')

memory.limit(size=9999999999999)

eight_cols <- c("#F8766D","#F99D1E","#FF61CC",
                "#C77CFF","#00BFC4","#74D33A",
                "#7C4728", "#5D6966")

## coucou 

################################################################################
######################### chargement packages ##################################
################################################################################

source("R_scripts/21_12_13_script1_base.R")
source("R_scripts/22_03_28_script7_function_extract_env.R")

library(sf)
library(maps)
library(rnaturalearthdata)
library(rnaturalearth)

library(raster)

library(ggpubr)
library(dplyr)
library(taxize)
library(stringr)

library(treemap)
library(ggrepel)

library(readxl)

library(data.table)
library(stars)

library(ncdf4)

library(mFD)

library(vegan)
library(NbClust)
library(GGally)
library(grid)

library("FactoMineR")
library("factoextra")
world <- ne_countries(scale = "medium", returnclass = "sf")
library(adespatial )


library(lme4)
library(brms)
library(bayesplot)
library(openxlsx)
library(worms)
library(rfishbase)
library(mgcv)

library(ade4)
library("rphylopic")

set.seed(23)


funct_my_sampling = function(x) {
  if (length(x) >= 2) {  x = x[sample(1:length(x), 2, replace = FALSE)]
  paste0(x, collapse = "|")
  } else { NA }

}

shp_ecoregions <- readOGR(dsn = "data/raw/others/ICES_ecoregions/ICES_ecoregions_20171207_erase_ESRI.shp")
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
# 
# raster_commun <- raster(ext = extent(-44, 69, 30, 85),
#                         res = c(1,0.5),
#                         crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# 
# letters_double <- paste0(rep(LETTERS, each = 26), rep(LETTERS, 26))[1:110]
# 
# raster_commun$my_spatial_id  <- as.factor(paste0(rep(letters_double, each = 113),
#                                                  rep(1:113, 110)))
# 
# raster_commun_df <-   as.data.frame(raster_commun, xy = TRUE)
# 
# 
# 
# shp_lme <- readOGR(dsn = "data/LME66/LMEs66.shp")
# shp_lme <- spTransform(shp_lme,
#                               CRS("+proj=longlat +ellps=WGS84 +datum=WGS84
#                                      +no_defs"))
# df_lme <- fortify(shp_lme)
# 
# maplme <- ggplot(data = world) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-20, 40), ylim = c( 31 , 63.5),expand = FALSE) +
#     geom_polygon(data = df_lme, fill = NA, size = 0.5, aes(x = long, y = lat,
#                      group =  interaction(id, piece)), color = 'black') +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   geom_label(aes(x = 20, y = 58, label = 'Baltic Sea'), color = '#F8766D') +
#   geom_label(aes(x = 2, y = 58, label = 'North Sea'), color = '#F89D18') +
#   geom_label(aes(x = -7, y = 50, label = 'Celtic-Biscay Shelf'), color = '#FF61CC') +
#   geom_label(aes(x = -10, y = 40, label = 'Iberian Coastal'), color = '#C77CFF') +
#   geom_label(aes(x = 20, y = 35, label = 'Mediterranean Sea'), color = '#00BFC4') +
# 
#   theme(axis.title = element_blank(),
#         panel.background = element_rect(color = 'black'))
# maplme
# ggsave(file = 'figures/taxo_fonctio/map_lme.jpg', plot = maplme,
#        width = 1.5,  height = 1.2, scale = 3)



################################################################################
######################### chargement des donnees  ###########################
################################################################################

load(file = 'data/processed/df_total2.Rdata')

sp_valid_names <- read.table(file = 'data/raw/bio/species_valid_names_complete.txt')


df_total2 <- merge(df_total2, sp_valid_names, by = 'genus_sp')
df_total2 <- df_total2 %>% 
  dplyr::select(- genus_sp) %>% 
  dplyr::mutate(genus_sp = valid_name) %>% 
  dplyr::select(- valid_name)

df_total2$Survey <- as.character(df_total2$Survey)
df_total2$Survey <- ifelse(df_total2$database == "MEDITS", "MEDITS", df_total2$Survey)
df_total2 <- df_total2 %>% 
  dplyr::filter(spawning.type != 'NA' & feeding.mode != 'NA'  & 
                  habitat != 'NA' & Year >= 1983 & 
                  abundance > 0)
head(df_total2)

df_total2$spawning.type <- factor(df_total2$spawning.type,
                                  exclude = TRUE)
df_total2$feeding.mode <- factor(df_total2$feeding.mode,
                                 exclude = TRUE)
df_total2$habitat <- factor(df_total2$habitat,
                            exclude = TRUE)

n_distinct(df_total2$genus_sp )
n_distinct(df_total2$family )
n_distinct(df_total2$genus )
n_distinct(df_total2$Survey )

n_distinct(df_total2$my_spatial_id )
n_distinct(df_total2$index )

range(df_total2$Year)


# dede <- df_total2 %>% 
#   dplyr::filter(genus_sp == "Dussumieria_elopsoides")
# mapdata <- ggplot(data = world) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-20, 40), ylim = c( 31 , 63.5),expand = FALSE) +
#   geom_point(data = dede, aes(x = Lon, y = Lat, color = Survey), size = 3) +
#   
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   theme(axis.title = element_blank())
# ggsave(file = 'figures/taxo_fonctio/map_Dussumieria_elopsoides.jpg', plot = mapdata,
#        width = 3.4,  height = 1.9, scale = 3)



# df_total2 <- df_total2 %>%
#   dplyr::filter(Year %between% c(1994, 2019)
#                 & Ecoregion != "Oceanic Northeast Atlantic")


# dede <- df_total2 %>%
#   dplyr::group_by(Year, Quarter, my_spatial_id) %>%
#   dplyr::summarise(n = n_distinct(index))
# 
# df_genus_sp_traits <- df_total2 %>%
#   dplyr::group_by(genus_sp, habitat, spawning.type, feeding.mode) %>%
#   dplyr::summarise(length.max = mean(length.max),
#                    length.infinity = mean(length.infinity),
#                    length.maturity = mean(length.maturity),
#                    age.maturity = mean(age.maturity),
#                    growth.coefficient = mean(growth.coefficient),
#                    tl = mean(tl),
#                    offspring.size = mean(offspring.size),
#                    fecundity = mean(fecundity))
# write.csv2(df_genus_sp_traits, file = "data/processed/species_traits.csv",
#            row.names = FALSE)

# df_map <- df_total2 %>%
#   dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic") %>%
#   dplyr::group_by(Lon, Lat, Survey ) %>%
#   dplyr::summarise(n_sp = n_distinct(genus_sp))
# 
# 
# df_ecoregion3 <- df_ecoregion3 %>%
#   dplyr::filter(Ecoregion %in% unique(df_total2$Ecoregion))
# 
# # df_ecoregion3_s <- df_ecoregion3 %>%
# #   dplyr::group_by(Ecoregion) %>%
# #   dplyr::summarise(lon_m = mean(long),
# #                    lat_m = mean(lat))
# # df_ecoregion3_s <- data.frame(df_ecoregion3_s)
# # df_ecoregion3_s$Ecoregion <- as.factor(df_ecoregion3_s$Ecoregion)
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Bay of Biscay and the Iberian Coast", 'lon_m'] <- -12
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Bay of Biscay and the Iberian Coast", 'lat_m'] <- 45.3
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Greater North Sea", 'lat_m'] <- 61
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Greater North Sea", 'lon_m'] <- 8
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Celtic Seas", 'lat_m'] <- 60.8
# # 
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Baltic Sea", 'lat_m'] <- 59
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Ionian Sea and the Central Mediterranean Sea", 'lon_m'] <- 17
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Ionian Sea and the Central Mediterranean Sea", 'lat_m'] <- 33.5
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Aegean-Levantine Sea", 'lon_m'] <- 29.5
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Aegean-Levantine Sea", 'lat_m'] <- 33
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Adriatic Sea", 'lon_m'] <- 16
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Adriatic Sea", 'lat_m'] <- 46
# # 
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Oceanic Northeast Atlantic", 'lat_m'] <- 49
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Oceanic Northeast Atlantic", 'lon_m'] <- -16
# # 
# # df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Western Mediterranean Sea", 'lat_m'] <- 37
# # 
# # 
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean\nLevantine Sea"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Oceanic Northeast Atlantic" ] <- "Oceanic\nNortheast\nAtlantic"
# # levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea"
# 
# 
# mapdata <- ggplot(data = world) +
#   geom_point(data = df_map, aes(x = Lon, y = Lat, color = Survey), size = 0.1) +
#   # geom_tile(data = raster_commun_df, aes(x = x, y = y), color = 'grey60', fill = NA) +
#   geom_sf(color = 'white', fill = 'grey70', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-20, 40), ylim = c( 31 , 63.5),expand = FALSE) +
#   guides(color = guide_legend(override.aes = list(size = 3))) +
#   # geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, aes(x = long, y = lat,
#   #                  group =  interaction(Ecoregion, piece)), color = 'black') +
#   # geom_label(data = df_ecoregion3_s, aes(x = lon_m, y = lat_m, label = Ecoregion)) +
#   theme(axis.title = element_blank(),
#         legend.position = "none")
# ggsave(file = 'figures/taxo_fonctio/map_Etienne.jpg', plot = mapdata,
#        width = 3.4,  height = 1.9, scale = 3)
# 
# df_bar <- df_total2 %>%
#   dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic") %>%
#   dplyr::group_by(Survey ) %>%
#   dplyr::summarise(n_haul = n_distinct(index))
# barplot <- ggplot(data = df_bar) +
#   geom_bar(aes(x = Survey, y = n_haul, fill = Survey), stat = "identity") +
#   theme_classic() + ylab('Total number of hauls') +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1))
# map_barplot <- ggarrange(mapdata , barplot, common.legend = TRUE,
#                             legend = 'none', ncol = 2, 
#                          widths = c(2, 1.5),
#                          labels = c("A", 'B'))
# 
# 
# ggsave(file = 'figures/taxo_fonctio/figure0.jpg', plot = map_barplot,
#        width = 3,  height = 1.4, scale = 3)
################################################################################
###############selection based on different criterion  ###########################
################################################################################

df_total2$yy_qq_carre = paste0(df_total2$Year, df_total2$Quarter, df_total2$my_spatial_id)
n_distinct(df_total2$yy_qq_carre)

df_small <- df_total2 %>% 
  dplyr::group_by(index,database,Survey,Country,Ship,
                  haul_number,Rect,swept_area,distance,hauling_duration,
                  Year,month,day,Quarter,Lon,Lat, Ecoregion, 
                  my_spatial_id, x_my_spatial_id , y_my_spatial_id, 
                  yy_qq_carre) %>% 
  dplyr::summarise(n_species = n_distinct(genus_sp))

n_distinct(df_small$yy_qq_carre)
df_small <- data.frame(df_small)

sel_rows = tapply(1:nrow(df_small),
                  list(df_small$Year, df_small$Quarter, df_small$my_spatial_id),
                  funct_my_sampling)
sel_rows = as.vector(sel_rows)
sel_rows = sel_rows[!is.na(sel_rows)]
sel_rows = unlist(strsplit(sel_rows, "\\|"))
sel_rows = as.numeric(sel_rows)

df_small_sampled <- df_small[sel_rows, ]


df_all_sampled1 <- df_total2 %>% 
  dplyr::filter(index %in% unique(df_small_sampled$index))

n_distinct(df_total2$genus_sp)
n_distinct(df_all_sampled1$genus_sp)

df_test <- df_all_sampled1 %>%
  dplyr::group_by(Year, Quarter, my_spatial_id) %>%
  dplyr::summarise(n_trawl = n_distinct(index))
table(df_test$n_trawl)

################################################################################
###############selection based on different criterion  ###########################
################################################################################

df_small1 <- df_all_sampled1 %>%
  dplyr::group_by(Year, Quarter, x_my_spatial_id, y_my_spatial_id,
                  my_spatial_id) %>%
  dplyr::summarise(n_sp = n_distinct(genus_sp),
                   abundance_m = mean(abundance),
                   abundance_tot = sum(abundance))  %>%
  mutate(Quarter = paste0("Q", Quarter))

df_small2 <- df_small1 %>%
  dplyr::group_by(Quarter, my_spatial_id, x_my_spatial_id, y_my_spatial_id) %>%
  dplyr::summarise(n_years = n_distinct(Year),
                   year_min = min(Year))

df_for_critere <- df_small2 %>%
  dplyr::filter(year_min <= 2005  & n_years >= 10)%>%
  dplyr::mutate(qq_carre = paste0(Quarter, my_spatial_id))

length(unique(df_for_critere$my_spatial_id))
table(df_for_critere$n_years)


df_all_sampled1$qq_carre <- paste0("Q",df_all_sampled1$Quarter,
                                   df_all_sampled1$my_spatial_id)

df_all_sampled2 <- df_all_sampled1 %>%
  dplyr::filter(qq_carre %in% unique(df_for_critere$qq_carre))

n_distinct(df_all_sampled2$my_spatial_id)
n_distinct(df_all_sampled2$Year)

dim(df_all_sampled2)

################################################################################
######################### chargement des donnees env ###########################
################################################################################

df_for_extract <- df_all_sampled2 %>%
  dplyr::mutate(Lon = x_my_spatial_id, 
                Lat = y_my_spatial_id) %>% 
  dplyr::group_by(Year, Quarter, Lon, Lat) %>% 
  dplyr::summarise(n_species = n_distinct(genus_sp))
df_for_extract$index <- c(1:dim(df_for_extract)[1])

df_for_extract2 <- data.frame(df_for_extract) %>%
  dplyr::mutate(x_my_spatial_id = Lon,
                y_my_spatial_id = Lat) %>%
  dplyr::select(Year, Quarter, x_my_spatial_id, y_my_spatial_id, index)

df_all_sampled2  <- df_all_sampled2 %>% 
  dplyr::mutate(index_old = index) %>% 
  dplyr::select(-index)
df_all_sampled2 <- merge(df_all_sampled2,df_for_extract2 ,
                         by = c("Year", "Quarter", "x_my_spatial_id", "y_my_spatial_id"))

liste_var_env <- c("chloro",  "mlotst", "ph",
                   "temp_bottom", "temp_surf",
                   "oxy_bottom", "oxy_surf",
                   "sal_bottom", "sal_surf",
                   "ugo_bottom", "vgo_bottom", 
                   "ugo_surf", "vgo_surf", "zo")
length(liste_var_env)

sapply(liste_var_env, fun_extract_env, res_temporal = 'year_quarter',
       rep_data = 'df_for_extract')


new_df_all <- merge(new_df_chloro_year_quarter, new_df_mlotst_year_quarter)
new_df_all <- merge(new_df_all, new_df_oxy_bottom_year_quarter)
new_df_all <- merge(new_df_all, new_df_oxy_surf_year_quarter)
new_df_all <- merge(new_df_all, new_df_sal_bottom_year_quarter)
new_df_all <- merge(new_df_all, new_df_sal_surf_year_quarter)
new_df_all <- merge(new_df_all, new_df_temp_bottom_year_quarter)
new_df_all <- merge(new_df_all, new_df_temp_surf_year_quarter)
new_df_all <- merge(new_df_all, new_df_ugo_bottom_year_quarter)
new_df_all <- merge(new_df_all, new_df_vgo_bottom_year_quarter)
new_df_all <- merge(new_df_all, new_df_ugo_surf_year_quarter)
new_df_all <- merge(new_df_all, new_df_vgo_surf_year_quarter)
new_df_all <- merge(new_df_all, new_df_zo_year_quarter)

new_df_all$curr_surf_mea <- sqrt(new_df_all$ugo_surf_mea^2 + new_df_all$vgo_surf_mea^2)
new_df_all$curr_bottom_mea <- sqrt(new_df_all$ugo_bottom_mea^2 + new_df_all$vgo_bottom_mea^2)


name_fich_topo <- paste0("C:/aurore/maestro/3-datasets/env/gebco_2021_from_drive.nc")
ncvar_topo <- nc_open(name_fich_topo)

lon <- ncvar_get(ncvar_topo, 'lon')
lat <- ncvar_get(ncvar_topo, 'lat')
toto <- ncvar_get(ncvar_topo, 'elevation')

dimnames(toto) <- list(lon, lat)

ras_t <- raster(rotate(rotate(rotate(toto))),
                xmn = min(lon),xmx = max(lon), ymn = min(lat), ymx = max(lat),
                CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# plot(ras_t)


df_for_extract2 <- data.frame(df_for_extract) %>%  dplyr::select(index, Lon, Lat)
names(df_for_extract2)[c(2,3)] <- c("x", "y")
coordinates(df_for_extract2) <- ~ x + y
proj4string(df_for_extract2) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

df_for_extract2$depth <- raster::extract(ras_t, df_for_extract2)

df_depth <- as.data.frame(reshape2::melt(toto))
names(df_depth ) <- c("lon", "lat", "depth")
df_depth <- df_depth %>% 
  dplyr::filter(depth <= 200)

df_depth$lon2 <- change.res(df_depth$lon, 0.1)
df_depth$lat2 <- change.res(df_depth$lat, 0.1)

df_depth2 <- df_depth %>% 
  dplyr::group_by(lon2, lat2) %>% 
  dplyr::summarize(depth_mean = mean(depth),
                   depth_min =  min(depth),
                   depth_max =  max(depth)) %>% 
  dplyr::mutate(depth_span = depth_max - depth_min)
save(df_depth2, file = 'data/raw/others/df_depth.Rdata')

a = ggplot(data = world) +
  geom_tile(data = df_depth2, aes(x = lon2, y = lat2, fill = -depth_mean, color = -depth_mean)) +
  geom_sf(color = 'white', fill = 'grey50', size = 0.01) +  theme_classic() +
  coord_sf(xlim = c(-15, 35), ylim = c( 32 , 65),expand = TRUE)+
  scale_fill_viridis_c(name = "Depth (m)", direction = -1) + 
  scale_color_viridis_c(name = "Depth (m)", direction = -1) +
  annotation_scale(width_hint = 0.1, location = "tr") +
  theme(axis.title = element_blank(),
        panel.background = element_rect(color = 'black', fill = NA, size = 1),
        axis.line.x.top = element_line(color = 'black',size = 1),
        axis.line.y.right = element_line(color = 'black',size = 1))
ggsave(file = 'figures/depth.jpg', plot = a,
       width = 1.7,  height = 1.5, scale = 3)



matrix_depth_span <- reshape2::acast(df_depth2, 
                                     lon2 ~ lat2, 
                                     value.var = "depth_span")

ras_depth_span <- raster(rotate(rotate(rotate(matrix_depth_span))),
                xmn = min(df_depth2$lon2),xmx = max(df_depth2$lon2), 
                ymn = min(df_depth2$lat2), ymx = max(df_depth2$lat2),
                CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
df_for_extract2$depth_span <- raster::extract(ras_depth_span, df_for_extract2)



df_topo <- as.data.frame(df_for_extract2)

df_topo1 <- df_topo %>% dplyr::select(index, depth, depth_span)

######################### correction positives

df_no_depth <- df_topo %>%
  dplyr::filter(depth >= 0) %>%
  dplyr::group_by(index,  x, y) %>%
  dplyr::summarise(n_species = length(unique(depth))) %>%
  dplyr::select(-n_species)
dim(df_no_depth)
df_no_depth <- data.frame(df_no_depth)
df_no_depth$depth2 <- NA

toto2 <- replace(toto, which(toto >= 0), NA)


liste_lat_ras <- as.num(colnames(toto2))
range(liste_lat_ras)
liste_lon_ras <- as.num(rownames(toto2))
range(liste_lon_ras)

for(i in 1:dim(df_no_depth)[1]){
  print(i)
  x_lon <- df_no_depth[i, 'x']
  x_lat <- df_no_depth[i, 'y']
  
  res_lat = colnames(toto2)[which(abs(liste_lat_ras - x_lat) == 
                                    min(abs(liste_lat_ras - x_lat)))][1]
  
  ligne_ras <- toto2[,which(colnames(toto2) == res_lat)]
  
  liste_lon_non_null <- ligne_ras[which(!is.na(ligne_ras))] 
  names(liste_lon_non_null) 
  
  liste2_lon <- as.num(names(liste_lon_non_null))
  range(liste2_lon)
  
  res_lon = names(liste_lon_non_null)[which(abs(liste2_lon - x_lon) == 
                                              min(abs(liste2_lon - x_lon)))][1]
  res_depth = liste_lon_non_null[which(abs(liste2_lon - x_lon) == 
                                         min(abs(liste2_lon - x_lon)))][1]
  
  df_topo1[df_topo1$index == df_no_depth[i, 'index'], 'depth'] <- res_depth
}

summary(df_topo1)

######################### all together

df_last <- merge(df_all_sampled2, df_topo1, by = 'index', all.x = TRUE)
df_last <- merge(df_last, new_df_all, by = 'index', all.x = TRUE)

df_last$depth <- -df_last$depth

df_last$log_abundance <- log10(df_last$abundance + 1)

save(df_last, file = "data/processed/df_last.Rdata")

################################################################################
######################### chargement des donnees last ##########################
################################################################################

load(file = "data/processed/df_last.Rdata")
n_distinct(df_last$genus_sp)
n_distinct(df_last$index_old)


# df_last <- df_all_sampled2
df_last$remove <- ifelse((df_last$Quarter == 4 & df_last$Year <= 2002 & 
                            df_last$Ecoregion != "Greater North Sea") |
                           ( df_last$Year < 1997 & 
                            df_last$Ecoregion == "Greater North Sea") |
                           (df_last$Quarter == 1 & 
                              df_last$Ecoregion == "Bay of Biscay and the Iberian Coast") | 
                           (df_last$Ecoregion == "Celtic Seas" &
                              df_last$Year <= 2002) , 
                         'yes', 'no')

table(df_last$remove)
df_last <- df_last %>% 
  dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic") %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = TRUE)) %>% 
  dplyr::filter(remove == 'no') %>% 
  dplyr::select(-remove)

n_distinct(df_last$genus_sp)
n_distinct(df_last$index)

unique(df_last$Survey)
n_distinct(df_last$my_spatial_id)
dede2 <- df_last %>%
  dplyr::filter(Survey != "BTS-VIII")
n_distinct(dede2$my_spatial_id)


head(df_last)
n_distinct(df_last$genus_sp) ## index = index unique par year, quarter, ICES
n_distinct(df_last$index_old) ## index_old = index unique par chalut

df_last %>% 
  dplyr::group_by(Survey) %>% 
  dplyr::summarise(min_yy = min(Year),
                   max_yy = max(Year))

df_total2$remove <- ifelse((df_total2$Quarter == 1 & 
                              df_total2$Ecoregion == "Bay of Biscay and the Iberian Coast")  , 
                         'yes', 'no')
df_number <- df_total2 %>% 
  dplyr::filter(remove == "no") %>% 
  dplyr::filter(!is.na(Ecoregion)) %>% 
  dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic") %>% 
  dplyr::group_by(Ecoregion, Year, Quarter) %>% 
  dplyr::summarise(n = n_distinct(index),
                   n_ices = n_distinct(my_spatial_id))
levels(df_number$Ecoregion)[levels(df_number$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and Central Med. Sea" 
levels(df_number$Ecoregion)[levels(df_number$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and Iberian Coast" 
levels(df_number$Ecoregion)[levels(df_number$Ecoregion) == "Western Mediterranean Sea" ] <-"Western Med. Sea" 
levels(df_number$Ecoregion)[levels(df_number$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine Sea"
levels(df_number$Ecoregion)[levels(df_number$Ecoregion) == "Greater North Sea" ] <- "Greater North Sea" 
df_number$sampling <- "Before sampling"


df_number2 <- df_last %>% 
  dplyr::group_by(Ecoregion, Year, Quarter) %>% 
  dplyr::summarise(n = n_distinct(index_old),
                   n_ices = n_distinct(my_spatial_id))
levels(df_number2$Ecoregion)[levels(df_number2$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and Central Med. Sea" 
levels(df_number2$Ecoregion)[levels(df_number2$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and Iberian Coast" 
levels(df_number2$Ecoregion)[levels(df_number2$Ecoregion) == "Western Mediterranean Sea" ] <-"Western Med. Sea" 
levels(df_number2$Ecoregion)[levels(df_number2$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine Sea"
levels(df_number2$Ecoregion)[levels(df_number2$Ecoregion) == "Greater North Sea" ] <- "Greater North Sea" 
df_number2$sampling <- "After sampling"

df_number3 <- rbind(df_number, df_number2)
df_number3$Quarter <- as.factor(df_number3$Quarter)

df_number3$sampling <- as.factor(df_number3$sampling)
df_number3$sampling <- factor(df_number3$sampling, levels = c( "Before sampling",
                                                               "After sampling" ))

  
  a = ggplot(df_number3,
           aes(x = Year, y = n, shape = Quarter , color = sampling)) + 
  geom_path() +geom_point() + facet_wrap( ~ Ecoregion, scales = 'free', ncol = 4) + 
  scale_color_manual(values = rev(two_cols), name = '') +
  ylab("Number of trawls") +
  # ggtitle("Before sampling" ) +
  theme_classic() + theme(   panel.background = element_rect(color = 'black'),
                          panel.grid.major = element_line(colour = "grey80"))
a

df_number2_m <- df_number2 %>% 
  dplyr::group_by(Ecoregion, Quarter) %>% 
  dplyr::summarise(n_m = mean(n),
                   sd_m = sd(n))
df_number2 <- merge(df_number2, df_number2_m)
df_number2$min <- df_number2$n_m - df_number2$sd_m
df_number2$max <- df_number2$n_m + df_number2$sd_m

# bbb <- ggarrange(plt_1, plt_2, plt_3, plt_4, ncol = 4, nrow =2, 
#                  plt_5, plt_6,plt_7, plt_8, legend = 'right', common.legend = TRUE)
# bbb
ggsave(file = 'figures/taxo_fonctio/figureSX_sampling.jpg', plot =a,
       width = 3.8,  height = 2.5, scale = 3)

# b = ggplot(df_number2, aes(x = Year, y = n, color = as.factor(Quarter) )) +
#   geom_path() + facet_wrap(Ecoregion~ ., scales = 'free', ncol = 4) +
#   ylab("Number of trawls") + ggtitle("After sampling" ) +
#   theme_classic() + theme(legend.position = "none",
#                           panel.background = element_rect(color = 'black'),
#                           panel.grid.major = element_line(colour = "grey80"))
# 
# c = ggarrange(a, b, nrow = 2)

# 
# b = ggplot(df_number3,
#            aes(x = Year, y = n_ices, shape = Quarter , color = sampling)) + 
#   geom_path() +geom_point() + facet_wrap( ~ Ecoregion, scales = 'free', ncol = 4) + 
#   scale_color_manual(values = two_cols, name = '') +
#   ylab("Number of spatial grid cell") +
#   # ggtitle("Before sampling" ) +
#   theme_classic() + theme(   panel.background = element_rect(color = 'black'),
#                              panel.grid.major = element_line(colour = "grey80"))
# b
# ggsave(file = 'figures/taxo_fonctio/figureSX_sampling_ices_grid_cell.jpg', plot = b,
#        width = 3.8,  height = 1.5, scale = 3)
# 
# 
# df_for_extract <- data.frame(df_last) %>% 
#   dplyr::group_by(my_spatial_id, x_my_spatial_id , y_my_spatial_id ) %>% 
#   dplyr::summarise(n_sp = n_distinct(genus_sp)) %>% 
#   dplyr::select(-n_sp)
# names(df_for_extract)[c(2,3)] <- c("x", "y")
# coordinates(df_for_extract) <- ~ x + y
# proj4string(df_for_extract) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# df_for_extract$lme <- raster::extract(shp_lme , df_for_extract)
# df_lme <-  data.frame(my_spatial_id = df_for_extract$my_spatial_id ,
#                       lme = df_for_extract$lme$LME_NAME)
# 
# df_lme_missing <- df_lme[is.na(df_lme$lme), ]
# 
# 
# 
# unique(df_lme$lme)
# df_complement_lme <- df_last %>% 
#   dplyr::filter(my_spatial_id %in% unique(df_lme_missing$my_spatial_id)) %>% 
#   dplyr::group_by(my_spatial_id, Ecoregion) %>% 
#   dplyr::summarise(n = n_distinct(index)) %>% 
#   dplyr::select(-n)
# df_complement_lme <- data.frame(df_complement_lme)
# df_complement_lme[df_complement_lme$my_spatial_id == "CF55", "Ecoregion"] <- "Greater North Sea"
# df_complement_lme[df_complement_lme$my_spatial_id == "DO61", "Ecoregion"] <- "Ionian Sea and the Central Mediterranean Sea"
# df_complement_lme[df_complement_lme$my_spatial_id == "DQ57", "Ecoregion"] <- "Ionian Sea and the Central Mediterranean Sea"
# 
# df_complement_lme <- df_complement_lme %>% 
#   dplyr::group_by(my_spatial_id, Ecoregion) %>% 
#   dplyr::summarise(n = n_distinct(Ecoregion)) %>% 
#   dplyr::select(-n)
# 
# df_complement_lme$lme <-  as.factor(df_complement_lme$Ecoregion)
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) ==  "Greater North Sea"] <- "North Sea"   
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) == "Celtic Seas"] <- "Celtic-Biscay Shelf"
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) == "Western Mediterranean Sea"] <-  "Mediterranean Sea" 
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) ==  "Ionian Sea and the Central Mediterranean Sea" ] <-  "Mediterranean Sea" 
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) ==  "Adriatic Sea"] <-  "Mediterranean Sea" 
# levels(df_complement_lme$lme)[levels(df_complement_lme$lme) ==   "Aegean-Levantine Sea"] <-  "Mediterranean Sea" 
# levels(df_complement_lme$lme)
# 
# 
# df_lme <- df_lme %>% dplyr::filter(!is.na(lme))
# df_lme <- rbind(df_lme, df_complement_lme[,names(df_lme)])
# liste_biscay <- unique(df_lme[df_lme$lme == "Bay of Biscay and the Iberian Coast", 'my_spatial_id'])
# 
# df_biscay <- df_last %>% 
#   dplyr::filter(my_spatial_id %in% liste_biscay) %>% 
#   dplyr::group_by(my_spatial_id, x_my_spatial_id, y_my_spatial_id) %>% 
#   dplyr::summarise(n = n_distinct(index)) %>% 
#   dplyr::select(-n)
#   unique(df_lme$lme)
#   
# df_lme[df_lme$my_spatial_id =="CW41", 'lme' ] <- "Celtic-Biscay Shelf"
# df_lme[df_lme$my_spatial_id =="CZ40", 'lme' ] <- "Celtic-Biscay Shelf"
# df_lme[df_lme$my_spatial_id =="DA41", 'lme' ] <- "Celtic-Biscay Shelf"
# df_lme[df_lme$my_spatial_id =="DB41", 'lme' ] <- "Celtic-Biscay Shelf"
# 
# df_lme[df_lme$my_spatial_id =="DF36", 'lme' ] <- "Iberian Coastal" 
# df_lme[df_lme$my_spatial_id =="DF40", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DF42", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DH36", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DJ36", 'lme' ] <-"Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DK36", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DP36", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DQ36", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DR36", 'lme' ] <- "Iberian Coastal"
# df_lme[df_lme$my_spatial_id =="DR37", 'lme' ] <- "Iberian Coastal"
# save(df_lme, file = "data/processed/df_lme.Rdata")
############# carte 
df_ecoregion3 <- df_ecoregion3 %>%
  dplyr::filter(Ecoregion %in% unique(df_last$Ecoregion))

df_ecoregion3_s <- df_ecoregion3 %>%
  dplyr::group_by(Ecoregion) %>%
  dplyr::summarise(lon_m = mean(long),
                   lat_m = mean(lat))
df_ecoregion3_s <- data.frame(df_ecoregion3_s)
df_ecoregion3_s$Ecoregion <- as.factor(df_ecoregion3_s$Ecoregion)

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Bay of Biscay and the Iberian Coast", 'lon_m'] <- -13
df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Bay of Biscay and the Iberian Coast", 'lat_m'] <- 46

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Greater North Sea", 'lat_m'] <- 61
df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Greater North Sea", 'lon_m'] <- 8

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Celtic Seas", 'lat_m'] <- 60.8


df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Baltic Sea", 'lat_m'] <- 59

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Ionian Sea and the Central Mediterranean Sea", 'lat_m'] <- 34

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Aegean-Levantine Sea", 'lon_m'] <- 25.6
df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Aegean-Levantine Sea", 'lat_m'] <- 43

df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Adriatic Sea", 'lon_m'] <- 16
df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Adriatic Sea", 'lat_m'] <- 46



df_ecoregion3_s[df_ecoregion3_s$Ecoregion ==  "Western Mediterranean Sea", 'lat_m'] <- 37


levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and Central Med. Sea"
levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast"
levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea"
levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean\nLevantine\nSea"
levels(df_ecoregion3_s$Ecoregion)[levels(df_ecoregion3_s$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea"

df_ecoregion3_s <- df_ecoregion3_s %>%
  dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic")
df_ecoregion3 <- df_ecoregion3 %>%
  dplyr::filter(Ecoregion != "Oceanic Northeast Atlantic")


df_map2 <- df_last %>%
  dplyr::group_by(Lon, Lat, Survey ) %>%
  dplyr::summarise(n_sp = n_distinct(genus_sp))
mapdata <- ggplot(data = world) +
  geom_tile(data = raster_commun_df, aes(x = x, y = y), color = 'grey60', fill = NA) +
  geom_sf(color = 'white', fill = 'grey50', size = 0.01) +  theme_classic() +
  coord_sf(xlim = c(-21, 31), ylim = c( 33.2 , 65),expand = TRUE) +
  geom_point(data = df_map2, aes(x = Lon, y = Lat, color = Survey),
             show.legend = FALSE, size = 0.15) +
  geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = '#006FC4', alpha = 0.5,
               aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
  scale_x_continuous(sec.axis = dup_axis()) +
  geom_label(data = df_ecoregion3_s, aes(x = lon_m, y = lat_m, label = Ecoregion)) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(color = 'black', fill = NA, size = 1),
        axis.line.x.top = element_line(color = 'black',size = 1),
        axis.line.y.right = element_line(color = 'black',size = 1))

ggsave(file = 'figures/taxo_fonctio/figure0_sampled.jpg', plot = mapdata,
       width = 2.2,  height = 2, scale = 3)
# 
# 
# # df_map_initial <- df_total2 %>%
#   dplyr::filter(Survey %in% c("FR-CGFS", "EVHOE")) %>%
#   dplyr::group_by(Lon, Lat, Survey) %>%
#   dplyr::summarise(n_sp = n_distinct(genus_sp))
# df_map_select1 <- df_all_sampled1 %>%
#   dplyr::filter(Survey %in% c("FR-CGFS", "EVHOE")) %>%
#   dplyr::group_by(Lon, Lat, Survey , Quarter) %>%
#   dplyr::summarise(n_sp = n_distinct(genus_sp))
# df_map_select2 <- df_all_sampled2 %>%
#   dplyr::filter(Survey %in% c("FR-CGFS", "EVHOE")) %>%
#   dplyr::group_by(Lon, Lat, Survey , Quarter) %>%
#   dplyr::summarise(n_sp = n_distinct(genus_sp))
# df_map_last <- df_last %>%
#   dplyr::filter(Survey %in% c("FR-CGFS", "EVHOE")) %>%
#   dplyr::group_by(Lon, Lat, Survey , Quarter) %>%
#   dplyr::summarise(n_sp = n_distinct(genus_sp))
# a = ggplot(data = world) + facet_grid(Quarter~ .) + ggtitle("Initial") +
#   geom_point(data = df_map_initial, aes(x = Lon, y = Lat, color = Survey), size = 1) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-21, 10), ylim = c( 40 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.3,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece)),
#                color = 'black') + theme(axis.title = element_blank())
# b = ggplot(data = world) +  facet_grid(Quarter~ .) + ggtitle("Filter 2 trawls per cell") +
#   geom_point(data = df_map_select1, aes(x = Lon, y = Lat, color = Survey), size = 1) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-21, 10), ylim = c( 40 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.3,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece)),
#                color = 'black') + theme(axis.title = element_blank())
# c = ggplot(data = world) +  facet_grid(Quarter~ .) + ggtitle("Plus filter more than 10 years and before 2005") +
#   geom_point(data = df_map_select2, aes(x = Lon, y = Lat, color = Survey), size = 1) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-21, 10), ylim = c( 40 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.3,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece)),
#                color = 'black') + theme(axis.title = element_blank())
# d = ggplot(data = world) +  facet_grid(Quarter~ .) + ggtitle("Plus filter with data env") +
#   geom_point(data = df_map_last, aes(x = Lon, y = Lat, color = Survey), size = 1) +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-21, 10), ylim = c( 40 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.3,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece)),
#                color = 'black') + theme(axis.title = element_blank())
# ggarrange(a, b, c, d, ncol = 4, common.legend = TRUE)
# 
# ggplot(data = world) + facet_wrap(~ Survey, ncol = 6) + 
#   geom_point(data = df_map_initial, aes(x = Lon, y = Lat), size = 1, color = 'red') +
#   geom_sf(color = 'white', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-21, 40), ylim = c( 30 , 62),expand = FALSE) +
#   theme(axis.title = element_blank(), axis.text = element_blank())

# n_distinct(df_last$family )
# n_distinct(df_last$genus )
# n_distinct(df_last$genus_sp )
# 
# n_distinct(df_last$Survey )
# range(df_last$Year )
# 
# n_distinct(df_last$my_spatial_id )
# n_distinct(df_last$index )
# 
# 
# 
# dede <- df_last %>%
#   dplyr::group_by(x_my_spatial_id, y_my_spatial_id , Year, Quarter) %>%
#   dplyr::summarise(n_t = n_distinct(index_old))
# a = ggplot(data = world) +
#   geom_tile(data = dede,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id), color = 'red') +
#   facet_grid(Quarter ~ Year ) +
#   geom_sf(color = NA, fill = 'grey80') +  theme_classic() +
#   coord_sf(xlim = c(-20, 30), ylim = c( 34.5 , 65),expand = FALSE) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(color = NA))
# 
# ggsave(file = 'figures/taxo_fonctio/sampling_map.jpg', plot = a,
#        width = 5,  height = 1, scale = 3)
# 
# dede_initial <- df_total2 %>%
#   dplyr::group_by(x_my_spatial_id, y_my_spatial_id , Year, Quarter) %>%
#   dplyr::summarise(n_t = n_distinct(index))
# b = ggplot(data = world) +
#   geom_tile(data = dede_initial,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id), color = 'red') +
#   facet_grid(Quarter ~ Year ) +
#   geom_sf(color = NA, fill = 'grey80') +  theme_classic() +
#   coord_sf(xlim = c(-20, 30), ylim = c( 34.5 , 65),expand = FALSE) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         strip.background = element_rect(color = NA))
# 
# ggsave(file = 'figures/taxo_fonctio/sampling_map_initial.jpg', plot = b,
#        width = 5,  height = 0.8, scale = 3)
# 
# 
# 
# df_total2_spatial <- df_total2 %>% 
#   dplyr::group_by(my_spatial_id, x_my_spatial_id, y_my_spatial_id, genus_sp) %>% 
#   dplyr::summarise(ab_m = mean(abundance))
# 
# 
# matrice_occurrences <- reshape2::dcast(df_total2_spatial, 
#                                        my_spatial_id +  x_my_spatial_id + y_my_spatial_id ~  genus_sp,
#                                        fun.aggregate = length)
# write.csv2(file = 'data/processed/matrice_occurrences.csv', matrice_occurrences, row.names = FALSE)


################################################################################
#################################### ACP taxo ##################################
################################################################################

trawl_sp_matrice <- reshape2::acast(data = df_last,
                                    index_old ~ genus_sp,
                                    value.var = 'log_abundance',
                                    fun.aggregate = mean)
trawl_sp_matrice[is.na(trawl_sp_matrice)] <- 0
trawl_sp_matrice_hell <- decostand(trawl_sp_matrice, method = 'hellinger')

df_complementary <- df_last %>%
  dplyr::group_by(index_old, Year, Quarter, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>%
  dplyr::summarise(n_sp = n_distinct(genus_sp),
                   abundance_tot = sum(abundance),
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro_mea = mean(chloro_mea, na.rm = TRUE),
                   mlotst_mea = mean(mlotst_mea, na.rm = TRUE),
                   oxy_bottom_mea = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf_mea = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom_mea = mean(temp_bottom_mea, na.rm = TRUE), 
                   temp_surf_mea = mean(temp_surf_mea, na.rm = TRUE),
                   curr_surf_mea   = mean(curr_surf_mea, na.rm = TRUE),
                   curr_bottom_mea = mean(curr_bottom_mea, na.rm = TRUE))

df_complementary_small1 <- data.frame(df_complementary) %>%
  dplyr::select( Year, Quarter, x_my_spatial_id, y_my_spatial_id,
                 depth,  11:19) %>%
  dplyr::mutate(lon_ICES = x_my_spatial_id, 
                lat_ICES = y_my_spatial_id) %>%
  dplyr::select(- x_my_spatial_id, -y_my_spatial_id)

df_complementary_small1_scalled <- scale(df_complementary_small1)


trawl_sp_matrice2 <- cbind(trawl_sp_matrice_hell, df_complementary_small1_scalled)
trawl_sp_matrice2 <- cbind(trawl_sp_matrice2, data.frame(df_complementary[,'Ecoregion']))


names(trawl_sp_matrice2)
n_column_taxo <- length(names(trawl_sp_matrice2))
n_column_taxo-12
n_column_taxo
names(trawl_sp_matrice2)[c((n_column_taxo-14):(n_column_taxo-1))]

res_acp_taxo <- PCA(trawl_sp_matrice2, scale.unit = FALSE, ncp = 50,
                    quanti.sup = c((n_column_taxo-14):(n_column_taxo-1)),
                    quali.sup = n_column_taxo,
                    graph = FALSE)

eig.val_taxo <- get_eigenvalue(res_acp_taxo)
head(eig.val_taxo)

############ results 

df_pca_taxo <-  data.frame(res_acp_taxo$ind$coord)
df_pca_taxo$index_old <- rownames(df_pca_taxo)

df_pca_taxo <- merge(df_complementary, df_pca_taxo)


df_pca_taxo2 <- df_pca_taxo %>% 
  dplyr::group_by(index_old,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   pco3_m  =  mean(Dim.3),
                   pco4_m  =  mean(Dim.4),
                   pco5_m  =  mean(Dim.5))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pco2_m),
                   pcoa3_m  =  mean(pco3_m),
                   pcoa4_m  =  mean(pco4_m),
                   pcoa5_m  =  mean(pco5_m))


################################################################################
#################################### ACP CWMs ##################################
################################################################################
df_genus_sp_traits <- df_last %>%
  dplyr::group_by(genus_sp, habitat, spawning.type, feeding.mode) %>%
  dplyr::summarise(length.maturity = mean((length.maturity)),
                   age.maturity = mean((age.maturity)),
                   growth.coefficient = mean((growth.coefficient)),
                   tl = mean(tl))
# 
# pl_tl <- ggplot(df_genus_sp_traits, aes(x = tl)) + geom_histogram() +
#   theme_classic() + xlab('Trophic level') +
#   theme(axis.title.y = element_blank())
# pl_length <- ggplot(df_genus_sp_traits, aes(x = length.maturity)) + geom_histogram() +
#   theme_classic() + xlab('Length at sexual maturity') +
#   theme(axis.title.y = element_blank())
# pl_age <- ggplot(df_genus_sp_traits, aes(x = age.maturity)) + geom_histogram() +
#   theme_classic() + xlab('Age at sexual maturity') +
#   theme(axis.title.y = element_blank())
# pl_growth <- ggplot(df_genus_sp_traits, aes(x = growth.coefficient)) + geom_histogram() +
#   theme_classic() + xlab('Growth coefficient') +
#   theme(axis.title.y = element_blank())
# 
# pl_spaw <- ggplot(df_genus_sp_traits, aes(x = spawning.type)) + geom_bar() +
#   theme_classic() + xlab('') + theme(axis.title = element_blank())
# pl_hab <- ggplot(df_genus_sp_traits, aes(x = habitat)) + geom_bar() +
#   theme_classic() + xlab('') + theme(axis.title = element_blank())
# pl_diet <- ggplot(df_genus_sp_traits, aes(x = feeding.mode)) + geom_bar() +
#   theme_classic() + xlab('') + theme(axis.title = element_blank())
# 
# # plot
# 
# df_genus_sp_traits <- df_last %>%
#   dplyr::group_by(genus_sp, habitat, spawning.type, feeding.mode) %>%
#   dplyr::summarise(length.maturity = mean(log10(length.maturity)),
#                    age.maturity = mean(log10(age.maturity)),
#                    growth.coefficient = mean(log10(growth.coefficient)),
#                    tl = mean(tl))
# 
# 
# names(df_genus_sp_traits)
# df_traits_alone <- df_genus_sp_traits[,c(2:8)]
# names(df_traits_alone) <- c('Vertical\nhabitat', 'Spawning\ntype', "Feeding\nmode",
#                             "Length at\nmaturity", "Age at\nmaturity", "Growth\ncoefficient",
#                             "Trophic\nlevel")
# corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlGn")[2:6]
# corColors
# 
# df_traits_alone <- data.frame(df_traits_alone)
# my_custom_cor <- function(data, mapping, color = I("black"), sizeRange = c(1, 5), ...) {
# 
#   # get the x and y data to use the other code
#   x <- GGally::eval_data_col(data, mapping$x)
#   y <- GGally::eval_data_col(data, mapping$y)
# 
#   ct <- cor.test(x,y)
# 
#   sig <- symnum(
#     ct$p.value, corr = FALSE, na = FALSE,
#     cutpoints = c(0, 0.001, 0.01, 0.05, 1),
#     symbols = c("***", "**", "*", " ")
#   )
# 
#   r <- unname(ct$estimate)
#   rt <- format(r, digits=2)[1]
#   tt <- as.character(rt)
# 
#   cex <- max(sizeRange)
#   # helper function to calculate a useable size
#   percent_of_range <- function(percent, range) {
#     percent * diff(range) + min(range, na.rm = TRUE)
#   }
# 
#   # plot the cor value
#   p <- ggally_text(
#     label = as.character(rt),
#     mapping = aes(),
#     xP = 0.5, yP = 0.5,
#     size = I(cex),
#     color = color,
#     ...
#   ) +
#     # add the sig stars
#     # geom_text(
#     #   aes_string(
#     #     x = 0.8,
#     #     y = 0.8
#     #   ),
#     #   label = sig,
#     #   size = I(cex),
#     #   color = color,
#     #   ...
#     # ) +
# 
#     theme(panel.background=element_rect(fill="white", color = "black", linetype = "dashed"),
#           panel.grid.minor=element_blank(),
#           panel.grid.major=element_blank())
# 
#   corColors <- RColorBrewer::brewer.pal(n = 7, name = "RdYlGn")[2:6]
# 
#   if (r <= -0.8) {
#     corCol <- corColors[1]
#   } else if (r <= -0.6) {
#     corCol <- corColors[2]
#   } else if (r < 0.6) {
#     corCol <- 'white'
#   } else if (r < 0.8) {
#     corCol <- corColors[4]
#   } else {
#     corCol <- corColors[5]
#   }
#   p <- p + theme(panel.background = element_rect(fill= corCol)
#   )
# 
#   p
# }
# 
# a = ggpairs(df_traits_alone,
#             mapping = ggplot2::aes(alpha = 0.8),
#             upper = list(continuous = my_custom_cor),
#             diag = list(continuous = wrap("densityDiag")),
#             lower = list(continuous = "smooth"))
# 
# ggsave(file = 'figures/taxo_fonctio/S3figure.jpg', plot = a,
#        width = 3.2,  height = 2.4, scale = 3)




########################
df_last$log_length.maturity <- log10(df_last$length.maturity)
df_last$log_age.maturity <- log10(df_last$age.maturity)
df_last$log_growth.coefficient<- log10(df_last$growth.coefficient)
df_last$log_tl <- log10(df_last$tl)

df_CMW_quanti <- df_last %>% 
  dplyr::select(index_old, genus_sp, log_abundance, log_length.maturity, log_age.maturity, 
                log_growth.coefficient, log_tl)

df_CMW_quanti2 <- df_CMW_quanti %>% 
  dplyr::group_by(index_old) %>% 
  dplyr::summarise(length.maturity = weighted.mean(log_length.maturity, log_abundance),
                   age.maturity = weighted.mean(log_age.maturity, log_abundance),
                   growth.coefficient = weighted.mean(log_growth.coefficient, log_abundance),
                   tl = weighted.mean(log_tl, log_abundance))
range(df_last$log_tl)
range(df_CMW_quanti2$tl)

#### spawning.type
df_CMW_quali <- df_last %>% 
  dplyr::select(index_old, genus_sp, log_abundance, 
                spawning.type, habitat, feeding.mode)
df_CMW_quali <- data.frame(df_CMW_quali)
df_CMW_quali_cast_spaw <- reshape2::dcast(df_CMW_quali, 
                                          index_old + genus_sp + log_abundance ~ spawning.type,
                                          length,
                                          value.var = 'spawning.type')

names(df_CMW_quali_cast_spaw)[names(df_CMW_quali_cast_spaw) == "non-guarder" ] <- "non_guarder" 
df_CMW_quali_cast_spaw$bearer <- as.num(df_CMW_quali_cast_spaw$bearer)
df_CMW_quali_cast_spaw$guarder <- as.num(df_CMW_quali_cast_spaw$guarder)
df_CMW_quali_cast_spaw$non_guarder <- as.num(df_CMW_quali_cast_spaw$non_guarder)

df_CMW_quali_cast_spaw2 <- df_CMW_quali_cast_spaw %>% 
  dplyr::group_by(index_old) %>% 
  dplyr::summarise(bearer = weighted.mean(bearer, log_abundance)*100,
                   guarder = weighted.mean(guarder, log_abundance)*100,
                   non_guarder = weighted.mean(non_guarder, log_abundance)*100)

#### habitat

df_CMW_quali_cast_hab <- reshape2::dcast(df_CMW_quali,
                                         index_old + genus_sp + log_abundance ~ habitat,
                                         length,
                                          value.var = "habitat")
names(df_CMW_quali_cast_hab)[names(df_CMW_quali_cast_hab) == "reef-associated" ] <- "reef_associated" 

df_CMW_quali_cast_hab$bathydemersal <- as.num(df_CMW_quali_cast_hab$bathydemersal)
df_CMW_quali_cast_hab$bathypelagic <- as.num(df_CMW_quali_cast_hab$bathypelagic)
df_CMW_quali_cast_hab$benthopelagic <- as.num(df_CMW_quali_cast_hab$benthopelagic)
df_CMW_quali_cast_hab$demersal <- as.num(df_CMW_quali_cast_hab$demersal)
df_CMW_quali_cast_hab$pelagic <- as.num(df_CMW_quali_cast_hab$pelagic)
df_CMW_quali_cast_hab$reef_associated <- as.num(df_CMW_quali_cast_hab$reef_associated)


df_CMW_quali_cast_hab2 <- df_CMW_quali_cast_hab %>% 
  dplyr::group_by(index_old) %>% 
  dplyr::summarise(bathydemersal = weighted.mean(bathydemersal, log_abundance)*100,
                   bathypelagic = weighted.mean(bathypelagic, log_abundance)*100,
                   benthopelagic = weighted.mean(benthopelagic, log_abundance)*100,
                   demersal = weighted.mean(demersal, log_abundance)*100,
                   pelagic = weighted.mean(pelagic, log_abundance)*100,
                   reef_associated = weighted.mean(reef_associated, log_abundance)*100)


#### feeding.mode

df_CMW_quali_cast_feed <- reshape2::dcast(df_CMW_quali,
                                          index_old + genus_sp + log_abundance ~ feeding.mode,
                                          length,
                                          value.var = "feeding.mode")
df_CMW_quali_cast_feed$benthivorous <- as.num(df_CMW_quali_cast_feed$benthivorous)
df_CMW_quali_cast_feed$generalist <- as.num(df_CMW_quali_cast_feed$generalist)
df_CMW_quali_cast_feed$herbivorous <- as.num(df_CMW_quali_cast_feed$herbivorous)
df_CMW_quali_cast_feed$piscivorous <- as.num(df_CMW_quali_cast_feed$piscivorous)
df_CMW_quali_cast_feed$planktivorous <- as.num(df_CMW_quali_cast_feed$planktivorous)

df_CMW_quali_cast_feed2 <- df_CMW_quali_cast_feed %>% 
  dplyr::group_by(index_old) %>% 
  dplyr::summarise(benthivorous = weighted.mean(benthivorous, log_abundance)*100,
                   generalist = weighted.mean(generalist, log_abundance)*100,
                   herbivorous = weighted.mean(herbivorous, log_abundance)*100,
                   piscivorous = weighted.mean(piscivorous, log_abundance)*100,
                   planktivorous = weighted.mean(planktivorous, log_abundance)*100)


#### all together

df_CMW_quali_all <- merge(df_CMW_quali_cast_spaw2, df_CMW_quali_cast_hab2)
df_CMW_quali_all <- merge(df_CMW_quali_all, df_CMW_quali_cast_feed2)

df_CMW_all <- merge(df_CMW_quanti2, df_CMW_quali_all)
df_CMW_all[,c(2:19)] <- scale(df_CMW_all[,c(2:19)])


df_CMW_all2 <- cbind(df_CMW_all, df_complementary_small1_scalled)
df_CMW_all2 <- cbind(df_CMW_all2, data.frame(df_complementary[,'Ecoregion']))


names(df_CMW_all2)
dim(df_CMW_all2)
n_column_taxo <- length(names(trawl_sp_matrice2))
n_column_taxo-13
n_column_taxo

names(df_CMW_all2[,c(2:34)])

result_mfa <- MFA(df_CMW_all2[,c(2:34)], 
                  group = c(1, 1, 1, 1, 
                            3, 6, 5, 
                            2, 1, 9, 2, 1), 
                  type = c("c", "c", "c", "c", "c", "c", "c", 
                           rep("c", 4), 'n' ),
                  name.group = c("length.maturity","age.maturity","growth.coefficient","tl", 
                                 "spawning",
                                 "habitat", 
                                 "diet",
                                 "Temporal","depth",
                                 "env",
                                 "spatial","Ecoregion"),
                  num.group.sup = c(8:12),
                  graph = FALSE)

eig.val <- get_eigenvalue(result_mfa)
head(eig.val)
# fviz_screeplot(result_mfa)


df_pca_fonctio <-  data.frame(result_mfa$ind$coord)
df_pca_fonctio$index_old <- df_CMW_all2$index_old

df_pca_fonctio <- merge(df_complementary, df_pca_fonctio)


df_pca_fonctio2 <- df_pca_fonctio %>% 
  dplyr::group_by(index_old,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   pco3_m  =  mean(Dim.3),
                   pco4_m  =  mean(Dim.4),
                   pco5_m  =  mean(Dim.5))

df_pca_fonctio3 <- df_pca_fonctio2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pco2_m),
                   pcoa3_m  =  mean(pco3_m),
                   pcoa4_m  =  mean(pco4_m),
                   pcoa5_m  =  mean(pco5_m))

# 
# fviz_pca_var(res_acp_taxo, 
#              col.var = "#ff7b7b",
#                     # col.ind = "#dbdbdb" ,
#                     select.var = list(c(contrib = 15)),
#                     col.quanti.sup = "blue", 
#                     # palette = c("#dadada", "#fff785"), 
#                     # habillage = trawl_sp_matrice2$Ecoregion,
#                     alpha.ind = 0.1, 
#                     axes = c(1, 2),
#                     addEllipses = TRUE, # Concentration ellipses
#                     ellipse.type = "convex",
#                     label = c( "quanti.var", "var")) + 
#   theme(legend.position = 'none')
# 
# fviz(result_mfa, 'quanti.var', label = "group.sup")

################################################################################
################################### 2 facets together ##################################

maptaxo <- ggplot(data = world) +  
  geom_tile(data = df_pca_taxo3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
            size = 1) + scale_fill_viridis_c(name = 'PC1', breaks = c(-0.4, 0, 0.4)) + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  geom_text(aes(x = 12, y = 50, label = 'PC1'), color = 'white') + 
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 



maptaxo2 <- ggplot(data = world) +  
  geom_tile(data = df_pca_taxo3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa2_m  ),
            size = 1) + scale_fill_viridis_c(name = 'PC2', breaks  = c(-0.5, -0.2, 0.1, 0.4)) + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  geom_text(aes(x = 12, y = 50, label = 'PC2'), color = 'white') + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 


mapfonctio <- ggplot(data = world) +  
  geom_tile(data = df_pca_fonctio3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
            size = 1) + 
  scale_fill_viridis_c(name = 'MF1',
                       # limits = c(-3.5, 3.3), 
                       na.value = "#F2E52C") + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  geom_text(aes(x = 12, y = 50, label = 'MF1'), color = 'white') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 

mapfonctio2 <- ggplot(data = world) +  
  geom_tile(data = df_pca_fonctio3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa2_m  ),
            size = 1) + 
  scale_fill_viridis_c(name = 'MF2',
                       # limits = c(-3.5, 3.3), 
                       na.value = "#F2E52C") + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  geom_text(aes(x = 12, y = 50, label = 'MF2'), color = 'white') + 
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 


################################################################################
################################### split periods ##################################
# 
# df_pca_taxo2$period <- cut(df_pca_taxo2$Year, breaks = c(1993, 2002, 2011, 2020),  
#                            include.lowest = TRUE)
# table(df_pca_taxo2$period, df_pca_taxo2$Year)
# levels(df_pca_taxo2$period)
# levels(df_pca_taxo2$period) <- c("1994-2002", "2003-2011", "2012-2019")
# df_pca_taxo3_bis  <- df_pca_taxo2 %>% 
#   dplyr::group_by(period, x_my_spatial_id, y_my_spatial_id, my_spatial_id,
#                   Ecoregion) %>% 
#   dplyr::summarise(pcoa1_m = mean(pca1_m),
#                    pcoa2_m  =  mean(pco2_m),
#                    pcoa3_m  =  mean(pco3_m),
#                    pcoa4_m  =  mean(pco4_m),
#                    pcoa5_m  =  mean(pco5_m))
# 
# df_pca_fonctio2$period <- cut(df_pca_fonctio2$Year, breaks = c(1993, 2002, 2011, 2020),  
#                               include.lowest = TRUE)
# table(df_pca_fonctio2$period, df_pca_fonctio2$Year)
# levels(df_pca_fonctio2$period)
# levels(df_pca_fonctio2$period) <- c("1994-2002", "2003-2011", "2012-2019")
# df_pca_fonctio3_bis <- df_pca_fonctio2 %>% 
#   dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
#                   Ecoregion) %>% 
#   dplyr::summarise(period, pcoa1_m = mean(pca1_m),
#                    pcoa2_m  =  mean(pco2_m),
#                    pcoa3_m  =  mean(pco3_m),
#                    pcoa4_m  =  mean(pco4_m),
#                    pcoa5_m  =  mean(pco5_m))
# 
# maptaxo <- ggplot(data = world) +  
#   geom_tile(data = df_pca_taxo3_bis,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
#             size = 1) + 
#   scale_fill_viridis_c(name = 'PC1', breaks = c(-0.4, 0, 0.4)) + 
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   facet_grid(~ period) +
#   geom_text(aes(x = 12, y = 50, label = 'PC1'), color = 'white') +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1)) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = 'bottom',
#         legend.background = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         legend.direction = 'horizontal',
#         plot.title = element_text(face = 'bold')) 
# 
# maptaxo2 <- ggplot(data = world) +  
#   geom_tile(data = df_pca_taxo3_bis,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa2_m  ),
#             size = 1) + scale_fill_viridis_c(name = 'PC2', breaks  = c(-0.5, -0.2, 0.1, 0.4)) + 
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   facet_grid(~ period) +
#   geom_text(aes(x = 12, y = 50, label = 'PC2'), color = 'white') + 
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1)) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = 'bottom',
#         legend.background = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         legend.direction = 'horizontal',
#         plot.title = element_text(face = 'bold')) 
# 
# 
# mapfonctio <- ggplot(data = world) +  
#   geom_tile(data = df_pca_fonctio3_bis,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
#             size = 1) +   facet_grid(~ period) +
#   scale_fill_viridis_c(name = 'MF1') + 
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_text(aes(x = 12, y = 50, label = 'MF1'), color = 'white') +
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1)) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = 'bottom',
#         legend.background = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         legend.direction = 'horizontal',
#         plot.title = element_text(face = 'bold')) 
# 
# mapfonctio2 <- ggplot(data = world) +  
#   geom_tile(data = df_pca_fonctio3_bis,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa2_m  ),
#             size = 1) +   facet_grid(~ period) +
#   scale_fill_viridis_c(name = 'MF2') + 
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_text(aes(x = 12, y = 50, label = 'MF2'), color = 'white') + 
#   guides(fill = guide_colorbar(barwidth = 7, barheight = 1)) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = 'bottom',
#         legend.background = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         legend.direction = 'horizontal',
#         plot.title = element_text(face = 'bold')) 
# 
# maps_periods <- ggarrange(maptaxo, maptaxo2, 
#                           mapfonctio, mapfonctio2,
#                           ncol = 2, nrow = 2)
# # maps_kmean
# 
# ggsave(file = 'figures/taxo_fonctio/figureSX_maps_periods.jpg', 
#        plot = maps_periods,
#        width = 3.5,  height = 1.7, scale = 3)
################################################################################
################################### classification ##################################

# res_nbclus_taxo <-NbClust(df_pca_taxo3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")],
#                           method = 'kmeans')
# res_nbclus_fonctio <-NbClust(df_pca_fonctio3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")],
#                           method = 'kmeans')
# 
# dist_taxo <- dist(df_pca_taxo3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")])
# res_hclust_taxo = hclust(dist_taxo, method = "complete")
# df_pca_taxo3$kmeans <- as.factor(cutree(res_hclust_taxo, k = 5))
# 
# dist_fonctio <- dist(df_pca_fonctio3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")])
# res_hclust_foncti = hclust(dist_fonctio, method = "complete")
# df_pca_fonctio3$kmeans <- as.factor(cutree(res_hclust_foncti, k = 5))
# 
# 
# df_pca_taxo_kk <- df_pca_taxo3 %>% 
#   dplyr::mutate(pc1_t = pcoa1_m, 
#                 pc2_t = pcoa2_m, 
#                 pc3_t = pcoa3_m) %>% 
#   dplyr::select(x_my_spatial_id, y_my_spatial_id, my_spatial_id , pc1_t, pc2_t, pc3_t)
# df_pca_fonctio_kk <- df_pca_fonctio3 %>% 
#   dplyr::mutate(pc1_f = pcoa1_m, 
#                 pc2_f = pcoa2_m, 
#                 pc3_f = pcoa3_m) %>% 
#   dplyr::select( my_spatial_id , pc1_f, pc2_f, pc3_f)
# df_pca_both <- merge(df_pca_taxo_kk, df_pca_fonctio_kk)
# 
# 
# dist_both <- dist(df_pca_both[,c("pc1_t", "pc2_t", "pc3_t", 
#                                     "pc1_f", "pc2_f", "pc3_f")])
# res_hclust_both = hclust(dist_both, method = "complete")
# df_pca_both$kmeans <- as.factor(cutree(res_hclust_both, k = 5))
# 
# maptaxo_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_taxo3,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   scale_fill_manual(values = c("#FF7175", "#FFCE4E","#8CDB5E",
#                                "#28618A",  "#E353D1")) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   ggtitle("Taxonomic") +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
# 
# mapfonctio_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_fonctio3,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   ggtitle("Functional") +
#   scale_fill_manual(values = c("#FF7175",  "#FFCE4E",   "#8CDB5E", 
#                                "#28618A", "#E353D1")) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
# 
# mapboth_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_both,
#             aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   ggtitle("Both") +
#   scale_fill_manual(values = c("#8CDB5E", "#E353D1", "#FFCE4E",
#                                "#FF7175", "#28618A")) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
# 
# maps_kmean <- ggarrange(maptaxo_kmeans, mapfonctio_kmeans, mapboth_kmeans, ncol = 3)
# # maps_kmean
# 
# ggsave(file = 'figures/taxo_fonctio/figureSX_kmeans.jpg', plot = maps_kmean,
#        width = 3,  height = 1, scale = 3)
# 
# 
# dist_taxo <- dist(df_pca_taxo3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")])
# res_hclust_taxo = hclust(dist_taxo, method = "complete")
# df_pca_taxo3$kmeans <- as.factor(cutree(res_hclust_taxo, k = 10))
# 
# dist_fonctio <- dist(df_pca_fonctio3[,c("pcoa1_m", "pcoa2_m", "pcoa3_m")])
# res_hclust_foncti = hclust(dist_fonctio, method = "complete")
# df_pca_fonctio3$kmeans <- as.factor(cutree(res_hclust_foncti, k = 10))
# 
# 
# df_pca_taxo_kk <- df_pca_taxo3 %>%
#   dplyr::mutate(pc1_t = pcoa1_m,
#                 pc2_t = pcoa2_m,
#                 pc3_t = pcoa3_m) %>%
#   dplyr::select(x_my_spatial_id, y_my_spatial_id, my_spatial_id , pc1_t, pc2_t, pc3_t)
# df_pca_fonctio_kk <- df_pca_fonctio3 %>%
#   dplyr::mutate(pc1_f = pcoa1_m,
#                 pc2_f = pcoa2_m,
#                 pc3_f = pcoa3_m) %>%
#   dplyr::select( my_spatial_id , pc1_f, pc2_f, pc3_f)
# df_pca_both <- merge(df_pca_taxo_kk, df_pca_fonctio_kk)
# 
# 
# dist_both <- dist(df_pca_both[,c("pc1_t", "pc2_t", "pc3_t",
#                                  "pc1_f", "pc2_f", "pc3_f")])
# res_hclust_both = hclust(dist_both, method = "complete")
# df_pca_both$kmeans <- as.factor(cutree(res_hclust_both, k = 10))
# 
# maptaxo_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_taxo3,aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   # scale_fill_manual(values = c("#FF7175", "#FFCE4E","#8CDB5E",
#   #                              "#28618A",  "#E353D1")) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   ggtitle("Taxonomic") +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
# 
# mapfonctio_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_fonctio3,aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   ggtitle("Functional") +
#   # scale_fill_manual(values = c("#FF7175",  "#FFCE4E",   "#8CDB5E",
#   #                              "#28618A", "#E353D1")) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
# 
# mapboth_kmeans <- ggplot(data = world) +
#   geom_tile(data = df_pca_both,aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = kmeans)) +
#   ggtitle("Both") +
#   # scale_fill_manual(values = c("#8CDB5E", "#E353D1", "#FFCE4E",
#   #                              "#FF7175", "#28618A")) +
#   geom_sf(color = 'grey80', fill = 'grey50', size = 0.05) +  theme_classic() +
#   coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
#   geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = 'grey20', alpha = 0.5,
#                aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
#   theme(axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.line = element_blank(),
#         legend.spacing.x = unit(0.05, 'cm'),
#         axis.ticks = element_blank(),
#         legend.title =  element_blank(),
#         legend.position = "none",
#         legend.background = element_blank(),
#         legend.direction = 'vertical',
#         plot.title = element_text(face = 'bold'))
######################################################################
################################### spatial patterns
df_pca_taxo$Ecoregion <- as.factor(df_pca_taxo$Ecoregion)
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_pca_fonctio$Ecoregion <- as.factor(df_pca_fonctio$Ecoregion)
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_text_taxo <- df_pca_taxo %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2))
df_text_taxo$group <- c('l','l',
                        't','t',
                        'r', 'r',
                        'b','b')


param_lineheit <- 0.6
levels(df_pca_taxo$Ecoregion)  
df_pca_taxo$Ecoregion <- factor(df_pca_taxo$Ecoregion,
                                           levels = c("Baltic Sea" , 
                                                      "Greater\nNorth Sea", 
                                                      "Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))

plt_space_taxo = ggplot(df_pca_taxo, aes(x = Dim.1, y = Dim.2, color = Ecoregion ))  +
  geom_point(size = 0.2, alpha = 0.03) +   scale_color_manual(values = eight_cols) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  stat_ellipse(aes(fill = Ecoregion),
               geom="polygon",level=0.8,alpha=0.1) +
  geom_point(data = df_text_taxo, aes(x = dim1_m, y = dim2_m,
                                      color = Ecoregion),
             size = 3, shape = 17) +  xlab('PC1 (21%)') + ylab('PC2 (12%)') +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group == 'r', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_x = 0.7, direction = "y", hjust = "left",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group == 't', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = 0.4, direction = "x",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group ==  'l', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_x = -0.7, direction = "y",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group ==  'b', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = -0.6, direction = "x",  
                   lineheight = param_lineheit) +
  coord_cartesian(xlim = c(-1.3, 1.85), ylim = c(-0.8, 0.8)) + theme_classic() + 
  theme(legend.position = 'none') + 
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))


levels(df_pca_fonctio$Ecoregion)  
df_pca_fonctio$Ecoregion <- factor(df_pca_fonctio$Ecoregion,
                                levels = c("Baltic Sea" , 
                                           "Greater\nNorth Sea", 
                                           "Celtic Seas", 
                                           "Bay of Biscay and\nIberian Coast",
                                           "Western\nMed. Sea",
                                           "Ionian Sea and\nCentral Med. Sea",
                                           "Adriatic Sea"    ,
                                           "Aegean-Levantine\nSea"))

df_text_fonctio <- df_pca_fonctio %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2))
df_text_fonctio$group <- c('l','b','r','t2',
                           't','b','l','t')


plt_space_functio = ggplot(df_pca_fonctio, aes(x = Dim.1, y = Dim.2, color =Ecoregion ))  +
  geom_point(size = 0.2, alpha = 0.03) +  scale_color_manual(values = eight_cols) +
  stat_ellipse(aes(fill= Ecoregion),
               geom="polygon",level=0.8,alpha=0.1) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(data = df_text_fonctio, aes(x = dim1_m, y = dim2_m,
                                         color = Ecoregion),
             size = 3, shape = 17) +xlab('MF1 (27%)') + ylab('MF2 (16%)') +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'r', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_x = 3.5, direction = "y", hjust = "left",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'l', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_x = -4, direction = "y",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 't2', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = 5, direction = "x",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 't', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = 3.5, direction = "x",  
                   lineheight = param_lineheit) +

  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'b', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = -4, direction = "x",  
                   lineheight = param_lineheit) +
  coord_cartesian(xlim = c(-6.5, 7), ylim = c(-5, 6.5)) + theme_classic() + 
  theme(legend.position = 'none') + 
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))

################################################################################
################################### espaces ##################################
################################################################################

################################### espaces

coord_sp <- as.data.frame(res_acp_taxo$var$coord)
coord_sp$sp <- rownames(coord_sp)
coord_sp$sp2 <- coord_sp$sp

# coord_sp$sp2 <- paste0('italic(', coord_sp$sp2 , ')')
coord_sp$sp2 <- gsub( "_", '\n', coord_sp$sp2)


quant_pc1_taxo1 <- quantile(coord_sp$Dim.1, 0.005)
quant_pc1_taxo2 <- quantile(coord_sp$Dim.1, 0.995)

quant_pc2_taxo1 <- quantile(coord_sp$Dim.2, 0.005)
quant_pc2_taxo2 <- quantile(coord_sp$Dim.2, 0.995)

coord_sp_sm <- coord_sp %>% 
  dplyr::filter((Dim.2 <= quant_pc2_taxo1 | Dim.2 >= quant_pc2_taxo2) |
                  (Dim.1 <= quant_pc1_taxo1| Dim.1 >= quant_pc1_taxo2))


# img_so <- pick_phylopic(name = "Merluccius merluccius", n = 5)
img_limande <- get_phylopic("b1d9363d-7ba2-4eda-9998-f45cd6d703f1")
img_sprat <- get_phylopic("381e7e19-bb44-4432-85ce-627f05d159a3")
img_merlu <- get_phylopic("9d48fc96-a097-4a14-890d-02e386a16894")
img_melano <- get_phylopic("e2eb38d0-cc5c-473c-a944-4d4dde610344")


coord_sp_sm[coord_sp_sm$sp2 == "Clupea\nharengus", "sp2"] <- "Clup.\nharengus"
coord_sp_sm[coord_sp_sm$sp2 == "Hippoglossoides\nplatessoides", "sp2"] <- "Hippo.\nplatessoides"
coord_sp_sm[coord_sp_sm$sp2 == "Limanda\nlimanda", "sp2"] <- "Lim.\nlimanda"
coord_sp_sm[coord_sp_sm$sp2 == "Melanogrammus\naeglefinus", "sp2"] <- "Melano.\naeglefinus"
coord_sp_sm[coord_sp_sm$sp2 == "Merluccius\nmerluccius","sp2" ] <- "Merlu.\nmerluccius"
coord_sp_sm[coord_sp_sm$sp2 == "Micromesistius\npoutassou","sp2" ] <- "Micro.\npoutassou"
coord_sp_sm[coord_sp_sm$sp2 == "Platichthys\nflesus", "sp2"] <- "Plat.\nflesus"
coord_sp_sm[coord_sp_sm$sp2 == "Pleuronectes\nplatessa", "sp2"] <- "Pleuro.\nplatessa"
coord_sp_sm[coord_sp_sm$sp2 == "Sprattus\nsprattus", "sp2"] <- "Sprat.\nsprattus"
coord_sp_sm[coord_sp_sm$sp2 == "Trachurus\ntrachurus", "sp2"] <- "Trach.\ntrachurus"
coord_sp_sm[coord_sp_sm$sp2 == "Trisopterus\nesmarkii","sp2" ] <- "Triso.\nesmarkii"
coord_sp_sm[coord_sp_sm$sp2 == "Merlangius\nmerlangus", "sp2"] <- "Merlan.\nmerlangus"

plt_sp_dim1_dim2 <- ggplot(coord_sp_sm, aes(x = Dim.1, y = Dim.2)) + 
  geom_point(data = coord_sp,  aes(x = Dim.1, y = Dim.2),
             color= 'red', size = 1, alpha = 0.1) + 
  
  # add_phylopic(img = img_limande, x = -0.165, y = 0.01, ysize = 0.008) +
  # add_phylopic(img = img_sprat, x = -0.1, y = -0.11, ysize = 0.008) +
  # add_phylopic(img = img_merlu, x = 0.1, y = 0.035, ysize = 0.008)  +
  # add_phylopic(img = img_melano, x = -0.05, y = 0.12, ysize = 0.008) +
  geom_hline(yintercept = 0, color = 'grey70') + 
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(color = 'red') + 
  geom_text_repel(aes(label = sp2), force = 5, max.overlaps = Inf,  
                  lineheight = param_lineheit,  fontface = "italic")  +
  ylab("PC2 (12%)") +  xlab("PC1 (21%)") + theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 


plt_sp_dim1_dim2

# quant_pc3_taxo1 <- quantile(coord_sp$Dim.3, 0.005)
# quant_pc3_taxo2 <- quantile(coord_sp$Dim.3, 0.995)
# 
# coord_sp_sm2 <- coord_sp %>% 
#   dplyr::filter((Dim.3 <= quant_pc3_taxo1 | Dim.3 >= quant_pc3_taxo2) |
#                   (Dim.1 <= quant_pc1_taxo1| Dim.1 >= quant_pc1_taxo2))
# plt_sp_dim1_dim3 <- ggplot(coord_sp_sm2, 
#                            aes(x = Dim.1, y = Dim.3)) + 
#   geom_hline(yintercept = 0, color = 'grey70') + 
#   geom_vline(xintercept = 0, color = 'grey70') +
#   geom_point(color = 'red') + 
#   geom_text_repel(aes(label = sp2), force = 3)  +
#   ylab("PC3") +  xlab("PC1 (21%)") + theme_classic()   +
#   theme(panel.background = element_rect(fill = 'white', color = 'black'),
#         panel.grid.major = element_line(color = 'grey90'))
# 


## en haut a gauche :  high growth coeff and plankti
df_genus_sp_traits %>% 
  dplyr::filter(feeding.mode == "planktivorous" &
                  growth.coefficient >= 1 )

# img_gob <- pick_phylopic(name = "Gobiidae", n = 5)
img_gobie <- get_phylopic("6feced70-f5c2-498d-a513-187dda3393b3")

## en haut a droite : generalis, high TL

quantile(df_genus_sp_traits$length.maturity)
df_genus_sp_traits %>% 
  dplyr::filter(feeding.mode == "generalist"  &
                  length.maturity >= 100 )
# img_re <- pick_phylopic(name = "Hexanchus Griseus", n = 5)
img_requin <- get_phylopic("a14db634-f3d1-437a-9d0b-21c03bfd0696")

## en haut au milieu : generalist, benthopelagic

quantile(df_genus_sp_traits$length.maturity)
# View(df_genus_sp_traits %>%
  # dplyr::filter(feeding.mode == "generalist"  &
                  # habitat == "benthopelagic" ))
# img_ze <- pick_phylopic(name = "Zenopsis nebulosa", n = 5)
img_zenop <- get_phylopic("b554ef96-1672-4ba9-83d7-028f08d542e1")


## en bas  : demersal benthivorus
# View(df_genus_sp_traits %>% 
#   dplyr::filter(feeding.mode == "benthivorous"  &
#                   habitat == "demersal" & 
#                   spawning.type == "non-guarder"))

# img_so <- pick_phylopic(name = "Solea solea", n = 5)
img_sole <- get_phylopic("5dae2b2d-412d-4da0-a5e7-52253a5c7c92")
img_sole <- rotate_phylopic(img = img_sole, angle = -90)
plot(img_sole)

## figure

coord_traits <- as.data.frame(result_mfa$quanti.var$coord)
coord_traits$traits <- rownames(coord_traits)

quant_pc1_f1 <- quantile(coord_traits$Dim.1, 0.1)
quant_pc1_f2 <- quantile(coord_traits$Dim.1, 0.9)

quant_pc2_f1 <- quantile(coord_traits$Dim.2, 0.1)
quant_pc2_f2 <- quantile(coord_traits$Dim.2, 0.9)

coord_traits_sm <- coord_traits %>% 
  dplyr::filter((Dim.2 <= quant_pc2_f1 | Dim.2 >= quant_pc2_f2) |
                  (Dim.1 <= quant_pc1_f1| Dim.1 >= quant_pc1_f2))
coord_traits_sm$traits[coord_traits_sm$traits == "length.maturity"] <- 'Maturity\nlength'
coord_traits_sm$traits[coord_traits_sm$traits == "growth.coefficient"] <- 'Growth\ncoefficient'
coord_traits_sm$traits[coord_traits_sm$traits == "tl"] <- 'Trophic\nlevel'

coord_traits_sm$traits[coord_traits_sm$traits == "benthopelagic"] <- 'Benthopelagic'
coord_traits_sm$traits[coord_traits_sm$traits == "demersal"] <- 'Demersal'
coord_traits_sm$traits[coord_traits_sm$traits == "benthivorous"] <- 'Benthivorous'
coord_traits_sm$traits[coord_traits_sm$traits == "generalist"] <- 'Generalist'
coord_traits_sm$traits[coord_traits_sm$traits == "planktivorous"] <- 'Planktivorous'

plt_trai_dim1_dim2 <- ggplot(coord_traits_sm,
                        aes(x = Dim.1, y = Dim.2)) + 
  geom_hline(yintercept = 0, color = 'grey70') + 
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(color = 'red') +  
  geom_text_repel(aes(label = traits),force = 3,  
                  lineheight = param_lineheit,  
                  nudge_x = ifelse(coord_traits_sm$traits %in% c("Planktivorous" , "Demersal"),
                                   0.5, ifelse(coord_traits_sm$traits %in% c("Maturity\nlength" ),
                                               -0.6, 0)),
                  nudge_y = ifelse(coord_traits_sm$traits %in% c("Benthivorous", "Demersal"),
                                   0.2, -0.3))  +
  ylab("MF2 (16%)") +  xlab("MF1 (27%)") +  theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 
  # add_phylopic(img = img_sole, x = -0.1, y = -0.9, ysize = 0.08) +
  # add_phylopic(img = img_requin, x = 0.67, y = -0.2, ysize = 0.08) +
  # add_phylopic(img = img_gobie, x = -0.7, y = -0.2, ysize = 0.08) +
  # add_phylopic(img = img_zenop, x = 0, y = 0.5, ysize = 0.15) 




# quant_pc3_f1 <- quantile(coord_traits$Dim.3, 0.1)
# quant_pc3_f2 <- quantile(coord_traits$Dim.3, 0.9)
# coord_traits_sm2 <- coord_traits %>% 
#   dplyr::filter((Dim.3 <= quant_pc3_f1 | Dim.3 >= quant_pc3_f2) |
#                   (Dim.1 <= quant_pc1_f1| Dim.1 >= quant_pc1_f2))
# plt_trai_dim1_dim3 <- ggplot(coord_traits_sm2,
#                              aes(x = Dim.1, y = Dim.3)) + 
#   geom_hline(yintercept = 0, color = 'grey70') +  
#   geom_vline(xintercept = 0, color = 'grey70') +
#   geom_point(color = 'red') +   
#   geom_text_repel(aes(label = traits),force = 3)  +
#   ylab("PC3") +  xlab("PC1") +  theme_classic()   +
#   theme(panel.background = element_rect(fill = 'white', color = 'black'),
#         panel.grid.major = element_line(color = 'grey90'))



maps_taxo = ggarrange(maptaxo, maptaxo2,  nrow = 2, align = 'hv',
               labels = c("C", 'D'),  
               font.label = list(size = 12, face = "plain"))

fig1_top  = ggarrange(plt_sp_dim1_dim2,plt_space_taxo, maps_taxo,  ncol = 3, nrow = 1, align = 'hv',
               widths  = c(1.2, 1.2, 0.8),
               labels = c("A", 'B', ''),  font.label = list(size = 12, face = "plain"))


maps_function = ggarrange(mapfonctio, mapfonctio2, nrow = 2, align = 'hv',
               labels = c('G', 'H'),  font.label = list(size = 12, face = "plain"))

fig1_bottom = ggarrange(plt_trai_dim1_dim2,plt_space_functio, maps_function,  ncol = 3, nrow = 1, align = 'hv',
                        widths  = c(1.2, 1.2, 0.8),
                        labels = c("E", 'F', ''),  font.label = list(size = 12, face = "plain"))


fig1_all = ggarrange(fig1_top, fig1_bottom, nrow = 2) 
 
fig1_all2 = annotate_figure(fig1_all, 
                            right = textGrob("Taxonomic", rot=-90,
                                               # vjust=-1, 
                                             hjust = 2.1,
                                        gp = gpar(cex = 1.3)))

fig1_all3 = annotate_figure(fig1_all2, 
                            right = textGrob("Functional", rot=-90,
                                             vjust=2.2,
                                             hjust = -1,
                                             gp = gpar(cex = 1.3)))

ggsave(file = 'figures/taxo_fonctio/figure1.jpg', plot = fig1_all3,
       width = 3,  height = 2, units = c("in"))


################################################################################
################################### 2 facets together ##################################

maptaxo3 <- ggplot(data = world) +  
  geom_tile(data = df_pca_taxo3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa3_m  ),
            size = 1) + scale_fill_viridis_c(name = 'PC3') + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  geom_text(aes(x = 12, y = 50, label = 'PC3'), color = 'white') + 
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 

mapfonctio3 <- ggplot(data = world) +  
  geom_tile(data = df_pca_fonctio3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa3_m  ),
            size = 1) + 
  scale_fill_viridis_c(name = 'MF3',
                       # limits = c(-3.5, 3.3), 
                       na.value = "#F2E52C") + 
  geom_sf(color = 'grey80', fill = 'grey50', size = 0.1) +  theme_classic() +
  coord_sf(xlim = c(-18, 31), ylim = c( 34.5 , 62),expand = FALSE) +
  guides(fill = guide_colorbar(barheight = 4, barwidth = 0.4)) +
  geom_text(aes(x = 12, y = 50, label = 'MF3'), color = 'white') +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.title =  element_blank(),
        legend.position = c(-0.08, 0.55),
        legend.background = element_blank(),
        legend.spacing.x = unit(0.05, 'cm'),
        legend.direction = 'vertical',
        plot.title = element_text(face = 'bold')) 

################################################################################
################################### spatial patterns
df_pca_taxo$Ecoregion <- as.factor(df_pca_taxo$Ecoregion)
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 
df_pca_taxo$Ecoregion <- factor(df_pca_taxo$Ecoregion,
                                   levels = c("Baltic Sea" , 
                                              "Greater\nNorth Sea", 
                                              "Celtic Seas", 
                                              "Bay of Biscay and\nIberian Coast",
                                              "Western\nMed. Sea",
                                              "Ionian Sea and\nCentral Med. Sea",
                                              "Adriatic Sea"    ,
                                              "Aegean-Levantine\nSea"))



df_text_taxo <- df_pca_taxo %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim3_m = mean(Dim.3))
df_text_taxo$group <- c('l','l',
                        't','t',
                        'r', 'r',
                        'b','b')
param_lineheit <- 0.6

plt_space_taxo2 = ggplot(df_pca_taxo, aes(x = Dim.1, y = Dim.3, color = Ecoregion ))  +
  geom_point(size = 0.2, alpha = 0.03) + scale_color_manual(values = eight_cols) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  stat_ellipse(aes(fill = Ecoregion),
               geom="polygon",level=0.8,alpha=0.1) +
  geom_point(data = df_text_taxo, aes(x = dim1_m, y = dim3_m,
                                      color = Ecoregion),
             size = 3, shape = 17) +  xlab('PC1 (21%)') + ylab('PC3 (9%)') +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group == 'r', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_x = 0.8, direction = "y", hjust = "left",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group == 't', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_y = 0.6, direction = "x",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group ==  'l', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_x = -1, direction = "y",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_taxo[df_text_taxo$group == 'b', ],
                   aes(x = dim1_m, y = dim2_m,    label = Ecoregion),
                   nudge_y = -4, direction = "x",  
                   lineheight = param_lineheit) +
  
  coord_cartesian(xlim = c(-2, 2.3), ylim = c(-0.8, 0.8)) + theme_classic() + 
  theme(legend.position = 'none') + 
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))

df_pca_fonctio$Ecoregion <- as.factor(df_pca_fonctio$Ecoregion)
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_pca_fonctio$Ecoregion <- factor(df_pca_fonctio$Ecoregion,
                                   levels = c("Baltic Sea" , 
                                              "Greater\nNorth Sea", 
                                              "Celtic Seas", 
                                              "Bay of Biscay and\nIberian Coast",
                                              "Western\nMed. Sea",
                                              "Ionian Sea and\nCentral Med. Sea",
                                              "Adriatic Sea"    ,
                                              "Aegean-Levantine\nSea"))

df_text_fonctio <- df_pca_fonctio %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim3_m = mean(Dim.3))
df_text_fonctio$group <- c('l','b','b','r',
                           't','t','l','r')

plt_space_functio2 = ggplot(df_pca_fonctio, aes(x = Dim.1, y = Dim.2, color =Ecoregion ))  +
  geom_point(size = 0.2, alpha = 0.03)+  scale_color_manual(values = eight_cols) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  stat_ellipse(aes(fill= Ecoregion),
               geom="polygon",level=0.8,alpha=0.1) +
  geom_point(data = df_text_fonctio, aes(x = dim1_m, y = dim3_m,
                                         color = Ecoregion),
             size = 3, shape = 17) + xlab('MF1 (27%)') + ylab('MF3 (11%)') +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'r', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_x = 3.5, direction = "y", hjust = "left",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'l', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_x = -4, direction = "y",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 't', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_y = 4, direction = "x",  
                   lineheight = param_lineheit) +
  geom_label_repel(data = df_text_fonctio[df_text_fonctio$group == 'b', ],
                   aes(x = dim1_m, y = dim3_m,    label = Ecoregion),
                   nudge_y = -4, direction = "x",  
                   lineheight = param_lineheit) +
  coord_cartesian(xlim = c(-7, 10), ylim = c(-6, 6)) + theme_classic() + 
  theme(legend.position = 'none') + 
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))



################################################################################
################################### espaces ##################################

coord_sp <- as.data.frame(res_acp_taxo$var$coord)
coord_sp$sp <- rownames(coord_sp)
coord_sp$sp2 <- coord_sp$sp
coord_sp$sp2 <- gsub( "_", '\n', coord_sp$sp2)


quant_pc1_taxo1 <- quantile(coord_sp$Dim.1, 0.005)
quant_pc1_taxo2 <- quantile(coord_sp$Dim.1, 0.995)

quant_pc3_taxo1 <- quantile(coord_sp$Dim.3, 0.005)
quant_pc3_taxo2 <- quantile(coord_sp$Dim.3, 0.995)

coord_sp_sm2 <- coord_sp %>%
  dplyr::filter((Dim.3 <= quant_pc3_taxo1 | Dim.3 >= quant_pc3_taxo2) |
                  (Dim.1 <= quant_pc1_taxo1| Dim.1 >= quant_pc1_taxo2))
plt_sp_dim1_dim3 <- ggplot(coord_sp_sm2,
                           aes(x = Dim.1, y = Dim.3)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(color = 'red') +
  geom_text_repel(aes(label = sp2), force = 5, max.overlaps = Inf,  
                  lineheight = param_lineheit,  fontface = "italic")  +
  ylab("PC3 (9%)") +  xlab("PC1 (21%)") + theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))


coord_traits <- as.data.frame(result_mfa$quanti.var$coord)
coord_traits$traits <- rownames(coord_traits)

quant_pc1_f1 <- quantile(coord_traits$Dim.1, 0.1)
quant_pc1_f2 <- quantile(coord_traits$Dim.1, 0.9)

quant_pc3_f1 <- quantile(coord_traits$Dim.3, 0.1)
quant_pc3_f2 <- quantile(coord_traits$Dim.3, 0.9)
coord_traits_sm2 <- coord_traits %>%
  dplyr::filter((Dim.3 <= quant_pc3_f1 | Dim.3 >= quant_pc3_f2) |
                  (Dim.1 <= quant_pc1_f1| Dim.1 >= quant_pc1_f2))
plt_trai_dim1_dim3 <- ggplot(coord_traits_sm2,
                             aes(x = Dim.1, y = Dim.3)) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(color = 'red') +
  geom_text_repel(aes(label = traits),force = 3)  +
  xlab('MF1 (27%)') + ylab('MF3 (11%)') +  theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90'))




fig1_top  = ggarrange(plt_sp_dim1_dim3,plt_space_taxo2, maptaxo3,  ncol = 3, nrow = 1, align = 'hv',
                      widths  = c(1.2, 1.2, 0.9),
                      labels = c("A", 'B', 'C'),  font.label = list(size = 12, face = "plain"))

fig1_bottom = ggarrange(plt_trai_dim1_dim3,plt_space_functio2, mapfonctio3,  ncol = 3, nrow = 1, align = 'hv',
                        widths  = c(1.2, 1.2, 0.9),
                        labels = c("D", 'E', 'F'),  font.label = list(size = 12, face = "plain"))


fig1_all = ggarrange(fig1_top, fig1_bottom, nrow = 2)

ggsave(file = 'figures/taxo_fonctio/figureS6.jpg', plot = fig1_all,
       width = 3.5,  height = 2.1, scale = 3)



################################################################################
################################### comparison taxo fonctio

head(df_pca_fonctio)
head(df_pca_taxo)

df_pca_fonctio_melt <- reshape2::melt(data = df_pca_fonctio, 
                                      id.vars = c('index_old', 'Year', 'Quarter', 
                                                  "x_my_spatial_id", "y_my_spatial_id", "my_spatial_id",
                                                  'depth'),
                                      measure.vars = c('Dim.1','Dim.2'),
                                      value.name = 'Functionnal')

df_pca_taxo_melt <- reshape2::melt(data = df_pca_taxo, 
                                      id.vars = c('index_old'),
                                      measure.vars = c('Dim.1','Dim.2'),
                                   value.name = 'Taxonomy')
df_pca_both <- merge(df_pca_taxo_melt, df_pca_fonctio_melt)
levels(df_pca_both$variable) <- c("PC1", "PC2", "PC3", "PC4", "PC5")


ee <- ggplot(df_pca_both, aes(x = Taxonomy, y = Functionnal)) + 
  geom_point(aes(color = variable), alpha = 0.2, size = 0.5) + 
  # geom_smooth() + 
  guides(color = guide_legend(override.aes = list(size = 1, alpha = 1))) +
  ylim(-9, 7) +  theme_classic() +
  theme(panel.grid.major = element_line(color = 'grey90'),
        legend.title = element_blank())
# ee
ggsave(file = 'figures/taxo_fonctio/figure_taxo_functio.jpg', plot = ee,
       width = 2.5,  height = 1.6, scale = 3)



################################################################################
################ influence external variables ##################################
################################################################################

var_quanti_taxo <- res_acp_taxo$quanti.sup$cor[,1:5]
df_var_quanti_taxo <- reshape2::melt(var_quanti_taxo)
df_var_quanti_taxo$value <- round(df_var_quanti_taxo$value, 1)
df_var_quanti_taxo <- df_var_quanti_taxo %>% 
  dplyr::filter(Var1 %ni% c('curr_surf_mea', 'oxy_surf_mea', 'curr_bottom_mea'))
df_var_quanti_taxo$Var1 <- as.factor(df_var_quanti_taxo$Var1)
df_var_quanti_taxo$Var1 <- factor(df_var_quanti_taxo$Var1, exclude = T)
levels(df_var_quanti_taxo$Var1)
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'depth'] <- 'Depth'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'chloro_mea'] <- 'Surface chloro'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'mlotst_mea'] <- 'Mixed layer'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'oxy_bottom_mea'] <- 'Bottom oxygen'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'temp_bottom_mea'] <- 'Bottom temp'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'temp_surf_mea'] <- 'Surface temp'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'lon_ICES'] <- 'Longitude'
levels(df_var_quanti_taxo$Var1)[levels(df_var_quanti_taxo$Var1) == 'lat_ICES'] <- 'Latitude'


plot_var_tax <- ggplot(data = df_var_quanti_taxo, 
                        aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = '#2e79a3', mid = "white", high = '#d44059',
                       midpoint = 0,
                       limit = c(-0.91,0.5),
                       name = 'Pearson\nCorrelation') +
  theme_minimal()+ xlab('') + ylab('') +
  geom_text(data = df_var_quanti_taxo,
            aes(Var1, Var2, label = value), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

var_quanti_fonctio <- result_mfa$quanti.var.sup$cor
df_var_quanti_fonc <- reshape2::melt(var_quanti_fonctio[,1:5])
df_var_quanti_fonc$value <- round(df_var_quanti_fonc$value, 1)
df_var_quanti_fonc <- df_var_quanti_fonc %>% 
  dplyr::filter(Var1 %ni% c('curr_surf_mea', 'oxy_surf_mea', 'curr_bottom_mea'))
df_var_quanti_fonc$Var1 <- as.factor(df_var_quanti_fonc$Var1)
df_var_quanti_fonc$Var1 <- factor(df_var_quanti_fonc$Var1, exclude = T)
levels(df_var_quanti_fonc$Var1)
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'depth'] <- 'Depth'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'chloro_mea'] <- 'Surface chloro'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'mlotst_mea'] <- 'Mixed layer'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'oxy_bottom_mea'] <- 'Bottom oxygen'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'temp_bottom_mea'] <- 'Bottom temp'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'temp_surf_mea'] <- 'Surface temp'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'lon_ICES'] <- 'Longitude'
levels(df_var_quanti_fonc$Var1)[levels(df_var_quanti_fonc$Var1) == 'lat_ICES'] <- 'Latitude'

plot_var_fon <- ggplot(data = df_var_quanti_fonc, 
                        aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = '#2e79a3', mid = "white", high = '#d44059',
                       midpoint = 0,
                       name = 'Pearson\nCorrelation') +
  theme_minimal()+ xlab('') + ylab('') +
  geom_text(aes(Var1, Var2, label = value), color = "black", size = 4) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

plot_var = ggarrange(plot_var_tax, plot_var_fon, nrow = 2, 
                    labels = c("(A)", '(B)'),  font.label = list(size = 12, face = "plain"))

# ggsave(file = 'figures/taxo_fonctio/old_figure3.jpg', plot = plot_var,
#        width = 2.5,  height = 2.5, scale = 3)

################################################################################
################################### temporal evolution

head(df_pca_fonctio)
head(df_pca_taxo)

df_pca_taxo_temporel <- df_pca_taxo %>% 
  dplyr::group_by(Year, Quarter) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   n = n_distinct(index_old)) %>% 
  dplyr::mutate(dim1_up = 1.9*dim1_sd/sqrt(n),
                dim2_up = 1.9*dim2_sd/sqrt(n),
                Quarter2 = paste0("Q", Quarter),
                yyqq = Year + Quarter/4)

ggplot(df_pca_taxo_temporel, aes(x = Year, y = n, color = as.factor(Quarter)))+
  geom_point()+geom_path()
df_pca_taxo_temporel$remove <- ifelse((df_pca_taxo_temporel$Quarter %in% c(1, 3) & 
                                        df_pca_taxo_temporel$Year %between% c(1994, 1996)) |
                                        (df_pca_taxo_temporel$Quarter %in% c(4) & 
                                           df_pca_taxo_temporel$Year %between% c(1994, 2002)),
                                      "yes", "no")
df_pca_taxo_temporel <-df_pca_taxo_temporel %>% 
  dplyr::filter(remove == 'no') %>% 
  dplyr::select(-remove)

# for(rep_q in 1:4){
#   
#   dede <- df_pca_taxo_temporel[df_pca_taxo_temporel$Quarter == rep_q, ]
#   res_nb <- NbClust(dede[,c("dim1_m", "dim2_m")], method = 'kmeans',
#                     min.nc = 2, max.nc = 4)
#   df_pca_taxo_temporel[df_pca_taxo_temporel$Quarter == rep_q, 'cluster'] <- res_nb$Best.partition
#   
# }
# df_pca_taxo_temporel$cluster <- as.factor(df_pca_taxo_temporel$cluster)

  
d1 = ggplot(df_pca_taxo_temporel,
            aes(x = dim1_m, y = dim2_m, color = Year)) + 
  geom_segment(aes(x = dim1_m - dim1_up, xend = dim1_m + dim1_up,
                   y = dim2_m, yend = dim2_m), alpha = 0.2) +
  geom_segment(aes(x = dim1_m , xend = dim1_m,
                   y = dim2_m - dim2_up, yend = dim2_m + dim2_up),
               alpha = 0.2) + xlab("PC1") +  ylab("PC2") +
  scale_color_viridis_c() + geom_path(alpha = 0.5) + 
  scale_shape(guide = 'none') + 
  geom_point(size = 2) + 
  facet_wrap(. ~ Quarter2, scales = 'free', ncol = 4) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        # legend.position = "none",
        axis.text = element_text(size=7))
d1

#######

df_pca_functio_temporel <- data.frame(df_pca_fonctio) %>% 
  dplyr::group_by(Year, Quarter) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   n = n_distinct(index_old)) %>% 
  dplyr::mutate(dim1_up = 1.9*dim1_sd/sqrt(n),
                dim2_up = 1.9*dim2_sd/sqrt(n),
                Quarter2 = paste0("Q", Quarter),
                yyqq = Year + Quarter/4)

df_pca_functio_temporel$remove <- ifelse((df_pca_functio_temporel$Quarter %in% c(1, 3) & 
                                            df_pca_functio_temporel$Year %between% c(1994, 1996)) |
                                        (df_pca_functio_temporel$Quarter %in% c(4) & 
                                           df_pca_functio_temporel$Year %between% c(1994, 2002)),
                                      "yes", "no")
df_pca_functio_temporel <-df_pca_functio_temporel %>% 
  dplyr::filter(remove == 'no') %>% 
  dplyr::select(-remove)

df_pca_functio_temporel$cluster <- NA

# for(rep_q in 1:4){
#   
#   dede <- df_pca_functio_temporel[df_pca_functio_temporel$Quarter == rep_q, ]
#   res_nb <- NbClust(dede[,c("dim1_m", "dim2_m")], method = 'kmeans',
#           min.nc = 2, max.nc = 4)
#   df_pca_functio_temporel[df_pca_functio_temporel$Quarter == rep_q, 'cluster'] <- res_nb$Best.partition
#   
# }
# df_pca_functio_temporel$cluster <- as.factor(df_pca_functio_temporel$cluster)


d1_f = ggplot(df_pca_functio_temporel,
              aes(x = dim1_m, y = dim2_m, color = Year)) + 
  geom_segment(aes(x = dim1_m - dim1_up, xend = dim1_m + dim1_up,
                   y = dim2_m, yend = dim2_m), alpha = 0.2) +
  geom_segment(aes(x = dim1_m , xend = dim1_m,
                   y = dim2_m - dim2_up, yend = dim2_m + dim2_up),
               alpha = 0.2) + xlab("MF1") +  ylab("MF2") +
  scale_color_viridis_c() + geom_path(alpha = 0.5) +  geom_point( size = 2) + 
  facet_wrap(. ~ Quarter2, scales = 'free', ncol = 4) + 
  theme_classic() +
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        # legend.position = "none",
        axis.text = element_text(size=7))



d3_last <-ggarrange(d1,d1_f, nrow = 2, common.legend = TRUE, 
                                labels = c("A", 'B'),  legend = 'right',
                    font.label = list(size = 12, face = "plain")) 

d3_last2 = annotate_figure(d3_last, 
                            right = textGrob("Taxonomic", rot=-90,
                                             # vjust=-1, 
                                             hjust = 2.1,
                                             gp = gpar(cex = 1.3)))

d3_last3 = annotate_figure(d3_last2, 
                            right = textGrob("Functional", rot=-90,
                                             vjust=2.2,
                                             hjust = -1,
                                             gp = gpar(cex = 1.3)))


ggsave(file = 'figures/taxo_fonctio/figure3.jpg', plot = d3_last3,
       width = 3.2,  height = 2, scale = 3)



################################################################################
################### temporal evolution by bioregion  ###########################
################################################################################
df_pca_fonctio$Ecoregion <- as.factor(df_pca_fonctio$Ecoregion)
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_fonctio$Ecoregion)[levels(df_pca_fonctio$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

head(df_pca_fonctio)

df_pca_fonctio_sm <- df_pca_fonctio %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   # dim4_m = mean(Dim.4),
                   # dim4_sd = sd(Dim.4),
                   n = n_distinct(index_old))

df_pca_fonctio_sm$Ecoregion <- as.factor(df_pca_fonctio_sm$Ecoregion)
df_pca_fonctio_sm_melt <- reshape2::melt(df_pca_fonctio_sm,
                                         id.vars = c("Ecoregion","Year", 'n'),
                                         measure.vars = c("dim1_m", "dim2_m",
                                                          "dim3_m"),
                                         variable.name = "dimension",
                                         value.name = "mean_value")
levels(df_pca_fonctio_sm_melt$dimension) <- c("MF1", "MF2", "MF3")
df_pca_fonctio_sm_melt2 <- reshape2::melt(df_pca_fonctio_sm,
                                         id.vars = c("Ecoregion","Year"),
                                         measure.vars = c("dim1_sd", "dim2_sd", 
                                                          "dim3_sd"),
                                         variable.name = "dimension",
                                         value.name = "sd_value")
levels(df_pca_fonctio_sm_melt2$dimension) <- c("MF1", "MF2", "MF3")
df_pca_fonctio_sm_melt <- merge(df_pca_fonctio_sm_melt, df_pca_fonctio_sm_melt2)

levels(df_pca_fonctio_sm_melt$Ecoregion)  
df_pca_fonctio_sm_melt$Ecoregion <- factor(df_pca_fonctio_sm_melt$Ecoregion,
                                           levels = c("Baltic Sea" , "Greater\nNorth Sea", 
                                                      "Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))


df_pca_fonctio_sm_melt$pval <- NA
df_pca_fonctio_sm_melt$coeff <- NA

for(rep_region in unique(df_pca_fonctio_sm_melt$Ecoregion)){
  for(rep_dim in unique(df_pca_fonctio_sm_melt$dimension)){
    # rep_region = "Celtic Seas"
    # rep_dim = "PC1"
    dede <- df_pca_fonctio_sm_melt %>% 
      dplyr::filter(Ecoregion == rep_region & dimension == rep_dim)
    
  # my_lm <- lme(mean_value ~ Year, data = dede,
  #              random=~1|dimension,correlation=corAR1())
  # 
  # df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
  #                          df_pca_fonctio_sm_melt$dimension == rep_dim,
  #                        "pval"] <- summary(my_lm)$tTable[2,5]
  # df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
  #                          df_pca_fonctio_sm_melt$dimension == rep_dim,
  #                        "coeff"] <- summary(my_lm)$tTable[2,1]
    
    my_lm <- lm(mean_value ~ Year, data = dede)
    df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
                             df_pca_fonctio_sm_melt$dimension == rep_dim,
                           "pval"] <- summary(my_lm)$coefficients[2,4]
    df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
                             df_pca_fonctio_sm_melt$dimension == rep_dim,
                           "coeff"] <- summary(my_lm)$coefficients[2,1]
  
    }
}

levels(df_pca_fonctio_sm_melt$Ecoregion)  
df_pca_fonctio_sm_melt$Ecoregion <- factor(df_pca_fonctio_sm_melt$Ecoregion,
                                           levels = c("Baltic Sea" , "Greater\nNorth Sea", 
                                                      "Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))
seuil_pval = 0.09

df_pca_fonctio_sm_melt_significant <- df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$pval < seuil_pval, ]
df_pca_fonctio_sm_melt_significant$Ecoregion <- factor(df_pca_fonctio_sm_melt_significant$Ecoregion,
                                                       exclude = TRUE)
levels(df_pca_fonctio_sm_melt_significant$Ecoregion)  
df_pca_fonctio_sm_melt_significant$Ecoregion <- factor(df_pca_fonctio_sm_melt_significant$Ecoregion,
                                           levels = c("Baltic Sea" , "Greater\nNorth Sea", 
                                                      "Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))

df_pca_fonctio_sm_melt_notsignificant <- df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$pval >= seuil_pval, ]
levels(df_pca_fonctio_sm_melt_notsignificant$Ecoregion)  
df_pca_fonctio_sm_melt_notsignificant$Ecoregion <- factor(df_pca_fonctio_sm_melt_notsignificant$Ecoregion,
                                           levels = c("Baltic Sea" , "Greater\nNorth Sea","Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))


df_text1 <- data.frame(year = rep(2021.5, 6), 
                       mean_value = c(rep(c(-Inf, Inf), 3)),
                       hjust_v = c(rep(c(-0.15, 1.2), 3)),
                       dimension = rep(c("MF1", "MF2","MF3"), each = 2),
                       label = c("Quick\ngrowth", "Hight TL,\nLate repro",
                                 "Benthivorous", "Generalist",
                                 "Non-guarder","Guarder"))
df_segm1 <- data.frame(x = rep(2020, 6), 
                       xend = rep(2020, 6),
                       y =    rep(c(-0.2, 0.2),3), 
                       yend = c(rep(c(-Inf, Inf), 3)),
                       dimension = rep(c("MF1", "MF2","MF3"), each = 2))



a <- ggplot(df_pca_fonctio_sm_melt_significant) +
  scale_fill_manual(drop = FALSE, values = eight_cols) + 
  scale_color_manual(drop = FALSE, values = eight_cols)  + 

  geom_ribbon(data = df_pca_fonctio_sm_melt_significant, 
              aes(x = Year,  ymin = mean_value - 1.9*sd_value/sqrt(n),
                   ymax = mean_value + 1.9*sd_value/sqrt(n),
                   fill = Ecoregion), alpha = 0.2, color = NA) +
  geom_point(data = df_pca_fonctio_sm_melt_significant,
             aes(x = Year, y = mean_value, color = Ecoregion), 
             size = 2) + 
  geom_path(data = df_pca_fonctio_sm_melt_significant,
            aes(x = Year, y = mean_value, color = Ecoregion), 
            alpha = 0.5) + 
  
  geom_ribbon(data = df_pca_fonctio_sm_melt_notsignificant,
              aes(x = Year, 
                  ymin = mean_value - 1.9*sd_value/sqrt(n),
                  ymax = mean_value + 1.9*sd_value/sqrt(n),
                  fill = Ecoregion),alpha = 0.03) +
  geom_point(data = df_pca_fonctio_sm_melt_notsignificant,
             aes(x = Year, y = mean_value, color = Ecoregion), alpha = 0.05, size = 2) + 
  geom_path(data = df_pca_fonctio_sm_melt_notsignificant,
            aes(x = Year, y = mean_value, group = Ecoregion, color = Ecoregion),
            alpha = 0.05) +

  facet_wrap(~ dimension, scales = 'free_y', ncol = 3) + 

  geom_text(data = df_text1,
            aes(x = year, y = mean_value, label = label, hjust = hjust_v),
           angle = 90, lineheight= 0.7) +
  geom_segment(data = df_segm1, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.3,"cm"))) +
  
  
  theme_classic() + xlim(1993, 2022.3) + ylab("Mean values") +

  theme(panel.background = element_rect(fill = NA, color = 'black'),
        legend.title = element_blank(),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey95'))
a


# df_result_trends_taxo <- df_pca_taxo_sm_melt %>% 
#   dplyr::group_by(Ecoregion , dimension  ) %>% 
#   dplyr::summarise(slope = mean(coeff),
#                    pval = mean(pval)) %>% 
#   dplyr::filter(pval <= seuil_pval) 
# 
# df_result_trends_fonctio <- df_pca_fonctio_sm_melt %>% 
#   dplyr::group_by(Ecoregion , dimension  ) %>% 
#   dplyr::summarise(slope = mean(coeff),
#                    pval = mean(pval)) %>% 
#   dplyr::filter(pval <= seuil_pval)
# data.frame(df_result_trends_taxo)
# data.frame(df_result_trends_fonctio)


############################################################################################################ 
head(df_pca_taxo)
df_pca_taxo$Ecoregion <- as.factor(df_pca_taxo$Ecoregion)
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_pca_taxo$Ecoregion)[levels(df_pca_taxo$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_pca_taxo_sm <- df_pca_taxo %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   n = n_distinct(index_old))

df_pca_taxo_sm$Ecoregion <- as.factor(df_pca_taxo_sm$Ecoregion)
df_pca_taxo_sm_melt <- reshape2::melt(df_pca_taxo_sm,
                                         id.vars = c("Ecoregion","Year", 'n'),
                                         measure.vars = c("dim1_m", "dim2_m",
                                                          "dim3_m"),
                                         variable.name = "dimension",
                                         value.name = "mean_value")
levels(df_pca_taxo_sm_melt$dimension) 
levels(df_pca_taxo_sm_melt$dimension) <- c("PC1", "PC2", "PC3")
df_pca_taxo_sm_melt2 <- reshape2::melt(df_pca_taxo_sm,
                                          id.vars = c("Ecoregion","Year"),
                                          measure.vars = c("dim1_sd", "dim2_sd",
                                                           "dim3_sd"),
                                          variable.name = "dimension",
                                          value.name = "sd_value")
levels(df_pca_taxo_sm_melt2$dimension) <-c("PC1", "PC2", "PC3")
df_pca_taxo_sm_melt <- merge(df_pca_taxo_sm_melt, df_pca_taxo_sm_melt2)

levels(df_pca_taxo_sm_melt$Ecoregion)  
df_pca_taxo_sm_melt$Ecoregion <- factor(df_pca_taxo_sm_melt$Ecoregion,
                                           levels = c("Baltic Sea" , "Greater\nNorth Sea", 
                                                      "Celtic Seas", 
                                                      "Bay of Biscay and\nIberian Coast",
                                                      "Western\nMed. Sea",
                                                      "Ionian Sea and\nCentral Med. Sea",
                                                      "Adriatic Sea"    ,
                                                      "Aegean-Levantine\nSea"))
seuil_pval = 0.09

for(rep_region in unique(df_pca_taxo_sm_melt$Ecoregion)){
  for(rep_dim in unique(df_pca_taxo_sm_melt$dimension)){
    dede <- df_pca_taxo_sm_melt %>% 
      dplyr::filter(Ecoregion == rep_region & dimension == rep_dim)
    
    # my_lm <- lme(mean_value ~ Year, data = dede,
    #              random=~1|dimension,correlation=corAR1())
    # 
    # df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
    #                          df_pca_fonctio_sm_melt$dimension == rep_dim,
    #                        "pval"] <- summary(my_lm)$tTable[2,5]
    # df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
    #                          df_pca_fonctio_sm_melt$dimension == rep_dim,
    #                        "coeff"] <- summary(my_lm)$tTable[2,1]
    
    my_lm <- lm(mean_value ~ Year, data = dede)
    df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$Ecoregion == rep_region &
                          df_pca_taxo_sm_melt$dimension == rep_dim,
                           "pval"] <- summary(my_lm)$coefficients[2,4]
    df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$Ecoregion == rep_region &
                          df_pca_taxo_sm_melt$dimension == rep_dim,
                           "coeff"] <- summary(my_lm)$coefficients[2,1]
    
  }
}


df_text2 <- data.frame(year = rep(2021.5, 6), 
                       mean_value = c(rep(c(-Inf, Inf), 3)),
                       hjust_v = c(rep(c(-0.15, 1.2), 3)),
                       dimension = rep(c("PC1", "PC2","PC3"), each = 2),
                       label = c("Limanda\nlimanda", "Merluccius\nmerluccius" ,
                                 "Sprattus\nsprattus", "Melanogrammus\naeglefinus" ,
                                 "Pleuronectes\nplatessa" ,"Clupea\nharengus" ))
df_segm2 <- data.frame(x = rep(2021.75, 6), 
                       xend = rep(2021.75, 6),
                       y =    rep(c(-0.04, 0.04),3), 
                       yend = c(rep(c(-Inf, Inf), 3)),
                       dimension = rep(c("PC1", "PC2","PC3"), each = 2))


c <- ggplot(df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$pval <= seuil_pval, ]) +
  scale_fill_manual(drop = FALSE, values = eight_cols) + 
  scale_color_manual(drop = FALSE, values = eight_cols)  + 
  
  geom_ribbon(aes(x = Year,
                  ymin = mean_value - 1.9*sd_value/sqrt(n),
                  ymax = mean_value + 1.9*sd_value/sqrt(n),
                  fill = Ecoregion), alpha = 0.2, color = NA) +
  geom_point(aes(x = Year, y = mean_value, color = Ecoregion),
             size = 2) +
  geom_path(aes(x = Year, y = mean_value, color = Ecoregion),
            alpha = 0.5) +
  
  geom_ribbon(data = df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$pval >= seuil_pval, ],
              aes(x = Year,
                  ymin = mean_value - 1.9*sd_value/sqrt(n),
                  ymax = mean_value + 1.9*sd_value/sqrt(n),
                  fill = Ecoregion), alpha = 0.03, color = NA) +
  geom_point(data = df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$pval >= seuil_pval, ],
             aes(x = Year, y = mean_value, color = Ecoregion),
             size = 2, alpha = 0.05) +
  geom_path(data = df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$pval >= seuil_pval, ],
            aes(x = Year, y = mean_value, color = Ecoregion),
            alpha = 0.05) +
  

  facet_wrap(~ dimension, scales = 'free_y', ncol = 4) + 
  
  geom_text(data = df_text2,
            aes(x = year, y = mean_value, label = label, hjust = hjust_v),
            angle = 90, lineheight= 1, fontface = "italic") +
  geom_segment(data = df_segm2, aes(x = x, xend = xend, y = y, yend = yend),
               arrow = arrow(length = unit(0.3,"cm"))) +
  
  theme_classic() + xlim(1993, 2023)  +ylab("Mean values") +
  
  theme(panel.background = element_rect(fill = NA, color = 'black'),
        legend.title = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_line(color = 'grey85'),
        panel.grid.minor = element_line(color = 'grey95'))
c

fig4_last <- ggarrange(c, a, nrow = 2, legend = "bottom",
                       common.legend = TRUE, heights = c(1, 1.05))

fig4_last2 <- annotate_figure(fig4_last, 
                            right = textGrob("Taxonomic", rot=-90,
                                             # vjust=-1, 
                                             hjust = 2.4,
                                             gp = gpar(cex = 1.3)))

fig4_last3 <- annotate_figure(fig4_last2, 
                            right = textGrob("Functional", rot=-90,
                                             vjust=2.2,
                                             hjust = -0.6,
                                             gp = gpar(cex = 1.3)))

ggsave(file = 'figures/taxo_fonctio/figure4.jpg', plot = fig4_last3,
       width = 3.5,  height = 2.4, scale = 3)


df_pca_taxo_sm_melt

df_pca_fonctio_sm_melt %>% 
  dplyr::filter(dimension != "MF3") %>% 
  dplyr::group_by( Ecoregion , dimension) %>% 
  dplyr::summarise(pval = mean(pval),
                   coeff = mean(coeff))

df_pca_taxo_sm_melt %>% 
  dplyr::filter(dimension != "PC3") %>% 
  dplyr::group_by( Ecoregion , dimension) %>% 
  dplyr::summarise(pval = mean(pval),
                   coeff = mean(coeff))

################################################################################
################################### GLS effet env ######################
################################################################################

head(df_pca_fonctio)
df_pca_fonctio$log_chloro <- log(df_pca_fonctio$chloro_mea + 1)
df_pca_fonctio$log_mlotst <- log(df_pca_fonctio$mlotst_mea + 1)

df_pca_fonctio$log_curr_bottom <- log(df_pca_fonctio$curr_bottom_mea) + 10
df_pca_fonctio$log_curr_surf <- log(df_pca_fonctio$curr_surf_mea) + 10


df_pca_fonctio_sm <- df_pca_fonctio %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   chloro = mean(log_chloro, na.rm = TRUE),
                   mlotst = mean(log_mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom_mea, na.rm = TRUE),
                   temp_surf = mean(temp_surf_mea, na.rm = TRUE),
                   curr_bottom = mean(log_curr_bottom, na.rm = TRUE),
                   curr_surf = mean(log_curr_surf, na.rm = TRUE),
                   n = n_distinct(index_old))
df_pca_fonctio_sm$chloro_scale <- scale(df_pca_fonctio_sm$chloro)
df_pca_fonctio_sm$mlotst_scale <- scale(df_pca_fonctio_sm$mlotst)
df_pca_fonctio_sm$oxy_bottom_scale <- scale(df_pca_fonctio_sm$oxy_bottom)
df_pca_fonctio_sm$oxy_surf_scale <- scale(df_pca_fonctio_sm$oxy_surf)
df_pca_fonctio_sm$temp_bottom_scale <- scale(df_pca_fonctio_sm$temp_bottom)
df_pca_fonctio_sm$temp_surf_scale <- scale(df_pca_fonctio_sm$temp_surf)
df_pca_fonctio_sm$curr_bottom_scale <- scale(df_pca_fonctio_sm$curr_bottom)
df_pca_fonctio_sm$curr_surf_scale <- scale(df_pca_fonctio_sm$curr_surf)

df_pca_fonctio_sm <- data.frame(df_pca_fonctio_sm)
liste_ecoregion <- unique(df_pca_fonctio_sm$Ecoregion )

for(rep_ecoR in liste_ecoregion){
  print(rep_ecoR)
        # rep_ecoR <- 'Celtic Seas'

  for(rep_dim in c('dim1_m','dim2_m','dim3_m')){
    # rep_dim = "dim2_m"
    print(rep_dim)
    
dede <- df_pca_fonctio_sm %>% 
  dplyr::filter(Ecoregion == rep_ecoR) %>% 
  dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                  !is.na(temp_bottom) &  !is.na(temp_surf) &  !is.na(curr_bottom))
ggplot(dede, aes(x = dim1_m, y = temp_surf)) + geom_point() + geom_smooth(method = 'lm')
# df_env <- dede[,c(9:16)]


gls1 = nlme::gls(formula(paste0(rep_dim, " ~  temp_bottom_scale + temp_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + curr_bottom_scale ")),
                 data = dede,
                 correlation = corGaus(form =~ Year))


summary(gls1)

library(MuMIn)
options(na.action = "na.fail")
dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
options(na.action = "na.omit")
aa = subset(dredge.glm1.BIC.p, delta < 2)

mod <- get.models(aa, 1)[[1]]
summary(mod)

table_coeff <- as.data.frame(summary(mod)$tTable)
table_coeff$variable_env <- rownames(table_coeff)
table_coeff <- reshape2::melt(table_coeff)
table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )

df_results <- table_coeff
df_results$Ecoregion <- rep_ecoR
df_results$dim <- rep_dim


size_tab <-dim(dredge.glm1.BIC.p)[1]
test1 <- as.data.frame(dredge.glm1.BIC.p[1:5, ])
test1$Ecoregion <- rep_ecoR
test1$Dimension <- rep_dim

if(rep_ecoR == liste_ecoregion[1] & rep_dim == 'dim1_m'){
  df_results_last <- df_results
  
  rm(test2)
  test2 <- test1

} else {
  df_results_last <- rbind(df_results_last, df_results)
    
  
  test2 <- rbind(test2, test1)
  rm(test1)
  
  }



}
  
}
names(df_results_last)[names(df_results_last) == "p-value"] <- 'p_val'
names(df_results_last)[names(df_results_last) == "Std.Error"] <- 'std_error'

df_results_last2 <- df_results_last %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(Ecoregion, dim)

df_aic_fonctio_regions <- test2



# write.csv2(file = "data/processed/AIC_GLS_spatial_all_models.csv", dredge.glm1.BIC.p,
#            row.names = FALSE)
# 
# top10 = subset(dredge.glm1.BIC.p, 1:11)
# 
# for(i in 1:11){
#   mod <- get.models(top10, i)[[1]]
# 
# sink(paste0("data/processed/modele_top", i, ".txt"))
# print(summary(mod))
# sink() 
# 
# 
# }



############ spatial 

df_pca_fonctio <-  data.frame(result_mfa$ind$coord)
df_pca_fonctio$index_old <- df_CMW_all2$index_old

df_pca_fonctio <- merge(df_complementary, df_pca_fonctio)
df_pca_fonctio$log_chloro <- log(df_pca_fonctio$chloro_mea + 1)
df_pca_fonctio$log_mlotst <- log(df_pca_fonctio$mlotst_mea + 1)

df_pca_fonctio$log_curr_bottom <- log(df_pca_fonctio$curr_bottom_mea) + 10
df_pca_fonctio$log_curr_surf <- log(df_pca_fonctio$curr_surf_mea) + 10


df_pca_fonctio2 <- df_pca_fonctio %>% 
  dplyr::group_by(index_old,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   pco3_m  =  mean(Dim.3),
                   pco4_m  =  mean(Dim.4),
                   pco5_m  =  mean(Dim.5),
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro = mean(log_chloro, na.rm = TRUE),
                   mlotst = mean(log_mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom_mea, na.rm = TRUE),
                   temp_surf = mean(temp_surf_mea, na.rm = TRUE),
                   curr_bottom = mean(log_curr_bottom, na.rm = TRUE),
                   curr_surf = mean(log_curr_surf, na.rm = TRUE),
                   n = n_distinct(index_old))


df_pca_fonctio3 <- df_pca_fonctio2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id) %>% 
  dplyr::summarise(dim1_m = mean(pca1_m),
                   dim2_m  =  mean(pco2_m),
                   dim3_m  =  mean(pco3_m),
                   dim4_m  =  mean(pco4_m),
                   dim5_m  =  mean(pco5_m),
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro = mean(chloro, na.rm = TRUE),
                   mlotst = mean(mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom, na.rm = TRUE),
                   temp_surf = mean(temp_surf, na.rm = TRUE),
                   curr_bottom = mean(curr_bottom, na.rm = TRUE),
                   curr_surf = mean(curr_surf, na.rm = TRUE))

df_pca_fonctio3$chloro_scale <- scale(df_pca_fonctio3$chloro)
df_pca_fonctio3$mlotst_scale <- scale(df_pca_fonctio3$mlotst)
df_pca_fonctio3$oxy_bottom_scale <- scale(df_pca_fonctio3$oxy_bottom)
df_pca_fonctio3$oxy_surf_scale <- scale(df_pca_fonctio3$oxy_surf)
df_pca_fonctio3$temp_bottom_scale <- scale(df_pca_fonctio3$temp_bottom)
df_pca_fonctio3$temp_surf_scale <- scale(df_pca_fonctio3$temp_surf)
df_pca_fonctio3$curr_bottom_scale <- scale(df_pca_fonctio3$curr_bottom)
df_pca_fonctio3$curr_surf_scale <- scale(df_pca_fonctio3$curr_surf)
df_pca_fonctio3$depth_scale <- scale(df_pca_fonctio3$depth)
df_pca_fonctio3$depth_span_scale <- scale(df_pca_fonctio3$depth_span)

df_pca_fonctio3 <- data.frame(df_pca_fonctio3)
df_pca_fonctio3$index <-  row_number(df_pca_fonctio3$x_my_spatial_id)


############ extract fishing 

for(rep_yy in 2012:2020){  print(rep_yy)
  df_yy_qq <- read.csv2(file = paste0('data/raw/others/GFW/df_yy_qq_', rep_yy,'_all.csv'))
  if(rep_yy == 2012){  df_yy_qq_last <- df_yy_qq} else {
    df_yy_qq_last <- rbind(df_yy_qq_last,df_yy_qq )
    rm(df_yy_qq) }
}

df_GFW <- df_yy_qq_last %>% 
  dplyr::group_by(lon2, lat2) %>%
  dplyr::summarise(mean_fishing = mean(sum_fishing_hours2)) %>% 
  dplyr::mutate(log_fishing = log(mean_fishing + 1))
df_GFW_matrix <- reshape2::acast(df_GFW, lon2 ~ lat2, value.var = "log_fishing")


ras_fishing <- raster(rotate(rotate(rotate(df_GFW_matrix))),
                      xmn = min(df_yy_qq_last$lon2),xmx = max(df_yy_qq_last$lon2),
                      ymn = min(df_yy_qq_last$lat2), ymx = max(df_yy_qq_last$lat2),
                      CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# plot(ras_fishing)

df_for_extract_f <- data.frame(df_pca_fonctio3) %>%  
  dplyr::select(index, x_my_spatial_id, y_my_spatial_id)
names(df_for_extract_f)[c(2,3)] <- c("x", "y")
coordinates(df_for_extract_f) <- ~ x + y
proj4string(df_for_extract_f) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

df_for_extract_f$fishing_hours <- raster::extract(ras_fishing, df_for_extract_f)
df_fishing_f <- as.data.frame(df_for_extract_f)

df_fishing_f <- df_fishing_f %>% dplyr::select(index, fishing_hours)
df_fishing_f$fishing_hours_scale <- scale(df_fishing_f$fishing_hours)

df_pca_fonctio3 <- merge(df_pca_fonctio3 ,df_fishing_f )

############ 

for(rep_dim in c('dim1_m','dim2_m','dim3_m')){
  print(rep_dim)
  
  dede <- df_pca_fonctio3 %>% 
    dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                    !is.na(temp_bottom) &  !is.na(temp_surf) &  
                    !is.na(depth) &  !is.na(curr_bottom) &  !is.na(fishing_hours_scale))
  
  # df_env <- dede[,c(9:16)]
  
  
  gls1 = nlme::gls(formula(paste0(rep_dim, " ~ depth_scale + depth_span_scale + 
                          temp_bottom_scale + temp_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + curr_bottom_scale +
                                  fishing_hours_scale")),
                   data = dede,
                   correlation = corGaus(form =~ x_my_spatial_id + y_my_spatial_id))

  summary(gls1)
  
  library(MuMIn)
  options(na.action = "na.fail")
  dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
  options(na.action = "na.omit")
  aa = subset(dredge.glm1.BIC.p, delta < 2)
  
  mod <- get.models(aa, 1)[[1]]
  summary(mod)
  
   table_coeff <- as.data.frame(summary(mod)$tTable)
  table_coeff$variable_env <- rownames(table_coeff)
  table_coeff <- reshape2::melt(table_coeff)
  table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
  table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )
  
  df_results <- table_coeff
  df_results$dim <- rep_dim
  
  
  size_tab <-dim(dredge.glm1.BIC.p)[1]
  test1 <- as.data.frame(dredge.glm1.BIC.p[1:5, ])
  # test1$Ecoregion <- rep_ecoR
  test1$Dimension <- rep_dim
  
  
  if( rep_dim == 'dim1_m'){
    df_results_last_spatial <- df_results
    
    rm(test2)
    test2 <- test1
    
  } else {
    df_results_last_spatial <- rbind(df_results_last_spatial, df_results)
    
    test2 <- rbind(test2, test1)
    rm(test1)
  }
  
}

names(df_results_last_spatial)[names(df_results_last_spatial) == "p-value"] <- 'p_val'
names(df_results_last_spatial)[names(df_results_last_spatial) == "Std.Error"] <- 'std_error'

df_results_last_spatial2 <- df_results_last_spatial %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(dim)

df_results_last_spatial2$Ecoregion <- 'Spatial distribution'

df_all_func <- rbind(df_results_last2, df_results_last_spatial2 )


test2$Ecoregion <- 'Spatial distribution'

df_aic_fonctio_spatial <- test2
df_aic_fonctio_regions


################################### GLS effet env taxo ######################
################################################################################

head(df_pca_taxo)
df_pca_taxo$log_chloro <- log(df_pca_taxo$chloro_mea + 1)
df_pca_taxo$log_mlotst <- log(df_pca_taxo$mlotst_mea + 1)

df_pca_taxo$log_curr_bottom <- log(df_pca_taxo$curr_bottom_mea) + 10
df_pca_taxo$log_curr_surf <- log(df_pca_taxo$curr_surf_mea) + 10


df_pca_taxo_sm <- df_pca_taxo %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   chloro = mean(log_chloro, na.rm = TRUE),
                   mlotst = mean(log_mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom_mea, na.rm = TRUE),
                   temp_surf = mean(temp_surf_mea, na.rm = TRUE),
                   curr_bottom = mean(log_curr_bottom, na.rm = TRUE),
                   curr_surf = mean(log_curr_surf, na.rm = TRUE),
                   n = n_distinct(index_old))
df_pca_taxo_sm$chloro_scale <- scale(df_pca_taxo_sm$chloro)
df_pca_taxo_sm$mlotst_scale <- scale(df_pca_taxo_sm$mlotst)
df_pca_taxo_sm$oxy_bottom_scale <- scale(df_pca_taxo_sm$oxy_bottom)
df_pca_taxo_sm$oxy_surf_scale <- scale(df_pca_taxo_sm$oxy_surf)
df_pca_taxo_sm$temp_bottom_scale <- scale(df_pca_taxo_sm$temp_bottom)
df_pca_taxo_sm$temp_surf_scale <- scale(df_pca_taxo_sm$temp_surf)
df_pca_taxo_sm$curr_bottom_scale <- scale(df_pca_taxo_sm$curr_bottom)
df_pca_taxo_sm$curr_surf_scale <- scale(df_pca_taxo_sm$curr_surf)

df_pca_taxo_sm <- data.frame(df_pca_taxo_sm)
liste_ecoregion <- unique(df_pca_taxo_sm$Ecoregion )

for(rep_ecoR in liste_ecoregion){
  print(rep_ecoR)
  # rep_ecoR <- 'Greater North Sea'
  
  for(rep_dim in c('dim1_m','dim2_m','dim3_m')){
    # rep_dim = "dim1_m"
    print(rep_dim)
    
    dede <- df_pca_taxo_sm %>% 
      dplyr::filter(Ecoregion == rep_ecoR) %>% 
      dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                      !is.na(temp_bottom) &  !is.na(temp_surf) &  !is.na(curr_bottom))
    ggplot(dede, aes(x = dim1_m, y = temp_surf)) + geom_point() + geom_smooth(method = 'lm')
    # df_env <- dede[,c(9:16)]
    
    
    gls1 = nlme::gls(formula(paste0(rep_dim, " ~  temp_bottom_scale + temp_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + curr_bottom_scale ")),
                     data = dede,
                     correlation = corGaus(form =~ Year))

    summary(gls1)
    
    library(MuMIn)
    options(na.action = "na.fail")
    dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
    options(na.action = "na.omit")
    aa = subset(dredge.glm1.BIC.p, delta < 2)
    
    mod <- get.models(aa, 1)[[1]]
    summary(mod)
    
    table_coeff <- as.data.frame(summary(mod)$tTable)
    table_coeff$variable_env <- rownames(table_coeff)
    table_coeff <- reshape2::melt(table_coeff)
    table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
    table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )
    
    df_results <- table_coeff
    df_results$Ecoregion <- rep_ecoR
    df_results$dim <- rep_dim
    
    size_tab <-dim(dredge.glm1.BIC.p)[1]
    test1 <- as.data.frame(dredge.glm1.BIC.p[1:5, ])
    test1$Ecoregion <- rep_ecoR
    test1$Dimension <- rep_dim
    
    
    if(rep_ecoR == liste_ecoregion[1] & rep_dim == 'dim1_m'){
      df_results_last_t <- df_results
      rm(test2)
      test2 <- test1
      
    } else {
      df_results_last_t <- rbind(df_results_last_t, df_results)
      test2 <- rbind(test2, test1)
      rm(test1)
    }
    
  }
  
}
names(df_results_last_t)[names(df_results_last_t) == "p-value"] <- 'p_val'
names(df_results_last_t)[names(df_results_last_t) == "Std.Error"] <- 'std_error'

df_results_last_t2 <- df_results_last_t %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(Ecoregion, dim)

df_aic_taxo_regions <- test2
df_aic_fonctio_regions
df_aic_fonctio_spatial


############ spatial ############ 

df_pca_taxo2 <- df_pca_taxo %>% 
  dplyr::group_by(index_old,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   pco3_m  =  mean(Dim.3),
                   pco4_m  =  mean(Dim.4),
                   pco5_m  =  mean(Dim.5),
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro = mean(log_chloro, na.rm = TRUE),
                   mlotst = mean(log_mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom_mea, na.rm = TRUE),
                   temp_surf = mean(temp_surf_mea, na.rm = TRUE),
                   curr_bottom = mean(log_curr_bottom, na.rm = TRUE),
                   curr_surf = mean(log_curr_surf, na.rm = TRUE),
                   n = n_distinct(index_old))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id) %>% 
  dplyr::summarise(dim1_m = mean(pca1_m),
                   dim2_m  =  mean(pco2_m),
                   dim3_m  =  mean(pco3_m),
                   dim4_m  =  mean(pco4_m),
                   dim5_m  =  mean(pco5_m),
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro = mean(chloro, na.rm = TRUE),
                   mlotst = mean(mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom, na.rm = TRUE),
                   temp_surf = mean(temp_surf, na.rm = TRUE),
                   curr_bottom = mean(curr_bottom, na.rm = TRUE),
                   curr_surf = mean(curr_surf, na.rm = TRUE))

df_pca_taxo3$chloro_scale <- scale(df_pca_taxo3$chloro)
df_pca_taxo3$mlotst_scale <- scale(df_pca_taxo3$mlotst)
df_pca_taxo3$oxy_bottom_scale <- scale(df_pca_taxo3$oxy_bottom)
df_pca_taxo3$oxy_surf_scale <- scale(df_pca_taxo3$oxy_surf)
df_pca_taxo3$temp_bottom_scale <- scale(df_pca_taxo3$temp_bottom)
df_pca_taxo3$temp_surf_scale <- scale(df_pca_taxo3$temp_surf)
df_pca_taxo3$curr_bottom_scale <- scale(df_pca_taxo3$curr_bottom)
df_pca_taxo3$curr_surf_scale <- scale(df_pca_taxo3$curr_surf)
df_pca_taxo3$depth_scale <- scale(df_pca_taxo3$depth)
df_pca_taxo3$depth_span_scale <- scale(df_pca_taxo3$depth_span)

df_pca_taxo3 <- data.frame(df_pca_taxo3)
df_pca_taxo3$index <-  row_number( df_pca_taxo3$x_my_spatial_id)

############ extract fishing 

df_for_extract <- data.frame(df_pca_taxo3) %>%  
  dplyr::select(index, x_my_spatial_id, y_my_spatial_id)
names(df_for_extract)[c(2,3)] <- c("x", "y")
coordinates(df_for_extract) <- ~ x + y
proj4string(df_for_extract) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

df_for_extract$fishing_hours <- raster::extract(ras_fishing, df_for_extract)
df_fishing <- as.data.frame(df_for_extract)

df_fishing <- df_fishing %>% dplyr::select(index, fishing_hours)
df_fishing$fishing_hours_scale <- scale(df_fishing$fishing_hours)

df_pca_taxo3 <- merge(df_pca_taxo3 ,df_fishing )

############ 

for(rep_dim in c('dim1_m','dim2_m','dim3_m')){
  print(rep_dim)
  
  dede <- df_pca_taxo3 %>% 
    dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                    !is.na(temp_bottom) &  !is.na(temp_surf) &  
                    !is.na(depth) &  !is.na(curr_bottom) &  !is.na(fishing_hours_scale))
  
  # df_env <- dede[,c(9:16)]
  
  
  gls1 = nlme::gls(formula(paste0(rep_dim, " ~ depth_scale + depth_span_scale + temp_bottom_scale +
                                           temp_surf_scale + mlotst_scale + oxy_bottom_scale  +
                                           chloro_scale + curr_bottom_scale +
                                  fishing_hours_scale")),
                   data = dede,
                   correlation = corGaus(form =~ x_my_spatial_id + y_my_spatial_id))
 
  summary(gls1)
  
  library(MuMIn)
  options(na.action = "na.fail")
  dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
  options(na.action = "na.omit")
  aa = subset(dredge.glm1.BIC.p, delta < 2)
  
  mod <- get.models(aa, 1)[[1]]
  summary(mod)
  
  table_coeff <- as.data.frame(summary(mod)$tTable)
  table_coeff$variable_env <- rownames(table_coeff)
  table_coeff <- reshape2::melt(table_coeff)
  table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
  table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )
  
  df_results <- table_coeff
  df_results$dim <- rep_dim
  
  size_tab <-dim(dredge.glm1.BIC.p)[1]
  test1 <- as.data.frame(dredge.glm1.BIC.p[1:5, ])
  # test1$Ecoregion <- rep_ecoR
  test1$Dimension <- rep_dim
  
  
  if( rep_dim == 'dim1_m'){
    df_results_last_spatial_t <- df_results
    rm(test2)
    test2 <- test1
    
  } else {
    df_results_last_spatial_t <- rbind(df_results_last_spatial_t, df_results)
    test2 <- rbind(test2, test1)
    rm(test1)
  }
  
}
names(df_results_last_spatial_t)[names(df_results_last_spatial_t) == "p-value"] <- 'p_val'
names(df_results_last_spatial_t)[names(df_results_last_spatial_t) == "Std.Error"] <- 'std_error'

df_results_last_spatial_t2 <- df_results_last_spatial_t %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(dim)

df_results_last_spatial_t2$Ecoregion <- 'Spatial distribution'

df_all_taxo <- rbind(df_results_last_t2, df_results_last_spatial_t2 )

df_aic_taxo_spatial <- test2
df_aic_taxo_spatial$Ecoregion <- 'Spatial distribution'


df_aic_taxo_regions
df_aic_fonctio_regions

####
df_aic_fonctio_spatial$facet <- 'Functional'
df_aic_taxo_spatial$facet <- 'Taxonomic'

df_aic_spatial <- rbind(df_aic_fonctio_spatial , df_aic_taxo_spatial)
df_aic_spatial <- df_aic_spatial %>% 
  dplyr::filter(Dimension %in% c("dim1_m", "dim2_m")) %>% 
  dplyr::select(facet, Dimension,
                1:10, 
                df   ,  logLik    ,  AICc  )
df_aic_spatial$Dimension <- as.factor(df_aic_spatial$Dimension)
# levels(df_aic_spatial$Dimension) <- c("PC1/MF1", "PC2/MF2")
# write.csv2(file = "data/processed/AIC_GLS_spatial.csv", df_aic_spatial,
#            row.names = FALSE)


####

df_aic_fonctio_regions$facet <- 'Functional'
df_aic_taxo_regions$facet <- 'Taxonomic'

df_aic_regions <- rbind(df_aic_fonctio_regions , df_aic_taxo_regions)
df_aic_regions <- df_aic_regions %>% 
  dplyr::filter(Dimension %in% c("dim1_m", "dim2_m")) %>% 
  dplyr::select(facet, Dimension,Ecoregion,
                1:7, 
                df   ,  logLik    ,  AICc  )
df_aic_regions$Dimension <- as.factor(df_aic_regions$Dimension)
levels(df_aic_regions$Dimension) <- c("PC1/MF1", "PC2/MF2")

df_aic_regions$Ecoregion <- factor(df_aic_regions$Ecoregion, exclude = T)
levels(df_aic_regions$Ecoregion)
levels(df_aic_regions$Ecoregion)[levels(df_aic_regions$Ecoregion) == "Aegean-Levantine\nSea"] <- "Aegean-Levantine Sea"
levels(df_aic_regions$Ecoregion)[levels(df_aic_regions$Ecoregion) == "Ionian Sea and\nCentral Med. Sea"] <- "Ionian Sea and the Central Mediterranean Sea"

levels(df_aic_regions$Ecoregion)[levels(df_aic_regions$Ecoregion) == "Greater\nNorth Sea"] <- "Greater North Sea"
levels(df_aic_regions$Ecoregion)[levels(df_aic_regions$Ecoregion) == "Bay of Biscay and\nIberian Coast"] <- "Bay of Biscay and the Iberian Coast"
levels(df_aic_regions$Ecoregion)[levels(df_aic_regions$Ecoregion) == "Western\nMed. Sea"] <- "Western Mediterranean Sea"

# write.csv2(file = "data/processed/AIC_GLS_regions.csv", df_aic_regions,
#            row.names = FALSE)



#######################  taxo et fonctio together

df_all_taxo$facet <- 'Taxonomic'
df_all_func$facet  <- 'Functional'


df_all_func$Ecoregion <- as.factor(df_all_func$Ecoregion)
levels(df_all_func$Ecoregion)
unique(df_all_func$Ecoregion)


df_all_taxo$Ecoregion <- as.factor(df_all_taxo$Ecoregion)
levels(df_all_taxo$Ecoregion)
levels(df_all_taxo$Ecoregion)[levels(df_all_taxo$Ecoregion) == "Aegean-Levantine\nSea"] <- "Aegean-Levantine Sea"
levels(df_all_taxo$Ecoregion)[levels(df_all_taxo$Ecoregion) == "Ionian Sea and\nCentral Med. Sea"] <- "Ionian Sea and the Central Mediterranean Sea"


df_all <- rbind(df_all_func, df_all_taxo)

levels(df_all$Ecoregion)
# df_all$Ecoregion <- factor(df_all$Ecoregion,
#                            levels = c("Spatial distribution" , "Baltic Sea"     , "Celtic Seas"   ,
#                                       "Bay of Biscay and\nIberian Coast"  ,
#                                       "Western\nMed. Sea" , "Ionian Sea and\nCentral Med. Sea",
#                                       "Adriatic Sea" ,    "Aegean-Levantine\nSea"))

levels(df_all$Ecoregion)[levels(df_all$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_all$Ecoregion)[levels(df_all$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <- "Bay of Biscay and\nIberian Coast" 
levels(df_all$Ecoregion)[levels(df_all$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_all$Ecoregion)[levels(df_all$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_all$Ecoregion)[levels(df_all$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_all$Ecoregion  <- factor(df_all$Ecoregion , exclude = TRUE)
levels(df_all$Ecoregion)
df_all$Ecoregion <- factor(df_all$Ecoregion, 
                           levels = c("Spatial distribution" , "Baltic Sea"     , "Celtic Seas"   ,
                                      "Bay of Biscay and\nIberian Coast",
                                      "Western\nMed. Sea" ,
                                      "Ionian Sea and\nCentral Med. Sea",
                                      "Adriatic Sea" ,    "Aegean-Levantine\nSea"   ))

df_all$variable_env <- as.factor(df_all$variable_env)

levels(df_all$variable_env)
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "chloro_scale"] <- "Surface chlorophyll" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "curr_bottom_scale"] <- "Bottom current" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "mlotst_scale"] <- "Mixed layer depth" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "oxy_bottom_scale"] <- "Bottom oxygen" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "temp_bottom_scale"] <- "Bottom temperature" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "temp_surf_scale"] <- "Surface temperature" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "fishing_hours_scale"] <- "Fishing effort" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "depth_span_scale"] <- "Depth span" 
levels(df_all$variable_env)[levels(df_all$variable_env) ==  "depth_scale"] <- "Depth" 

df_all$dim <- as.factor(df_all$dim)
levels(df_all$dim)
levels(df_all$dim) <- c("PC1/MF1", "PC2/MF2", "PC3/MF3")




df_pval_temporel_taxo <- data.frame(df_pca_taxo_sm_melt) %>% 
  dplyr::group_by(Ecoregion, dimension) %>% 
  dplyr::summarise(pval_tred = mean(pval)) %>% 
  dplyr::mutate(facet = "Taxonomic", 
                dim = dimension) %>% 
  dplyr::select(-dimension)
df_pval_temporel_taxo$dim <- as.factor(df_pval_temporel_taxo$dim)
levels(df_pval_temporel_taxo$dim)
levels(df_pval_temporel_taxo$dim) <- c("PC1/MF1", "PC2/MF2", "PC3/MF3")

df_pval_temporel_fonctio <- data.frame(df_pca_fonctio_sm_melt) %>% 
  dplyr::group_by(Ecoregion, dimension) %>% 
  dplyr::summarise(pval_tred = mean(pval)) %>% 
  dplyr::mutate(facet = "Functional" , 
                dim = dimension) %>% 
  dplyr::select(-dimension)
df_pval_temporel_fonctio$dim <- as.factor(df_pval_temporel_fonctio$dim)
levels(df_pval_temporel_fonctio$dim)
levels(df_pval_temporel_fonctio$dim) <- c("PC1/MF1", "PC2/MF2", "PC3/MF3")



df_pval_temporel_all <- rbind(df_pval_temporel_taxo, df_pval_temporel_fonctio)
df_pval_temporel_all <- data.frame(df_pval_temporel_all)

df_pval_spatial <- data.frame(expand.grid(dim = unique(df_pval_temporel_all$dim),
                                           facet = unique(df_pval_temporel_all$facet),
                                           Ecoregion = 'Spatial distribution'))
df_pval_spatial$pval_tred <- 0.00002
df_pval_temporel_all2 <- rbind(df_pval_temporel_all, df_pval_spatial)

levels(df_all$Ecoregion)
levels(df_pval_temporel_all2$Ecoregion)
df_all2 <- merge(df_all, df_pval_temporel_all2)

df_all2 <- df_all2 %>% 
  dplyr::filter(pval_tred <= seuil_pval )


row1 <- c("Baltic Sea", "Celtic Seas", "Bay of Biscay and\nIberian Coast" ,
          "Greater North Sea"   , "Spatial distribution") 

df_all2$facet <- as.factor(df_all2$facet)
levels(df_all2$facet)
df_all2$facet <- factor(df_all2$facet, levels = c( "Taxonomic" ,"Functional"  ))

aaa1 = ggplot(df_all2[df_all2$Ecoregion %in% row1, ], 
             aes(x = Value, y = variable_env, color = dim)) + 
  geom_pointrange(aes(xmin =Value -  std_error, xmax = Value + std_error,
                      fill = dim), shape=21, fatten = 1.1, size = 1,
                  position = position_dodge2(width = 0.35)) + 
  scale_color_manual(values = three_cols) +
  scale_fill_manual(values = three_cols) +
  theme_classic() + geom_vline(xintercept = 0 ) +xlim(-1, 1) +
  facet_grid(facet ~ Ecoregion, scales = 'free', space='free_y') + 
  xlab("Slope") +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = 'grey90'),
        legend.title = element_blank(),
        axis.title = element_blank())
  
df_all2_row2 <- df_all2[df_all2$Ecoregion %ni% row1, ]
levels(df_all2_row2$facet) <- c("Taxonomic"  , "F.")

aaa2 = ggplot(df_all2_row2, 
              aes(x = Value, y = variable_env, color = dim)) + 
  geom_pointrange(aes(xmin =Value -  std_error, xmax = Value + std_error,
                      fill = dim), shape=21, fatten = 1.1, size = 1,
                  position = position_dodge2(width = 0.35)) + 
  scale_color_manual(values = three_cols, drop = FALSE) +
  scale_fill_manual(values = three_cols, drop = FALSE) +
  theme_classic() + geom_vline(xintercept = 0 ) +
  facet_grid(facet ~ Ecoregion, space='free_y', scales = 'free') + xlim(-1, 1) +
  xlab("Slope") +
  theme(panel.background = element_rect(color = "black", fill = NA),
        panel.grid.major.y = element_line(color = 'grey90'),
        legend.title = element_blank(),
        axis.title.y = element_blank())

aaa3 <- ggarrange(aaa1, aaa2, nrow = 2, 
                  common.legend = TRUE,
                  legend = 'right', heights = c(2,1.2))
aaa3

ggsave(file = 'figures/taxo_fonctio/figure5.jpg', plot = aaa3,
       width = 3,  height = 2, scale = 3)


################################################################################
#######################  save CSV 

df_all$variable_env2 <- paste0(df_all$variable_env,' (', round(df_all$Value, 2),')')

df_all2 <- df_all %>% 
  dplyr::group_by(Ecoregion, dim) %>% 
  dplyr::summarise(test = paste0(variable_env2, collapse = ", "))


df_all2_cast <- reshape2::dcast(df_all2, Ecoregion ~ dim)

write.csv2(df_all2_cast, file = "data/processed/results_gls.csv", row.names = FALSE)
write.csv2(df_all, file = "data/processed/results_gls2.csv", row.names = FALSE)

# I(chloro^2) + I(mlotst^2) + I(oxy_bottom^2)  +
# I(temp_bottom^2) + I(temp_surf^2) +  I(curr_bottom^2)


################################################################################
############################ timeseries abundance species ######################
################################################################################
unique(df_last$genus_sp)[order(unique(df_last$genus_sp))]
unique(df_last$Ecoregion)

# "#F8766D" "Baltic Sea"
# "#F99D1E" "Greater North Sea"
# "#FF61CC" "Celtic Seas"
# "#C77CFF" "Bay of Biscay and the Iberian Coast"
# "#00BFC4" "Western Mediterranean Sea"
# "#74D33A" "Ionian Sea and the Central Mediterranean Sea"
# "#7C4728" "Adriatic Sea"
# "#5D6966" "Aegean-Levantine Sea"

df_limanda <- df_last %>% 
  dplyr::filter(genus_sp == "Limanda_limanda") %>% 
  dplyr::filter(Ecoregion %in% c("Baltic Sea", "Celtic Seas",
                                 "Bay of Biscay and the Iberian Coast"))

df_limanda_small <- df_limanda %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_limanda_small$Ecoregion)


a = ggplot(df_limanda_small, aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  # coord_cartesian(ylim = c(0.2, 2.5)) +
  scale_color_manual(values = c("#F8766D",  "#FF61CC")) +
  ggtitle("Limanda limanda") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank(),
                     plot.title = element_text(face = 'italic'))

a


df_merlu <- df_last %>% 
  dplyr::filter(genus_sp %in% c("Merluccius_merluccius")) %>% 
  dplyr::filter(Ecoregion %in% c( "Western Mediterranean Sea",
                                 "Ionian Sea and the Central Mediterranean Sea",
                                 "Aegean-Levantine Sea"))
df_merlu_small <- df_merlu %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_merlu_small$Ecoregion)
levels(df_merlu_small$Ecoregion) <- c("Aegean-Levantine\nSea",
                                      "Ionian Sea and the\nCentral Med. Sea",
                                      "Western Med. Sea" )
b = ggplot(df_merlu_small, aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  scale_color_manual(values = c( "#5D6966","#74D33A", "#00BFC4"  )) +
  ggtitle("Merluccius merluccius") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank(),
                     plot.title = element_text(face = 'italic'))

b



df_pleuro <- df_last %>% 
  dplyr::filter(genus_sp == "Pleuronectes_platessa") %>% 
  dplyr::filter(Ecoregion %in% c("Western Mediterranean Sea", 
                                 "Celtic Seas", 
                                 "Greater North Sea",
                                 "Adriatic Sea"))

df_pleuro_small <- df_pleuro %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))

c = ggplot(df_pleuro_small, aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  scale_color_manual(values = c( "#F8766D","#FF61CC","#F99D1E"  )) +
  ggtitle("Pleuronectes platessa") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank(),
                     plot.title = element_text(face = 'italic'))

c


df_sprat  <- df_last %>% 
  dplyr::filter(genus_sp == "Sprattus_sprattus") %>% 
  dplyr::filter(Ecoregion %in% c( "Greater North Sea",
                                 "Western Mediterranean Sea"))

df_sprat_small <- df_sprat %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_sprat_small$Ecoregion)
levels(df_sprat_small$Ecoregion) <- c("Greater North Sea",
                                      "Western Med. Sea" )
c_sprat = ggplot(df_sprat_small, aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  scale_color_manual(values = c( "#F99D1E","#00BFC4"  )) +
  ggtitle("Sprattus sprattus") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank(),
                     plot.title = element_text(face = 'italic'))




fig_sX <- ggarrange(a,c,b,c_sprat, nrow = 4)


ggsave(fig_sX,  width = 2,  height = 2.5, scale = 3,
       file = 'figures/taxo_fonctio/timeseries_abundance_taxo.jpg')




df_req  <- df_last %>% 
  dplyr::filter(genus_sp == "Hexanchus_griseus") %>% 
  dplyr::filter(Ecoregion %in% c(  "Bay of Biscay and the Iberian Coast"))
df_req_s <- df_req %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old ))
ggplot(df_req_s, aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  scale_color_manual(values = c( "#F99D1E","#00BFC4"  )) +
  ggtitle("Sprattus sprattus") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank(),
                     plot.title = element_text(face = 'italic'))


df_celtic_biscay <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_tl$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Bay of Biscay and the Iberian Coast"))

df_celtic_biscay_small <- df_celtic_biscay %>% 
  dplyr::group_by(Ecoregion, Year, genus_sp) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 10 & Year != 2018) %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
ggplot(df_celtic_biscay_small,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n), group = Ecoregion),
              fill = 'grey90', color = NA) + facet_wrap(~ genus_sp)+
   geom_point() + geom_path() + geom_smooth(method = 'lm') +
  ggtitle("High TL species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank())

"Lepidorhombus_whiffiagonis"
"Merluccius_merluccius"
"Zeus_faber"
"Galeus_melastomus" 
"Leucoraja_naevus" 

df_sp_tl[df_sp_tl$genus_sp %in% c("Lepidorhombus_whiffiagonis",
                                  "Merluccius_merluccius",
                                  "Leucoraja_naevus"), ]


quantile_age <- quantile(df_sp_traits$length.maturity, 1-0.2)
df_sp_age <-df_sp_traits %>% 
  dplyr::filter(length.maturity >= quantile_age)
dim(df_sp_age)
df_celtic_biscay <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_age$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Bay of Biscay and the Iberian Coast"))

df_celtic_biscay_small <- df_celtic_biscay %>% 
  dplyr::group_by(Ecoregion, Year, genus_sp) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 10 & Year != 2018) %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
ggplot(df_celtic_biscay_small,
       aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n), group = Ecoregion),
              fill = 'grey90', color = NA) + facet_wrap(~ genus_sp)+
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  ggtitle("High TL species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank())
df_sp_age[df_sp_age$genus_sp %in% c("Galeus_melastomus",
                                  "Leucoraja_naevus",
                                  "Zeus_faber"), ]

################################################################################
################################### timeseries abundance traits######################
################################################################################
unique(df_last$genus_sp)[order(unique(df_last$genus_sp))]
unique(df_last$Ecoregion)

# "#F8766D" "Baltic Sea"
# "#F99D1E" "Greater North Sea"
# "#FF61CC" "Celtic Seas"
# "#C77CFF" "Bay of Biscay and the Iberian Coast"
# "#00BFC4" "Western Mediterranean Sea"
# "#74D33A" "Ionian Sea and the Central Mediterranean Sea"
# "#7C4728" "Adriatic Sea"
# "#5D6966" "Aegean-Levantine Sea"


df_sp_traits <- df_last %>% 
  dplyr::group_by(genus_sp, spawning.type,habitat , feeding.mode) %>% 
  dplyr::summarise(length.maturity = mean(length.maturity),
                   age.maturity = mean(age.maturity),
                   growth.coefficient = mean(growth.coefficient),
                   tl = mean(tl))

limit_basse = 0.2

quantile_K1 <- quantile(df_sp_traits$growth.coefficient, limit_basse)
quantile_K2 <- quantile(df_sp_traits$growth.coefficient, 1-limit_basse)
df_sp_K1 <- df_sp_traits %>% 
  dplyr::filter(growth.coefficient >= 0.31)
df_sp_K2 <- df_sp_traits %>% 
  dplyr::filter(growth.coefficient >= 0.52)
dim(df_sp_K2)

df_baltic_north1 <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_K1$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Baltic Sea"))
df_baltic_north2 <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_K2$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Greater North Sea"))
df_baltic_north <- rbind(df_baltic_north1, df_baltic_north2)


df_baltic_north_small <- df_baltic_north %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_baltic_north_small$Ecoregion)
levels(df_baltic_north_small$Ecoregion) <- c( "Baltic Sea","Greater North\nSea")
# 
# df_baltic_north_small$remove <- ifelse(df_baltic_north_small$Ecoregion == "Baltic Sea" &
#                                          df_baltic_north_small$Year %in% c(1998, 2010,2013,2014,2015),
#                                        "yes","no")

a = ggplot(df_baltic_north_small,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) +
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  # coord_cartesian(ylim = c(1.8, 3.7)) +
  scale_color_manual(values = c("#F8766D", "#F99D1E")) +
  ggtitle("Fast-growing species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank())

a



liste_sp_demersal <- unique(df_sp_traits[df_sp_traits$feeding.mode  == 'benthivorous', 
                                         "genus_sp"])$genus_sp

df_baltic_north2 <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_demersal) %>% 
  dplyr::filter(Ecoregion %in% c("Baltic Sea", "Greater North Sea"))

df_baltic_north_small2 <- df_baltic_north2 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 10) %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_baltic_north_small2$Ecoregion)
levels(df_baltic_north_small2$Ecoregion) <- c(  "Baltic Sea",
                                               "Greater North\nSea")

b = ggplot(df_baltic_north_small2,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#F8766D", "#F99D1E")) +
  ggtitle("Benthivorous species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.y= element_blank(),
                     axis.title.x= element_blank())
b





liste_sp_guarder <- unique(df_sp_traits[df_sp_traits$spawning.type  %in% c("guarder"), 
                                        "genus_sp"])$genus_sp

df_baltic_north3 <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_guarder) %>% 
  dplyr::filter(Ecoregion %in% c("Baltic Sea", "Greater North Sea"))

df_baltic_north_small3 <- df_baltic_north3 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_baltic_north_small3$Ecoregion)
levels(df_baltic_north_small3$Ecoregion) <- c( "Baltic Sea",  
                                               "Greater North\nSea")


c = ggplot(df_baltic_north_small3,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#F8766D", "#F99D1E")) +
  ggtitle("Guarder species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.y= element_blank(),
                     axis.title.x= element_blank())

c

pl_north <- ggarrange(a,b,c, ncol = 3, common.legend = TRUE, legend = 'right')
pl_north



############# celtic et biscay 
limit_basse = 0.1

quantile_tl1 <- quantile(df_sp_traits$tl, limit_basse)
quantile_tl2 <- quantile(df_sp_traits$tl, 1-limit_basse)
df_sp_tl <-df_sp_traits %>% 
  dplyr::filter(tl >= quantile_tl2)
dim(df_sp_tl)

df_celtic_biscay <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_tl$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Bay of Biscay and the Iberian Coast"))

df_celtic_biscay_small <- df_celtic_biscay %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 10 & Year != 2018) %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_celtic_biscay_small$Ecoregion)
levels(df_celtic_biscay_small$Ecoregion) <- c("Bay of Biscay\nand the\nIberian Coast")


d = ggplot(df_celtic_biscay_small,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n), group = Ecoregion),
              fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#C77CFF","#FF61CC" )) +
  ggtitle("High TL species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.x= element_blank())

d


# quantile_length2 <- quantile(df_sp_traits$length.maturity, 
#                              1-limit_basse)
# df_sp_le <-df_sp_traits %>% 
#   dplyr::filter(length.maturity >= quantile_length2)
# 
# df_celtic_biscay <- df_last %>% 
#   dplyr::filter(genus_sp %in% unique(df_sp_le$genus_sp)) %>% 
#   dplyr::filter(Ecoregion %in% c("Celtic Seas","Bay of Biscay and the Iberian Coast"))
# df_celtic_biscay_small <- df_celtic_biscay %>% 
#   dplyr::group_by(Ecoregion, Year) %>% 
#   dplyr::summarise(abundance_m = mean(log_abundance),
#                    abundance_sd = sd(log_abundance),
#                    n = n_distinct(index_old )) %>% 
#   # dplyr::filter(n > 4 & Year != 2018) %>%
#   dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
# levels(df_celtic_biscay_small$Ecoregion)
# levels(df_celtic_biscay_small$Ecoregion) <- c("Bay of Biscay\nand the\nIberian Coast", 
#                                               "Celtic Seas")
# 
# dddd = ggplot(df_celtic_biscay_small,
#            aes(x = Year , y = abundance_m , color= Ecoregion)) + 
#   geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
#                   ymax = abundance_m + 1.9*abundance_sd/sqrt(n), group = Ecoregion),
#               fill = 'grey90', color = NA) + 
#   geom_point() + geom_path() + geom_smooth(method = 'lm') +
#   scale_color_manual(values = three_cols) +
#   ggtitle("Hight TL species") +   ylab('Mean abundance') + 
#   theme_bw() + theme(legend.title = element_blank(),
#                      axis.title.x= element_blank())


liste_sp_benthi <- unique(df_sp_traits[df_sp_traits$feeding.mode  == 'benthivorous', 
                                         "genus_sp"])$genus_sp


df_celtic_biscay2 <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_benthi) %>% 
  dplyr::filter(Ecoregion %in% c("Celtic Seas", "Bay of Biscay and the Iberian Coast"))

df_celtic_biscay2 <- df_celtic_biscay2 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>%
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_celtic_biscay2$Ecoregion)
levels(df_celtic_biscay2$Ecoregion) <- c("Bay of Biscay\nand the\nIberian Coast",
                                         "Celtic Seas")

e = ggplot(df_celtic_biscay2,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#C77CFF","#FF61CC" )) +
  ggtitle("Benthivorous species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.y= element_blank(),
                     axis.title.x= element_blank())
e




liste_sp_non_guarder <- unique(df_sp_traits[df_sp_traits$spawning.type  == 'non-guarder', 
                                           "genus_sp"])$genus_sp
df_celtic_biscay3 <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_non_guarder) %>% 
  dplyr::filter(Ecoregion %in% c("Celtic Seas"))

df_celtic_biscay3 <- df_celtic_biscay3 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 10)

f = ggplot(df_celtic_biscay3,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#FF61CC" )) +
  ggtitle("Non-guarder species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.y= element_blank(),
                     axis.title.x= element_blank())


# "#F8766D" "Baltic Sea"
# "#F99D1E" "Greater North Sea"
# "#FF61CC" "Celtic Seas"
# "#C77CFF" "Bay of Biscay and the Iberian Coast"
# "#00BFC4" "Western Mediterranean Sea"
# "#74D33A" "Ionian Sea and the Central Mediterranean Sea"
# "#7C4728" "Adriatic Sea"
# "#5D6966" "Aegean-Levantine Sea"


liste_sp_guarder <- unique(df_sp_traits[df_sp_traits$spawning.type  == 'guarder', 
                                            "genus_sp"])$genus_sp
df_celtic_biscay4 <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_guarder) %>% 
  dplyr::filter(Ecoregion %in% c("Bay of Biscay and the Iberian Coast"))

df_celtic_biscay4 <- df_celtic_biscay4 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) 

g = ggplot(df_celtic_biscay4,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion), fill = 'grey90', color = NA) + 
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#C77CFF" )) +
  ggtitle("Guarder species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title.y= element_blank(),
                     axis.title.x= element_blank(),
                     legend.position = "none")


pl_west <- ggarrange(d,e,f,g, ncol = 4, common.legend = TRUE, legend = 'right')
pl_west

############# med

df_ad_med2 <- df_last %>% 
  dplyr::filter(genus_sp %in% unique(df_sp_K1$genus_sp)) %>% 
  dplyr::filter(Ecoregion %in% c("Western Mediterranean Sea"))

df_ad_med_small2 <- df_ad_med2 %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 1) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
# levels(df_ad_med_small2$Ecoregion)
levels(df_ad_med_small2$Ecoregion) <- c("Western\nMed. Sea")


h = ggplot(df_ad_med_small2,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion),
              fill = 'grey90', color = NA) + 
  geom_point() + geom_path() +  geom_smooth(method = 'lm') +
  # coord_cartesian(ylim = c(0.2, 2.5)) +
  scale_color_manual(values = c("#00BFC4")) +
  ggtitle("Fast-growing species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = 'none',
                     axis.title.x= element_blank())
h



liste_sp_non_guarder <- unique(df_sp_traits[df_sp_traits$spawning.type  == 'non-guarder', 
                                            "genus_sp"])$genus_sp
df_west_med <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_non_guarder) %>% 
  dplyr::filter(Ecoregion %in% c("Western Mediterranean Sea"))

df_west_med_small <- df_west_med %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 5) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))
levels(df_west_med_small$Ecoregion) <- c("Western\nMed. Sea")

i = ggplot(df_west_med_small,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion),fill = 'grey90', color = NA) +
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#00BFC4")) +
  coord_cartesian() +
  ggtitle("Non-guarder species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = 'none',
                     axis.title= element_blank())


df_adri <- df_last %>% 
  dplyr::filter(genus_sp %in% liste_sp_non_guarder) %>% 
  dplyr::filter(Ecoregion %in% c("Adriatic Sea") &
                  feeding.mode == "benthivorous" )

df_adri_small <- df_adri %>% 
  dplyr::group_by(Ecoregion, Year) %>% 
  dplyr::summarise(abundance_m = mean(log_abundance),
                   abundance_sd = sd(log_abundance),
                   n = n_distinct(index_old )) %>% 
  dplyr::filter(n > 5) %>% 
  dplyr::mutate(Ecoregion = factor(Ecoregion, exclude = T))

i_adri = ggplot(df_adri_small,
           aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion),fill = 'grey90', color = NA) +
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#7C4728")) +
  coord_cartesian() +
  ggtitle("Benthivorous species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     legend.position = 'none',
                     axis.title= element_blank())



pl_med <- ggarrange(h,i,i_adri,ncol = 3, widths = c(1,1,1),
                    common.legend = TRUE, legend = 'right')



pl_all <- ggarrange(pl_north,pl_west,pl_med, nrow = 3, common.legend = FALSE)
pl_all


ggsave(pl_all,  width = 3,  height = 2.4, scale = 3,
       file = 'figures/taxo_fonctio/timeseries_abundance_.jpg')


i_adri = ggplot(df_adri_small,
                aes(x = Year , y = abundance_m , color= Ecoregion)) + 
  geom_ribbon(aes(x = Year, ymin = abundance_m - 1.9*abundance_sd/sqrt(n),
                  ymax = abundance_m + 1.9*abundance_sd/sqrt(n),
                  group = Ecoregion),fill = 'grey90', color = NA) +
  geom_point() + geom_path() + geom_smooth(method = 'lm') +
  scale_color_manual(values = c("#7C4728")) +
  coord_cartesian() +
  ggtitle("Benthivorous species") +   ylab('Mean abundance') + 
  theme_bw() + theme(legend.title = element_blank(),
                     axis.title= element_blank())
ggsave(i_adri,  width = 1.2,  height = 1, scale = 3,
       file = 'figures/taxo_fonctio/timeseries_abundance_forlegend.jpg')






################################################################################
################## test effet chalut donnees initiales ##################################
################################################################################

# https://pmarchand1.github.io/ECL7102/notes_cours/14-Analyses_multivariees_Partie2.pdf



df_number_cell_yy_qq <- df_last %>% 
  dplyr::group_by(Year, Quarter, my_spatial_id, yy_qq_carre) %>% 
  dplyr::summarise(n_survey = n_distinct(Survey)) %>% 
  dplyr::filter(n_survey >= 2) 
table(df_number_cell_yy_qq$n_survey)

df_last_text_survey <- df_last %>% 
  dplyr::filter(yy_qq_carre %in% unique(df_number_cell_yy_qq$yy_qq_carre))


trawl_sp_matrice <- reshape2::acast(data = df_last_text_survey,
                                    index_old ~ genus_sp,
                                    value.var = 'log_abundance',
                                    fun.aggregate = mean)
trawl_sp_matrice[is.na(trawl_sp_matrice)] <- 0
trawl_sp_matrice_hell <- decostand(trawl_sp_matrice, method = 'hellinger')

res_acp <- PCA(trawl_sp_matrice_hell, scale.unit = FALSE, ncp = 50)

df_coord_chaluts <- data.frame(res_acp$ind$coord)[,c(1:13)]
df_coord_chaluts$index_old <- rownames(df_coord_chaluts)

df_complementary <- data.frame(df_last_text_survey) %>%
  dplyr::group_by(index_old, Year, Quarter, Lon, Lat,depth, Survey) %>%
  dplyr::summarise(n_sp = n_distinct(genus_sp)) %>% 
  dplyr::select(-n_sp)

trawl_sp_matrice_hell2 <- data.frame(trawl_sp_matrice_hell)
trawl_sp_matrice_hell2$index_old <- rownames(trawl_sp_matrice_hell2)

df_coord_chaluts_PCA <- merge(df_coord_chaluts, df_complementary)
df_coord_chaluts_raw <- merge(trawl_sp_matrice_hell2, df_complementary)

## methode 1 = PERMANOVA 

res_ado <- adonis2(trawl_sp_matrice ~ Year + Quarter + Lon + Lat +  depth + Survey, 
                   df_complementary, sqrt.dist = TRUE)


## methode 2 = ACP puis projection des var 
mat_dist <- dist(trawl_sp_matrice_hell)
pcoa_hell <- wcmdscale(mat_dist, eig = TRUE, add = TRUE)
envf <- envfit(pcoa_hell, select(data.frame(df_complementary), Year , Quarter , Lon ,
                                 Lat ,  depth , Survey))
envf


## methode 3 = RDA : methode sous contrainte, mais basée sur distance euclidienne ...

my_rda <- rda(trawl_sp_matrice_hell ~ Year + Quarter + Lon + Lat +  depth + Survey,
               data = df_complementary)
anova.cca(my_rda)
RsquareAdj(my_rda)  # Total variance explained by the RDA
anova.cca(my_rda, by = "terms")  # Test which terms are significant


## methode 4 = db-RDA 

res_acp <- PCA(trawl_sp_matrice_hell, scale.unit = FALSE, ncp = 50)

my_dbrda <- rda(res_acp$ind$coord ~ Year + Quarter + Lon + Lat +  depth + Survey,
              data = df_complementary)
anova.cca(my_dbrda, by = "terms")  # Test which terms are significant


# or

dbRDA = capscale(trawl_sp_matrice ~ Year + Quarter + Lon + Lat +  depth + Survey,
                 df_complementary, dist="bray", sqrt.dist=TRUE)
plot(dbRDA) # use base plot, might be done with ggplot2
anova.cca(dbRDA, by = "terms")  # Test which terms are significant


#1 a distance matrix is calculated using the distance measure of choice. 
#2 a principle coordinates analysis (PCoA) is done on the matrix. 
#3 th eigenvalues obtained in the PCoA are plugged into an RDA. 


################################################################################
################################### prc ##################################
################################################################################


if(!require(vegan)){install.packages("vegan"); library(vegan)}
if(!require(quantreg)){install.packages("quantreg"); library(quantreg)}

df_Q1 <- df_last[df_last$Quarter == 1, ]

nbofyears<-length(levels(as.factor(df_Q1$Year)))
nbofsites<-length(levels(as.factor(df_Q1$my_spatial_id)))

siteslabels<-seq(1,nbofsites,1)
Year<-df_Q1$Year

sitesfunction<-gl(nbofsites,nbofyears,labels=siteslabels)
time_period<-df_Q1$Year

time_period<-'1994-2019'

# cleaning table of NA items
species<-df_Q1[,'m_abundance']
length(species)

time<-as.factor(df_Q1$Year)
length(time)


sites<-as.factor(df_Q1$my_spatial_id)
length(sites)


# PRC
PRCresults <- prc(species,time,sites)
PRCresults

PRCresults
plot(PRCresults,ylim=c(-1,1))
# test for PRC significance
testPRC<-anova.cca(PRCresults,step=1000)
testPRC

################################################################################
################################### statis ##################################
################################################################################

#################### taxo #################### 
head(df_pca_taxo)
df_pca_taxo$Year <- as.num(df_pca_taxo$Year)
df_pca_taxo$Quarter <- as.num(df_pca_taxo$Quarter)

df_statis_taxo_response <- df_pca_taxo[,c( "Dim.1","Dim.2","Dim.3","Dim.4","Dim.5")]
df_statis_taxo_variables_bloc1 <- df_pca_taxo[,c( "Year","Quarter")]
df_statis_taxo_variables_bloc2 <- df_pca_taxo[,c("x_my_spatial_id", "y_my_spatial_id")]
df_statis_taxo_variables_bloc3 <- df_pca_taxo[,c( "chloro_mea","mlotst_mea",
                                                  "oxy_bottom_mea","oxy_surf_mea",
                                                  "temp_bottom_mea","temp_surf_mea",
                                                  "curr_bottom_mea","curr_surf_mea")]

ktabX_chik <- list(df_statis_taxo_variables_bloc1 ,
                   df_statis_taxo_variables_bloc2, 
                   df_statis_taxo_variables_bloc3)
names(ktabX_chik) <- c('temporel','spatial','env')




dudi1 <- dudi.pca(df_statis_taxo_response, scann = FALSE, scal = FALSE)


wit1_spatial <- wca(dudi1, as.factor(df_pca_taxo$my_spatial_id), scann = FALSE)
kta3_spatial <- ktab.within(wit1_spatial)
statis3_spatial <- statis(kta3_spatial, scann = FALSE)
RV_spatial <- statis3_spatial$RV
diag(RV_spatial)=NA

mean(RV_spatial, na.rm = T)
sd(RV_spatial, na.rm = T)


# RV.rtest(data.frame(RV_spatial))
# coeffRV(RV_spatial)



wit1_year <- wca(dudi1, as.factor(df_pca_taxo$Year), scann = FALSE)
kta3_year  <- ktab.within(wit1_year)
statis3_year  <- statis(kta3_year, scann = FALSE)
RV_year  <- statis3_year$RV
diag( )=NA
rv_yy_taxo <- mean(RV_year, na.rm = T)
mean(RV_year, na.rm = T)
sd(RV_year, na.rm = T)

wit1_quarter <- wca(dudi1, as.factor(df_pca_taxo$Quarter), scann = FALSE)
kta3_quarter  <- ktab.within(wit1_quarter)
statis3_quarter  <- statis(kta3_quarter, scann = FALSE)
RV_quarter  <- statis3_quarter$RV
RV_quarter2 <- RV_quarter[RV_quarter != 1]
rv_qq_taxo <- mean(RV_quarter2)

# 
# if(adegraphicsLoaded()) {
#   s.arrow(statis3_year$C.li, pgrid.text.cex = 0)
#   kplot(statis3_year, traj = TRUE, arrow = FALSE, plab.cex = 0, psub.cex = 3, ppoi.cex = 3)
# } else {
#   s.arrow(statis3_year$C.li, cgrid = 0)
#   kplot(statis3_year, traj = TRUE, arrow = FALSE, unique = TRUE, 
#         clab = 0, csub = 3, cpoi = 3)
# }
# 
# statis3_year


liste_cells <- unique(df_pca_taxo$my_spatial_id)

df_spatial_id <- data.frame(my_spatial_id = rep(liste_cells, each = 8))
df_spatial_id$voisin <- NA

for(rep_cell in liste_cells){
  print(rep_cell)
  # rep_cell = "AC3"
  let1 <- substr(rep_cell, 1, 1); position_let1 = which(LETTERS == let1)
  let2 <- substr(rep_cell, 2, 2); position_let2 = which(LETTERS == let2)
  num <- substr(rep_cell, 3, 4)
  
    rep_gauche <- paste0(let1, let2,as.num(num)-1 )
    rep_droite <- paste0(let1, let2,as.num(num)+1 )
    
  if(let2 %ni% c('A', "Z")){
  rep_haut <- paste0(let1, LETTERS[position_let2-1], as.num(num) )
  rep_bas <- paste0(let1, LETTERS[position_let2+1], as.num(num) )
  
  rep_haut_gauche <- paste0(let1, LETTERS[position_let2-1], as.num(num)-1 )
  rep_haut_droite <-  paste0(let1, LETTERS[position_let2-1], as.num(num)+1 )
  
  rep_bas_gauche <-paste0(let1, LETTERS[position_let2+1], as.num(num)-1 )
  rep_bas_droite <- paste0(let1, LETTERS[position_let2+1], as.num(num)+1 )
    
  } else if(let2 == "Z"){
    rep_haut <- paste0(let1, LETTERS[position_let2-1], as.num(num) )
    rep_bas <- paste0(LETTERS[position_let1+1],'A', as.num(num) )
    
    rep_haut_gauche <- paste0(let1, LETTERS[position_let2-1], as.num(num) - 1)
    rep_haut_droite <- paste0(let1, LETTERS[position_let2-1], as.num(num) + 1)
    
    rep_bas_gauche <- paste0(LETTERS[position_let1+1],'A', as.num(num) - 1)
    rep_bas_droite <- paste0(LETTERS[position_let1+1],'A', as.num(num) + 1)
    
    
  } else if(let2 == "A"){
    rep_haut <- paste0(LETTERS[position_let1-1], "Z", as.num(num) )
    rep_bas <- paste0(let1, LETTERS[position_let2+1], as.num(num) )
    
    rep_haut_gauche <- paste0(LETTERS[position_let1-1], "Z", as.num(num) - 1)
    rep_haut_droite <- paste0(LETTERS[position_let1-1], "Z", as.num(num) + 1)
    
    rep_bas_gauche <- paste0(let1, LETTERS[position_let2+1], as.num(num) - 1)
    rep_bas_droite <- paste0(let1, LETTERS[position_let2+1], as.num(num) + 1)
    
  }
    
    
    
    df_spatial_id[df_spatial_id$my_spatial_id == rep_cell, 'voisin'] <- c(rep_haut_gauche,
                                                                          rep_haut, 
                                                                          rep_haut_droite,
                                                                          rep_gauche,
                                                                          rep_droite,
                                                                           rep_bas_gauche,
                                                                          rep_bas,
                                                                          rep_bas_droite)
    
  }
  
  
df_RV_spatial <- as.data.frame(reshape2::melt(RV_spatial))
names(df_RV_spatial)<-c("my_spatial_id", "voisin" , "RV")

dede <- merge(df_spatial_id, df_RV_spatial, all.x = TRUE)
dede2 <- dede %>% 
  dplyr::group_by(my_spatial_id) %>% 
  dplyr::summarise(RV_m = mean(RV, na.rm = T))

df_coords <- df_pca_taxo %>% 
  dplyr::group_by(my_spatial_id,x_my_spatial_id,y_my_spatial_id) %>% 
  dplyr::summarise(y_m = mean(Year))
dede2 <- merge(dede2,df_coords )

mapRV <- ggplot(data = world) +  
  geom_sf(color = NA, fill = 'grey80') +  theme_classic() + 
  
  geom_tile(data = dede2,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill =1- RV_m  ),
            size = 1) + scale_fill_viridis_c(name = 'RV') + 
  coord_sf(xlim = c(-20, 30), ylim = c( 34.5 , 65),expand = FALSE) + 
  guides(fill = guide_colorbar(barwidth = 1, barheight = 8))  +
  theme(axis.title = element_blank())
mapRV
ggsave(file = 'figures/taxo_fonctio/mapRV.jpg', plot = mapRV,
       width = 3,  height = 2, scale = 3)


#################### fonctio #################### 

head(df_pca_fonctio)
df_pca_fonctio$Year <- as.num(df_pca_fonctio$Year)
df_pca_fonctio$Quarter <- as.num(df_pca_fonctio$Quarter)

df_statis_func_response <- df_pca_fonctio[,c( "Dim.1","Dim.2","Dim.3","Dim.4","Dim.5")]



dudi1_f <- dudi.pca(df_statis_func_response, scann = FALSE, scal = FALSE)


wit1_spatial_f <- wca(dudi1_f, as.factor(df_pca_fonctio$my_spatial_id), scann = FALSE)
kta3_spatial_f <- ktab.within(wit1_spatial_f)
statis3_spatial_f <- statis(kta3_spatial_f, scann = FALSE)
RV_spatial_f <- statis3_spatial_f$RV
diag(RV_spatial_f)=NA

mean(RV_spatial_f, na.rm = T)
sd(RV_spatial_f, na.rm = T)



wit1_year_f <- wca(dudi1_f, as.factor(df_pca_fonctio$Year), scann = FALSE)
kta3_year_f  <- ktab.within(wit1_year_f)
statis3_year_f  <- statis(kta3_year_f, scann = FALSE)
RV_year_f  <- statis3_year_f$RV
diag(RV_year_f)=NA

mean(RV_year_f, na.rm = T)
sd(RV_year_f, na.rm = T)



wit1_quarter <- wca(dudi1, as.factor(df_pca_fonctio$Quarter), scann = FALSE)
kta3_quarter  <- ktab.within(wit1_quarter)
statis3_quarter  <- statis(kta3_quarter, scann = FALSE)
RV_quarter  <- statis3_quarter$RV
RV_quarter2 <- RV_quarter[RV_quarter != 1]
rv_qq_functio <- mean(RV_quarter2)



mean(RV_spatial, na.rm = T)
sd(RV_spatial, na.rm = T)

mean(RV_spatial_f, na.rm = T)
sd(RV_spatial_f, na.rm = T)


mean(RV_year, na.rm = T)
sd(RV_year, na.rm = T)

mean(RV_year_f, na.rm = T)
sd(RV_year_f, na.rm = T)



df_statis <- data.frame(facet = c("Taxonomy", "Functional"), 
                         RV_spatial = NA, RV_yearly = NA)
df_statis[df_statis$facet == "Taxonomy", "RV_spatial"] <- round(rv_spatial_taxo, 2)
df_statis[df_statis$facet == "Taxonomy", "RV_yearly"] <- round(rv_yy_taxo, 2)
df_statis[df_statis$facet == "Functional", "RV_spatial"] <- round(rv_spatial_functio, 2)
df_statis[df_statis$facet == "Functional", "RV_yearly"] <- round(rv_yy_functio, 2)

write.csv2(file = "data/processed/df_statis.csv", df_statis, row.names = FALSE)

################################################################################
######################### relations taxo fonctio ###############################
################################################################################

head(df_for_pcoa2)
df_pca_taxo_melt <- reshape2::melt(data = df_pca_taxo, 
                                    id.vars = c("index_unique"),
                                    measure.vars = c("Dim.1","Dim.2","Dim.3","Dim.4", "Dim.5"),
                                   value.name = 'Taxonomy',
                                   variable.name = 'Dimension')
df_pca_taxo_melt$Dimension <- substr(df_pca_taxo_melt$Dimension, 5, 5)


df_pca_func_melt <- reshape2::melt(data = df_for_pcoa2, 
                                   id.vars = c("index_unique"),
                                   measure.vars = c("pcoa1","pcoa2","pcoa3","pcoa4", "pcoa5"),
                                   value.name = 'Functional',
                                   variable.name = 'Dimension')
df_pca_func_melt$Dimension <- substr(df_pca_func_melt$Dimension, 5, 5)


df_both <- merge(df_pca_func_melt, df_pca_taxo_melt)

ggplot(df_both, aes(x = Taxonomy, y = Functional)) + 
  geom_point(alpha = 0.1) + geom_smooth(method = 'lm') + 
  # facet_wrap(~ Dimension) + 
  theme_bw()

mylm <- lm(data = df_both, pcoa1 ~ Dim.1)
summary(mylm)


################################################################################
################################### var env ##################################
################################################################################

####### SPatial taxo

df_pca_taxo2 <- df_pca_taxo %>% 
  dplyr::group_by(index_unique,  Year, Quarter, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   depth  =  mean(depth, na.rm = T),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   oxy_surf_mea  =  mean(oxy_surf_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   temp_surf_mea  =  mean(temp_surf_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T),
                   curr_surf_mea  =  mean(curr_surf_mea, na.rm = T))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pco2_m),
                   depth  =  mean(depth, na.rm = T),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T)) %>% 
  dplyr::filter(!is.na(depth) & !is.na(chloro_mea) & !is.na(mlotst_mea) &
                  !is.na(oxy_bottom_mea) & !is.na(temp_bottom_mea) & !is.na(curr_bottom_mea)  )

df_pca_taxo3[,c("depth" , "chloro_mea" , "mlotst_mea" ,
                "oxy_bottom_mea" , "temp_bottom_mea" , 
                "curr_bottom_mea")] <- scale(df_pca_taxo3[,c("depth" , "chloro_mea" , 
                                                             "mlotst_mea" , "oxy_bottom_mea" ,
                                                             "temp_bottom_mea" ,  "curr_bottom_mea")])



my_rda2 <- rda(df_pca_taxo3[,c("pcoa1_m",'pcoa2_m')] ~ depth + chloro_mea +
                 mlotst_mea + oxy_bottom_mea + temp_bottom_mea +  curr_bottom_mea,
               data = df_pca_taxo3)
anova.cca(my_rda2)
RsquareAdj(my_rda2)  # Total variance explained by the RDA
anova.cca(my_rda2, by = "axis")  # Test which axis are significant

anova.cca(my_rda2, by = "terms")  # Test which terms are significant

sqrt(vif.cca(my_rda2))

global_r2 <- RsquareAdj(my_rda2)$adj.r.squared

env <- as.matrix(df_pca_taxo3[,c("curr_bottom_mea","temp_bottom_mea","mlotst_mea","depth", 
                                "chloro_mea","oxy_bottom_mea")])

fs <- forward.sel(df_pca_taxo3[,c("pcoa1_m",'pcoa2_m')], # Y matrix
                  env, # X matrix
                  adjR2thresh = global_r2, # Set the adj.R2 threshold
                  alpha = 0.001, # Set alpha level
                  nperm = 999  )   
df_taxo_spatial  <- fs[,c(1,2,3)]
df_taxo_spatial$pattern <- 'Spatial'
df_taxo_spatial$facet <- "Taxonomy"
df_taxo_spatial


####### SPatial fonctio

df_for_pcoa3 <- df_for_pcoa2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pcoa1),
                   pcoa2_m  =  mean(pcoa1),
                   depth  =  mean(depth, na.rm = T),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T))%>% 
  dplyr::filter(!is.na(depth) & !is.na(chloro_mea) & !is.na(mlotst_mea) &
                  !is.na(oxy_bottom_mea) & !is.na(temp_bottom_mea) & !is.na(curr_bottom_mea)  )

df_for_pcoa3[,c("depth" , "chloro_mea" , "mlotst_mea" ,
                "oxy_bottom_mea" , "temp_bottom_mea" , 
                "curr_bottom_mea")] <- scale(df_for_pcoa3[,c("depth" , "chloro_mea" , 
                                                             "mlotst_mea" , "oxy_bottom_mea" ,
                                                             "temp_bottom_mea" ,  "curr_bottom_mea")])

my_rda_func <- rda(df_for_pcoa3[,c("pcoa1_m",'pcoa2_m')] ~ depth + chloro_mea +
                 mlotst_mea + oxy_bottom_mea + temp_bottom_mea +  curr_bottom_mea,
               data = df_for_pcoa3)
anova.cca(my_rda_func)
RsquareAdj(my_rda_func)  # Total variance explained by the RDA
anova.cca(my_rda_func, by = "axis")  # Test which axis are significant

anova.cca(my_rda_func, by = "terms")  # Test which terms are significant

sqrt(vif.cca(my_rda_func))

global_r2_func <- RsquareAdj(my_rda_func)$adj.r.squared

env_func <- as.matrix(df_for_pcoa3[,c("curr_bottom_mea","temp_bottom_mea",
                                 "chloro_mea","oxy_bottom_mea","depth", "mlotst_mea")])

fs_func <- forward.sel(df_for_pcoa3[,c("pcoa1_m",'pcoa2_m')], # Y matrix
                       env_func, # X matrix
                  adjR2thresh = global_r2_func, # Set the adj.R2 threshold
                  alpha = 0.001, # Set alpha level
                  nperm = 999) 
df_fonctio_spatial  <- fs_func[,c(1,2,3)]
df_fonctio_spatial$pattern <- 'Spatial'
df_fonctio_spatial$facet <- "Functional"
df_fonctio_spatial


####### temporel taxo

df_pca_taxo2 <- df_pca_taxo %>% 
  dplyr::group_by(index_unique,  Year, Quarter, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   depth  =  mean(depth, na.rm = T),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   oxy_surf_mea  =  mean(oxy_surf_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   temp_surf_mea  =  mean(temp_surf_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T),
                   curr_surf_mea  =  mean(curr_surf_mea, na.rm = T))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pco2_m),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T)) %>% 
  dplyr::filter(!is.na(chloro_mea) & !is.na(mlotst_mea) &
                  !is.na(oxy_bottom_mea) & !is.na(temp_bottom_mea) & !is.na(curr_bottom_mea)  )
df_pca_taxo3[,c( "chloro_mea" , "mlotst_mea" ,
                "oxy_bottom_mea" , "temp_bottom_mea" , 
                "curr_bottom_mea")] <- scale(df_pca_taxo3[,c( "chloro_mea" , 
                                                             "mlotst_mea" , "oxy_bottom_mea" ,
                                                             "temp_bottom_mea" ,  "curr_bottom_mea")])

my_rda3 <- rda(df_pca_taxo3[,c("pcoa1_m",'pcoa2_m')] ~ chloro_mea +
                 mlotst_mea + oxy_bottom_mea + temp_bottom_mea +  curr_bottom_mea,
               data = df_pca_taxo3)
anova.cca(my_rda3)
RsquareAdj(my_rda3)  # Total variance explained by the RDA
anova.cca(my_rda3, by = "axis")  # Test which axis are significant

anova.cca(my_rda3, by = "terms")  # Test which terms are significant

sqrt(vif.cca(my_rda3))

global_r3 <- RsquareAdj(my_rda3)$r.squared

env <- as.matrix(df_pca_taxo3[,c( "temp_bottom_mea",
                                    "chloro_mea","mlotst_mea" ,
                                    "curr_bottom_mea",  "oxy_bottom_mea")])

fs3 <- forward.sel(df_pca_taxo3[,c("pcoa1_m",'pcoa2_m')], # Y matrix
                  env, # X matrix
                  adjR2thresh = 0.99, # Set the adj.R2 threshold
                  alpha = 0.001, # Set alpha level
                  nperm = 999  )   
df_taxo_tempo <- fs3[,c(1,2,3)]
df_taxo_tempo$pattern <- 'Temporal'
df_taxo_tempo$facet <- "Taxonomy"
df_taxo_tempo

####### temporel fonctio

df_for_pcoa3 <- data.frame(df_for_pcoa2) %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(pcoa1_m = mean(pcoa1),
                   pcoa2_m  =  mean(pcoa1),
                   chloro_mea  =  mean(chloro_mea, na.rm = T),
                   mlotst_mea  =  mean(mlotst_mea, na.rm = T),
                   oxy_bottom_mea  =  mean(oxy_bottom_mea, na.rm = T),
                   temp_bottom_mea  =  mean(temp_bottom_mea, na.rm = T),
                   curr_bottom_mea  =  mean(curr_bottom_mea, na.rm = T))%>% 
  dplyr::filter(!is.na(chloro_mea) & !is.na(mlotst_mea) &
                  !is.na(oxy_bottom_mea) & !is.na(temp_bottom_mea) & !is.na(curr_bottom_mea)  )

df_for_pcoa3[,c( "chloro_mea" , "mlotst_mea" ,
                 "oxy_bottom_mea" , "temp_bottom_mea" , 
                 "curr_bottom_mea")] <- scale(df_for_pcoa3[,c( "chloro_mea" , 
                                                               "mlotst_mea" , "oxy_bottom_mea" ,
                                                               "temp_bottom_mea" ,  "curr_bottom_mea")])

  
  my_rda3 <- rda(df_for_pcoa3[,c("pcoa1_m",'pcoa2_m')] ~ chloro_mea +
                   mlotst_mea + oxy_bottom_mea + temp_bottom_mea +  curr_bottom_mea,
                 data = df_for_pcoa3)
  anova.cca(my_rda3)
  RsquareAdj(my_rda3)  # Total variance explained by the RDA
  anova.cca(my_rda3, by = "axis")  # Test which axis are significant
  
  anova.cca(my_rda3, by = "terms")  # Test which terms are significant
  
  sqrt(vif.cca(my_rda3))
  
  global_r3 <- RsquareAdj(my_rda3)$adj.r.squared
  
  env <- as.matrix(df_for_pcoa3[,c( "temp_bottom_mea", "curr_bottom_mea",
                                      "oxy_bottom_mea","mlotst_mea", "chloro_mea")])
  
  fs3 <- forward.sel(df_for_pcoa3[,c("pcoa1_m",'pcoa2_m')], # Y matrix
                     env, # X matrix
                     adjR2thresh = 0.9, # Set the adj.R2 threshold
                     alpha = 0.1, # Set alpha level
                     nperm = 999)   
  df_func_tempo <- fs3[,c(1,2,3)]
  
  df_func_tempo$pattern <- 'Temporal'
  df_func_tempo$facet <- "Functional"
  df_func_tempo



####### Figure

df_result_rda <- rbind(df_taxo_spatial[1:4,], df_fonctio_spatial, 
                       df_taxo_tempo, df_func_tempo)


df_result_rda$variables <- as.factor(df_result_rda$variables)
levels(df_result_rda$variables)
levels(df_result_rda$variables) <- c('Chlorophyll', "Bottom current",
                                     'Depth', 'Mixed layer depth', 'Bottom oxygen',
                                     'Bottom temperature')


a = ggplot(df_result_rda[df_result_rda$variables != 'Chlorophyll', ], 
           aes(x = variables, y = R2*100)) + 
  geom_bar(stat = 'identity' )+ facet_grid(pattern~ facet) + theme_classic() + 
  coord_flip() + xlab("") + ylab("Deviance explained") + 
  theme(panel.background = element_rect(fill = "white", color = 'black'),
        panel.grid.major.y = element_line(color = 'grey90'))
a
ggsave(file = 'figures/taxo_fonctio/figure4.jpg', plot = a,
       width = 2,  height = 1, scale = 3)



################################################################################
################################ by subregion ##################################
################################################################################
df_last[df_last$Ecoregion == "Oceanic Northeast Atlantic", 'Ecoregion'] <- "Celtic Seas"

df_last$Ecoregion_quarter <- paste0(df_last$Ecoregion, '_', df_last$Quarter)
liste_region_quarter <- unique(df_last$Ecoregion_quarter)

for(rep_region_quarter in liste_region_quarter){
  print(rep_region_quarter)
  
  temp <- df_last %>% dplyr::filter(Ecoregion_quarter == rep_region_quarter)
  
  if(n_distinct(temp$my_spatial_id) >= 10){
    trawl_sp_matrice <- reshape2::acast(data = temp,
                                    index_unique ~ genus_sp,
                                    value.var = 'log_abundance',
                                    fun.aggregate = mean)
trawl_sp_matrice[is.na(trawl_sp_matrice)] <- 0
trawl_sp_matrice_hell <- decostand(trawl_sp_matrice, method = 'hellinger')

df_complementary <- temp %>%
  dplyr::group_by(index_unique, Year, Quarter, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>%
  dplyr::summarise(n_sp = n_distinct(genus_sp),
                   abundance_tot = sum(m_abundance),
                   depth = mean(depth, na.rm = TRUE),
                   chloro_mea = mean(chloro_mea, na.rm = TRUE),
                   mlotst_mea = mean(mlotst_mea, na.rm = TRUE),
                   oxy_bottom_mea = mean(oxy_bottom_mea, na.rm = TRUE),
                   oxy_surf_mea = mean(oxy_surf_mea, na.rm = TRUE),
                   temp_bottom_mea = mean(temp_bottom_mea, na.rm = TRUE), 
                   temp_surf_mea = mean(temp_surf_mea, na.rm = TRUE),
                   curr_surf_mea   = mean(curr_surf_mea, na.rm = TRUE),
                   curr_bottom_mea = mean(curr_bottom_mea, na.rm = TRUE))

df_complementary_small1 <- data.frame(df_complementary) %>%
  dplyr::select( Year, Quarter, x_my_spatial_id, y_my_spatial_id,
                 depth,  11:18) %>%
  dplyr::mutate(lon_ICES = x_my_spatial_id, 
                lat_ICES = y_my_spatial_id) %>%
  dplyr::select(- x_my_spatial_id, -y_my_spatial_id)

df_complementary_small1_scalled <- scale(df_complementary_small1)


trawl_sp_matrice2 <- cbind(trawl_sp_matrice_hell, df_complementary_small1_scalled)
trawl_sp_matrice2 <- cbind(trawl_sp_matrice2, data.frame(df_complementary[,'Ecoregion']))


names(trawl_sp_matrice2)
n_column_taxo <- length(names(trawl_sp_matrice2))
iii <- which(names(trawl_sp_matrice2) == "Year")

res_acp_taxo <- PCA(trawl_sp_matrice2, scale.unit = FALSE, ncp = 50,
                    quanti.sup = c(iii:(n_column_taxo-1)),
                    quali.sup = n_column_taxo)

############ results 

df_pca_taxo <-  data.frame(res_acp_taxo$ind$coord)
df_pca_taxo$index_unique <- rownames(df_pca_taxo)

df_pca_taxo <- merge(df_complementary, df_pca_taxo)


df_pca_taxo2 <- df_pca_taxo %>% 
  dplyr::group_by(index_unique,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),
                   pco3_m  =  mean(Dim.3),
                   pco4_m  =  mean(Dim.4),
                   pco5_m  =  mean(Dim.5))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pco2_m),
                   pcoa3_m  =  mean(pco3_m),
                   pcoa4_m  =  mean(pco4_m),
                   pcoa5_m  =  mean(pco5_m))
min_lat <- min(df_pca_taxo3$y_my_spatial_id) - 2
max_lat <- max(df_pca_taxo3$y_my_spatial_id) + 2

min_lon <- min(df_pca_taxo3$x_my_spatial_id) - 2
max_lon <- max(df_pca_taxo3$x_my_spatial_id) + 2


maptaxo <- ggplot(data = world) +  
  geom_tile(data = df_pca_taxo3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
            size = 1) + scale_fill_viridis_c(name = 'PC1') + 
  geom_sf(color = NA, fill = 'grey80') +  theme_classic() + ggtitle('Taxonomy') +
  coord_sf(xlim = c(min_lon, max_lon), ylim = c( min_lat , max_lat),expand = FALSE) + 
  # guides(fill = guide_colorbar(barwidth = 6, barheight = 0.5)) +
  theme(axis.title = element_blank(),
        legend.position = 'right',
        plot.title = element_text(face = 'bold')) 

coord_sp <- as.data.frame(res_acp_taxo$var$coord)
coord_sp$sp <- rownames(coord_sp)
quant_pc1_taxo1 <- quantile(coord_sp$Dim.1, 0.05)
quant_pc1_taxo2 <- quantile(coord_sp$Dim.1, 0.95)

plt_sp <- ggplot(coord_sp[(coord_sp$Dim.1 <= quant_pc1_taxo1| 
                             coord_sp$Dim.1 >= quant_pc1_taxo2), ],
                 aes(x = 1, y = Dim.1)) + xlim(0.96, 1.005) +
  geom_point(color = 'red') +
  geom_text_repel(aes(label = sp),
                  force        = 6,
                  nudge_x      = -0.1,
                  direction    = "y",
                  hjust        = 1)  +
  ylab("PC1") +
  theme_classic() + theme(axis.text.x =element_blank(),
                          axis.line.x =element_blank(),
                          axis.ticks.x =element_blank(),
                          axis.title.x =element_blank()) 
### temporal 


df_pca_taxo_temporel <- df_pca_taxo %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   n = n_distinct(index_unique)) %>% 
  dplyr::mutate(dim1_up = 1.9*dim1_sd/sqrt(n),
                dim2_up = 1.9*dim2_sd/sqrt(n))


head(coord_sp)

df_pca_taxo_temporel$species <- NA
for(rep_yy in unique(df_pca_taxo_temporel$Year)){
  
  dede1 <- as.matrix(coord_sp[,c( "Dim.1","Dim.2" )])
  dede2 <- as.matrix(df_pca_taxo_temporel[df_pca_taxo_temporel$Year == rep_yy,
                                          c( "dim1_m","dim2_m" )])
  dede3 <- rbind(dede2, dede1)
  rep_dist <- as.matrix(dist(dede3))
  rep_dist <- rep_dist[1, 2:dim(rep_dist)[1]]
  
  df_pca_taxo_temporel[df_pca_taxo_temporel$Year == rep_yy, 'species'] <-   names(rep_dist)[rep_dist == min(rep_dist)] 
}

df_repel <- df_pca_taxo_temporel[1, ]
  
for(n in 2:dim(df_pca_taxo_temporel)[1]){
    if(df_pca_taxo_temporel[n, 'species'] != df_pca_taxo_temporel[n-1, 'species'] ){
      df_repel <- rbind(df_repel, df_pca_taxo_temporel[n, ])
    }
  }
 
df_repel$species_yy <- paste0(df_repel$species, ' (', df_repel$Year, ')')

df_pca_taxo_temporel2 <- merge(df_pca_taxo_temporel, df_repel, all = T)


d1 = ggplot(df_pca_taxo_temporel,aes(x = dim1_m, y = dim2_m, color = Year)) + 
  geom_segment(aes(x = dim1_m - dim1_up, xend = dim1_m + dim1_up,
                   y = dim2_m, yend = dim2_m), alpha = 0.5) +
  geom_segment(aes(x = dim1_m , xend = dim1_m,
                   y = dim2_m - dim2_up, yend = dim2_m + dim2_up), 
               alpha = 0.5) + xlab("PC1") +  ylab("PC2") +
  scale_color_viridis_c() + geom_path() +  geom_point() + 
  theme_classic()  + theme(legend.position = 'none',
                           legend.title = element_blank())

d2 = ggplot(df_pca_taxo_temporel2,
            aes(x = 1, y = Year, color = Year, label = species_yy)) + 
  scale_color_viridis_c() + geom_path() +  geom_point() + 
  geom_label(aes(x = 1, y = Year, color = Year, label = species_yy),
             box.padding = 0.5, max.overlaps = Inf,
             fill = "white") +
  theme_classic() + theme(axis.text.x =element_blank(),
                          axis.line.x =element_blank(),
                          axis.ticks.x =element_blank(),
                          axis.title.x =element_blank(),
                          legend.position = 'none') 

################################# ACP fonctio ##################################

sp_traits1 <- data.frame(df_traits_last) %>% 
  dplyr::filter(!is.na(spawning.type) & !is.na(length.infinity) &
                  spawning.type != "NA" & 
                  genus_sp  %in% unique(temp$genus_sp)) %>% 
  dplyr::mutate(diet = factor(feeding.mode, exclude = TRUE),
                K.growth = growth.coefficient,
                spawning = factor(spawning.type, exclude = TRUE),
                habitat = factor(habitat, exclude = TRUE)) %>% 
  dplyr::arrange(genus_sp)
rownames(sp_traits1) <- sp_traits1$genus_sp

sp_traits <- sp_traits1 %>% 
  dplyr::select(length.maturity, age.maturity,K.growth, 
                tl, spawning, habitat, diet)

sp_traits_cat <- data.frame(trait_name = names(sp_traits),
                            trait_type = c(rep("Q", 4), rep("N", 3)))


# 3. Computing distances between species based on functional traits
sp_dist<- mFD::funct.dist(
  sp_tr         = sp_traits,
  tr_cat        = sp_traits_cat,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = TRUE)
as.matrix(sp_dist)[1:10, 1:5]

# 4. Computing functional spaces & their quality. Computing functional spaces & their qualit
fspaces_quality_fruits <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 10,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

sp_faxes_coord <- fspaces_quality_fruits$"details_fspaces"$"sp_pc_coord"
df_sp_faxes_coord <- data.frame(sp_faxes_coord[,c(1:5)])
df_sp_faxes_coord$genus_sp <- rownames(df_sp_faxes_coord)

rep_traits <- vegan::envfit(sp_faxes_coord[,c(1:5)], sp_traits, perm = 999)
coords_traits_num <- as.data.frame(rep_traits$vectors$arrows)
coords_traits_num$trait <- row.names(coords_traits_num)

coords_traits_fact <- as.data.frame(rep_traits$factors$centroids)
coords_traits_fact$trait <- row.names(coords_traits_fact)
coords_traits_fact$trait <- ifelse(substr(coords_traits_fact$trait, 1, 1) == "s",
                                   paste0(substr(coords_traits_fact$trait, 9, 25), " (spawning)"), 
                                   ifelse(substr(coords_traits_fact$trait, 1, 1) == "h",
                                          paste0(substr(coords_traits_fact$trait, 8, 25), " (habitat)"), 
                                          paste0(substr(coords_traits_fact$trait, 5, 25), " (diet)")))

coords_traits_all <- rbind(coords_traits_num, coords_traits_fact)
coords_traits_all

######################### add dimensions all data

df_for_pcoa <- merge(data.table(temp), 
                     data.table(df_sp_faxes_coord),
                     by = 'genus_sp')

df_for_pcoa2 <- df_for_pcoa %>% 
  dplyr::group_by(index_unique,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pcoa1 = weighted.mean(PC1, log_abundance),
                   pcoa2 = weighted.mean(PC2, log_abundance),
                   pcoa3 = weighted.mean(PC3, log_abundance),
                   pcoa4 = weighted.mean(PC4, log_abundance),
                   pcoa5 = weighted.mean(PC5, log_abundance),
                   depth = mean(depth),
                   chloro_mea = mean(chloro_mea),
                   mlotst_mea = mean(mlotst_mea),
                   oxy_bottom_mea = mean(oxy_bottom_mea),
                   oxy_surf_mea = mean(oxy_surf_mea),
                   temp_surf_mea = mean(temp_surf_mea),
                   temp_bottom_mea = mean(temp_bottom_mea),
                   curr_surf_mea = mean(curr_surf_mea),
                   curr_bottom_mea = mean(curr_bottom_mea),
                   zo_mea = mean(zo_mea))

############ results 

df_for_pcoa3 <- df_for_pcoa2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pcoa1),
                   pcoa2_m  =  mean(pcoa1),
                   pcoa3_m  =  mean(pcoa1),
                   pcoa4_m  =  mean(pcoa1),
                   pcoa5_m  =  mean(pcoa1))

mapfunctio <- ggplot(data = world) +  
  geom_tile(data = df_for_pcoa3,
            aes(x = x_my_spatial_id, y = y_my_spatial_id, fill = pcoa1_m  ),
            size = 1) + scale_fill_viridis_c(name = 'PC1') + 
  geom_sf(color = NA, fill = 'grey80') +  theme_classic() + 
  ggtitle('Functional') +
  coord_sf(xlim = c(min_lon, max_lon), ylim = c( min_lat , max_lat),expand = FALSE) + 
  theme(axis.title = element_blank(),
        legend.position = 'right',
        plot.title = element_text(face = 'bold')) 


quant_pc1 <- quantile(coords_traits_all$PC1)
quant_pc1_fonctio1 <- quantile(coords_traits_all$PC1, 0.1)
quant_pc1_fonctio2 <- quantile(coords_traits_all$PC1, 0.9)

plt_traits <- ggplot(coords_traits_all[(coords_traits_all$PC1 <= quant_pc1_fonctio1 | 
                                          coords_traits_all$PC1 >= quant_pc1_fonctio2), ],
                     aes(x = 1, y = PC1)) + xlim(0.96, 1.005) +
  geom_point(color = 'red') + geom_text_repel(aes(label = trait),
                                              force        = 8,
                                              nudge_x      = -0.1,
                                              direction    = "y",
                                              hjust        = 1)  +
  theme_classic() + theme(axis.text.x =element_blank(),
                          axis.line.x =element_blank(),
                          axis.ticks.x =element_blank(),
                          axis.title.x =element_blank()) 
### temporal 


df_pca_functio_temporel  <- (df_for_pcoa2)  %>% 
  dplyr::group_by(Year) %>% 
  dplyr::summarise(dim1_m = mean(pcoa1),
                   dim1_sd = sd(pcoa1),
                   dim2_m = mean(pcoa2),
                   dim2_sd = sd(pcoa2),
                   n = n_distinct(index_unique)) %>% 
  dplyr::mutate(dim1_up = 1.9*dim1_sd/sqrt(n),
                dim2_up = 1.9*dim2_sd/sqrt(n))



head(coords_traits_all)
df_pca_functio_temporel <- data.frame(df_pca_functio_temporel)
df_pca_functio_temporel$traits <- NA
for(rep_yy in unique(df_pca_functio_temporel$Year)){
  
  dede1 <- as.matrix(coords_traits_all[,c( "PC1","PC2" )])
  dede2 <- as.matrix(df_pca_functio_temporel[df_pca_functio_temporel$Year == rep_yy,
                                          c( "dim1_m","dim2_m" )])
  dede3 <- rbind(dede2, dede1)
  rep_dist <- as.matrix(dist(dede3))
  rep_dist <- rep_dist[1, 2:dim(rep_dist)[1]]
  
  df_pca_functio_temporel[df_pca_functio_temporel$Year == rep_yy, 'traits'] <-   names(rep_dist)[rep_dist == min(rep_dist)] 
}

df_repel_f <- df_pca_functio_temporel[1, ]

for(n in 2:dim(df_pca_functio_temporel)[1]){
  if(df_pca_functio_temporel[n, 'traits'] != df_pca_functio_temporel[n-1, 'traits'] ){
    df_repel_f <- rbind(df_repel_f, df_pca_functio_temporel[n, ])
  }
}

df_repel_f$traits_yy <- paste0(df_repel_f$traits, ' (', df_repel_f$Year, ')')

df_pca_functio_temporel2 <- merge(df_pca_functio_temporel, df_repel_f, all = T)


d1_f = ggplot(df_pca_functio_temporel,aes(x = dim1_m, y = dim2_m, color = Year)) + 
  geom_segment(aes(x = dim1_m - dim1_up, xend = dim1_m + dim1_up,
                   y = dim2_m, yend = dim2_m), alpha = 0.5) +
  geom_segment(aes(x = dim1_m , xend = dim1_m,
                   y = dim2_m - dim2_up, yend = dim2_m + dim2_up), 
               alpha = 0.5) + xlab("PC1") +  ylab("PC2") +
  scale_color_viridis_c() + geom_path() +  geom_point() + 
  theme_classic()  + theme(legend.position = 'none',
                           legend.title = element_blank())

d2_f = ggplot(df_pca_functio_temporel2,
            aes(x = 1, y = Year, color = Year, label = traits_yy)) + 
  scale_color_viridis_c() + geom_path() +  geom_point() + 
  geom_label(aes(x = 1, y = Year, color = Year, label = traits_yy),
             box.padding = 0.5, max.overlaps = Inf,
             fill = "white") +
  theme_classic() + theme(axis.text.x =element_blank(),
                          axis.line.x =element_blank(),
                          axis.ticks.x =element_blank(),
                          axis.title.x =element_blank(),
                          legend.position = 'none') 



##################
taxo_top  = ggarrange(maptaxo, plt_sp, ncol = 2, 
                      widths = c(2, 1.7),
               labels = c(rep_region_quarter, '','',''), hjust = -2, 
               font.label = list(size = 12, face = "plain"))
taxo_bottom  = ggarrange(d1, d2, ncol = 2, 
                      widths = c(2, 1.7))
taxo_all = ggarrange(taxo_top, taxo_bottom , nrow = 2)

func_top  = ggarrange(mapfunctio, plt_traits, ncol = 2, 
                      widths = c(2, 1.7))
func_bottom  = ggarrange(d1_f, d2_f, ncol = 2, 
                         widths = c(2, 1.7))
func_all = ggarrange(func_top, func_bottom , nrow = 2)




a2 = ggarrange(taxo_all, func_all , ncol = 2)

ggsave(file = paste0('figures/taxo_fonctio/figure1_',rep_region_quarter,'.jpg'), 
       plot = a2, width = 5,  height = 2.5, scale = 3)


  }


}
