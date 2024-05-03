library(rgdal)
library(ggplot2)
library(dplyr)
library(data.table)
library(rnaturalearthdata)
library(rnaturalearth)

as.num <- function(x) { x <- as.numeric(as.character(x))}

world <- ne_countries(scale = "medium", returnclass = "sf")


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


############# load biological data ############# 

df_env <- read.csv2(file = 'data/raw/df_env.csv')

head(df_env)

############# carte 
df_ecoregion3 <- df_ecoregion3 %>%
  dplyr::filter(Ecoregion %in% unique(df_env$Ecoregion))

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

### please note that on this version, the longitude and latitude are the center of the ICES rectangles
# whereas in the published version, the map is done with the exact trawl coordinates
mapdata <- ggplot(data = world) +
  geom_sf(color = 'white', fill = 'grey50', size = 0.01) +  theme_classic() +
  coord_sf(xlim = c(-21, 31), ylim = c( 33.2 , 65),expand = TRUE) +
  geom_point(data = df_env, aes(x = x_my_spatial_id, y = y_my_spatial_id, color = Survey),
             show.legend = FALSE, size = 1) +
  geom_polygon(data = df_ecoregion3, fill = NA, size = 0.5, color = '#006FC4', alpha = 0.5,
               aes(x = long, y = lat, group =  interaction(Ecoregion, piece))) +
  scale_x_continuous(sec.axis = dup_axis()) +
  geom_label(data = df_ecoregion3_s, aes(x = lon_m, y = lat_m, label = Ecoregion)) +
  theme(axis.title = element_blank(),
        panel.background = element_rect(color = 'black', fill = NA, size = 1),
        axis.line.x.top = element_line(color = 'black',size = 1),
        axis.line.y.right = element_line(color = 'black',size = 1))

ggsave(file = 'figures/figure1.jpg', plot = mapdata,
       width = 2.2,  height = 2, scale = 3)

