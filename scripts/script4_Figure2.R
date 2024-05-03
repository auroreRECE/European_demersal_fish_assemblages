library(rgdal)
library(ggplot2)
library(dplyr)
library(data.table)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggrepel)
library(ggpubr)
library(grid)

as.num <- function(x) { x <- as.numeric(as.character(x))}

world <- ne_countries(scale = "medium", returnclass = "sf")

eight_cols <- c("#F8766D","#F99D1E","#FF61CC",
                "#C77CFF","#00BFC4","#74D33A",
                "#7C4728", "#5D6966")


############# load and process taxonomic data ############# 
df_env <- read.csv2(file = 'data/raw/df_env.csv')
df_env <- df_env %>% dplyr::select(-X)

df_PCA <- read.csv2(file = 'data/raw/df_PCA.csv')
df_PCA <- df_PCA %>% dplyr::select(-X)

df_PCA <- merge(df_PCA, df_env)

df_PCA_small1 <- df_PCA %>% 
  dplyr::group_by(ID_unique,  Year, Quarter, depth, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pca2_m  =  mean(Dim.2),
                   pca3_m  =  mean(Dim.3))

df_PCA_small2 <- df_PCA_small1 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pca2_m),
                   pcoa3_m  =  mean(pca3_m))

load(file = "C:/aurore/European_demersal_fish_assemblages/data/raw/res_acp_taxo.Rdata")
res_acp_taxo

############# load and process functional  data ############# 

df_MFA <- read.csv2(file = 'data/raw/df_MFA.csv')
df_MFA <- df_MFA %>% dplyr::select(-X)

df_MFA <- merge(df_MFA, df_env)

df_MFA_small1 <- df_MFA %>% 
  dplyr::group_by(ID_unique,  Year, Quarter, depth, Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pca2_m  =  mean(Dim.2),
                   pca3_m  =  mean(Dim.3))

df_MFA_small2 <- df_MFA_small1 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id, my_spatial_id,
                  Ecoregion) %>% 
  dplyr::summarise(pcoa1_m = mean(pca1_m),
                   pcoa2_m  =  mean(pca2_m),
                   pcoa3_m  =  mean(pca3_m))

load(file = "C:/aurore/European_demersal_fish_assemblages/data/raw/result_mfa.Rdata")
result_mfa

############# create the Figure ############# 
############# spatial part

df_PCA$Ecoregion <- as.factor(df_PCA$Ecoregion)
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_text_taxo <- df_PCA %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2))
df_text_taxo$group <- c('l','l',
                        't','t',
                        'r', 'r',
                        'b','b')
param_lineheit <- 0.6

df_PCA$Ecoregion <- factor(df_PCA$Ecoregion,
                                levels = c("Baltic Sea" , 
                                           "Greater\nNorth Sea", 
                                           "Celtic Seas", 
                                           "Bay of Biscay and\nIberian Coast",
                                           "Western\nMed. Sea",
                                           "Ionian Sea and\nCentral Med. Sea",
                                           "Adriatic Sea"    ,
                                           "Aegean-Levantine\nSea"))

plt_space_taxo = ggplot(df_PCA, aes(x = Dim.1, y = Dim.2, color = Ecoregion ))  +
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


df_MFA$Ecoregion <- as.factor(df_MFA$Ecoregion)
levels(df_MFA$Ecoregion)
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_MFA$Ecoregion <- factor(df_MFA$Ecoregion,
                                   levels = c("Baltic Sea" , 
                                              "Greater\nNorth Sea", 
                                              "Celtic Seas", 
                                              "Bay of Biscay and\nIberian Coast",
                                              "Western\nMed. Sea",
                                              "Ionian Sea and\nCentral Med. Sea",
                                              "Adriatic Sea"    ,
                                              "Aegean-Levantine\nSea"))

df_text_fonctio <- df_MFA %>% 
  dplyr::group_by(Ecoregion) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim2_m = mean(Dim.2))
df_text_fonctio$group <- c('l','b','r','t2',
                           't','b','l','t')


plt_space_functio = ggplot(df_MFA, aes(x = Dim.1, y = Dim.2, color =Ecoregion ))  +
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


############# spaces part

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
  geom_hline(yintercept = 0, color = 'grey70') + 
  geom_vline(xintercept = 0, color = 'grey70') +
  geom_point(color = 'red') + 
  geom_text_repel(aes(label = sp2), force = 5, max.overlaps = Inf,  
                  lineheight = param_lineheit,  fontface = "italic")  +
  ylab("PC2 (12%)") +  xlab("PC1 (21%)") + theme_classic()   +
  theme(panel.background = element_rect(fill = 'white', color = 'black'),
        panel.grid.major = element_line(color = 'grey90')) 


## 

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

plt_trai_dim1_dim2 <- ggplot(coord_traits_sm,aes(x = Dim.1, y = Dim.2)) + 
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


## all together 

fig_top  = ggarrange(plt_space_taxo, plt_sp_dim1_dim2,ncol = 2, nrow = 1, align = 'hv',
                      widths  = c(1.2, 1.2),
                      labels = c("A", 'B'),  font.label = list(size = 12, face = "plain"))

fig_bottom = ggarrange(plt_space_functio, plt_trai_dim1_dim2,  ncol = 2, nrow = 1, align = 'hv',
                        widths  = c(1.2, 1.2),
                        labels = c("C", 'D'),  font.label = list(size = 12, face = "plain"))

fig_all = ggarrange(fig_top, fig_bottom, nrow = 2) 

fig_all2 = annotate_figure(fig_all, 
                            right = textGrob("Taxonomic", rot=-90,
                                             # vjust=-1, 
                                             hjust = 2.1,
                                             gp = gpar(cex = 1.3)))

fig_all3 = annotate_figure(fig_all2, 
                            right = textGrob("Functional", rot=-90,
                                             vjust=2.2,
                                             hjust = -1,
                                             gp = gpar(cex = 1.3)))

ggsave(file = 'figures/figure2.jpg', plot = fig_all3,
       width = 3,  height = 2, units = c("in"), scale = 3)
ggsave(file = 'figures/figure1.jpg', plot = mapdata,
       width = 2.2,  height = 2, scale = 3)
