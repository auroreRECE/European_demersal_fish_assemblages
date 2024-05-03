
library(dplyr)

as.num <- function(x) { x <- as.numeric(as.character(x))}

############# load and process taxonomic data ############# 
df_env <- read.csv2(file = 'data/raw/df_env.csv')
df_env <- df_env %>% dplyr::select(-X)

df_PCA <- read.csv2(file = 'data/raw/df_PCA.csv')
df_PCA <- df_PCA %>% dplyr::select(-X)
df_PCA <- merge(df_PCA, df_env)

############# load and process functional  data ############# 

df_MFA <- read.csv2(file = 'data/raw/df_MFA.csv')
df_MFA <- df_MFA %>% dplyr::select(-X)
df_MFA <- merge(df_MFA, df_env)

################################################################################
################### temporal evolution by bioregion  ###########################
################################################################################

df_MFA$Ecoregion <- as.factor(df_MFA$Ecoregion)
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_MFA$Ecoregion)[levels(df_MFA$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

head(df_MFA)

df_pca_fonctio_sm <- df_MFA %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   n = n_distinct(ID_unique))

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

    temp <- df_pca_fonctio_sm_melt %>% 
      dplyr::filter(Ecoregion == rep_region & dimension == rep_dim)

    my_lm <- lm(mean_value ~ Year, data = temp)
    df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
                             df_pca_fonctio_sm_melt$dimension == rep_dim,
                           "pval"] <- summary(my_lm)$coefficients[2,4]
    df_pca_fonctio_sm_melt[df_pca_fonctio_sm_melt$Ecoregion == rep_region &
                             df_pca_fonctio_sm_melt$dimension == rep_dim,
                           "coeff"] <- summary(my_lm)$coefficients[2,1]

  }
}

write.csv2(df_pca_fonctio_sm_melt,
           file = 'data/intermediate/df_temporal_trends_MFA.csv')

############################################################################################################ 
head(df_PCA)
df_PCA$Ecoregion <- as.factor(df_PCA$Ecoregion)
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Ionian Sea and the Central Mediterranean Sea" ] <- "Ionian Sea and\nCentral Med. Sea" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Bay of Biscay and the Iberian Coast" ] <-"Bay of Biscay and\nIberian Coast" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Western Mediterranean Sea" ] <-"Western\nMed. Sea" 
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Aegean-Levantine Sea"] <- "Aegean-Levantine\nSea"
levels(df_PCA$Ecoregion)[levels(df_PCA$Ecoregion) == "Greater North Sea" ] <- "Greater\nNorth Sea" 

df_pca_taxo_sm <- df_PCA %>% 
  dplyr::group_by(Ecoregion, Year ) %>% 
  dplyr::summarise(dim1_m = mean(Dim.1),
                   dim1_sd = sd(Dim.1),
                   dim2_m = mean(Dim.2),
                   dim2_sd = sd(Dim.2),
                   dim3_m = mean(Dim.3),
                   dim3_sd = sd(Dim.3),
                   n = n_distinct(ID_unique))

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

for(rep_region in unique(df_pca_taxo_sm_melt$Ecoregion)){
  for(rep_dim in unique(df_pca_taxo_sm_melt$dimension)){
    dede <- df_pca_taxo_sm_melt %>% 
      dplyr::filter(Ecoregion == rep_region & dimension == rep_dim)
    
    my_lm <- lm(mean_value ~ Year, data = dede)
    df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$Ecoregion == rep_region &
                          df_pca_taxo_sm_melt$dimension == rep_dim,
                        "pval"] <- summary(my_lm)$coefficients[2,4]
    df_pca_taxo_sm_melt[df_pca_taxo_sm_melt$Ecoregion == rep_region &
                          df_pca_taxo_sm_melt$dimension == rep_dim,
                        "coeff"] <- summary(my_lm)$coefficients[2,1]
    
  }
}

write.csv2(df_pca_taxo_sm_melt,
           file = 'data/intermediate/df_temporal_trends_PCA.csv')
