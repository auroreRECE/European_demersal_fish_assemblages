
library(dplyr)
library(nlme)
library(MuMIn)

as.num <- function(x) { x <- as.numeric(as.character(x))}

############# load and process taxonomic data ############# 
df_env <- read.csv2(file = 'data/raw/df_env.csv')
df_env <- df_env %>% dplyr::select(-X)

df_PCA <- read.csv2(file = 'data/raw/df_PCA.csv')
df_PCA <- df_PCA %>% dplyr::select(-X)
df_PCA <- merge(df_PCA, df_env)

load(file = "C:/aurore/European_demersal_fish_assemblages/data/raw/res_acp_taxo.Rdata")
# res_acp_taxo

############# load and process functional  data ############# 

df_MFA <- read.csv2(file = 'data/raw/df_MFA.csv')
df_MFA <- df_MFA %>% dplyr::select(-X)
df_MFA <- merge(df_MFA, df_env)

load(file = "C:/aurore/European_demersal_fish_assemblages/data/raw/result_mfa.Rdata")
# result_mfa

################################################################################
################################### GLS effet env ######################
################################################################################

head(df_MFA)
df_MFA$log_chloro <- log(df_MFA$chloro_mea + 1)
df_MFA$log_mlotst <- log(df_MFA$mlotst_mea + 1)

df_MFA$log_curr_bottom <- log(df_MFA$curr_bottom_mea) + 10
df_MFA$log_curr_surf <- log(df_MFA$curr_surf_mea) + 10


df_MFA_small <- df_MFA %>% 
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
                   sal_surf = mean(sal_surf_mea, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE),
                   
                   n = n_distinct(ID_unique))


df_MFA_small$chloro_scale <- scale(df_MFA_small$chloro)
df_MFA_small$mlotst_scale <- scale(df_MFA_small$mlotst)
df_MFA_small$oxy_bottom_scale <- scale(df_MFA_small$oxy_bottom)
df_MFA_small$oxy_surf_scale <- scale(df_MFA_small$oxy_surf)
df_MFA_small$temp_bottom_scale <- scale(df_MFA_small$temp_bottom)
df_MFA_small$temp_surf_scale <- scale(df_MFA_small$temp_surf)
df_MFA_small$curr_bottom_scale <- scale(df_MFA_small$curr_bottom)
df_MFA_small$curr_surf_scale <- scale(df_MFA_small$curr_surf)
df_MFA_small$sal_surf_scale <- scale(df_MFA_small$sal_surf)

df_MFA_small$chloro_std_scale <- scale(df_MFA_small$chloro_std)
df_MFA_small$mlotst_std_scale <- scale(df_MFA_small$mlotst_std)
df_MFA_small$oxy_bottom_std_scale <- scale(df_MFA_small$oxy_bottom_std)
df_MFA_small$oxy_surf_std_scale <- scale(df_MFA_small$oxy_surf_std)
df_MFA_small$temp_bottom_std_scale <- scale(df_MFA_small$temp_bottom_std)
df_MFA_small$temp_surf_std_scale <- scale(df_MFA_small$temp_surf_std)
df_MFA_small$curr_surf_std_scale <- scale(df_MFA_small$curr_surf_std)
df_MFA_small$curr_bottom_std_scale <- scale(df_MFA_small$curr_bottom_std)
df_MFA_small$sal_surf_std_scale <- scale(df_MFA_small$sal_surf_std)


df_MFA_small <- data.frame(df_MFA_small)
liste_ecoregion <- unique(df_MFA_small$Ecoregion )

for(rep_ecoR in liste_ecoregion){
  print(rep_ecoR)

  for(rep_dim in c('dim1_m','dim2_m')){
    print(rep_dim)
    
    temp <- df_MFA_small %>% 
      dplyr::filter(Ecoregion == rep_ecoR) %>% 
      dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                      !is.na(temp_bottom) &  !is.na(temp_surf) &  !is.na(curr_bottom))

      gls1 = nlme::gls(formula(paste0(rep_dim, " ~  temp_bottom_scale + temp_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + 
                      
                      temp_bottom_std_scale + temp_surf_std_scale + mlotst_std_scale +
                      oxy_bottom_std_scale + chloro_std_scale ")),
                       data = temp,
                       correlation = corGaus(form =~ Year)
      )
 
    options(na.action = "na.fail")
    dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
    options(na.action = "na.omit")
    aa = subset(dredge.glm1.BIC.p, delta < 2)
    
    mod <- get.models(aa, 1)[[1]]

    table_coeff <- as.data.frame(summary(mod)$tTable)
    table_coeff$variable_env <- rownames(table_coeff)
    table_coeff <- reshape2::melt(table_coeff)
    table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
    table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )
    
    df_results <- table_coeff
    df_results$Ecoregion <- rep_ecoR
    df_results$dim <- rep_dim
    
    vect_obc <- temp[,rep_dim]
    temp$pred <- predict(object = mod, newdata = temp)
    SSE1 <- sum(c(vect_obc - temp$pred) ^ 2, na.rm = TRUE)
    SSTot1 <-  sum(c(vect_obc - mean(vect_obc)) ^ 2)
    my_rsquared1 <- round((1 - SSE1/SSTot1) * 100, 3)

    df_results$R2 <- my_rsquared1
    
    if(rep_ecoR == liste_ecoregion[1] & rep_dim == 'dim1_m'){
      df_results_last <- df_results
    } else {
      df_results_last <- rbind(df_results_last, df_results)
    }
    
  }
  
}
names(df_results_last)[names(df_results_last) == "p-value"] <- 'p_val'
names(df_results_last)[names(df_results_last) == "Std.Error"] <- 'std_error'

df_results_temporel_functio <- df_results_last %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(Ecoregion, dim)


############ spatial 

df_pca_fonctio <-  data.frame(result_mfa$ind$coord)
df_pca_fonctio$ID_unique <- df_CMW_all2$ID_unique

df_pca_fonctio <- merge(df_complementary, df_pca_fonctio)
df_pca_fonctio$log_chloro <- log(df_pca_fonctio$chloro_mea + 1)
df_pca_fonctio$log_mlotst <- log(df_pca_fonctio$mlotst_mea + 1)

df_pca_fonctio$log_curr_bottom <- log(df_pca_fonctio$curr_bottom_mea + 1)
df_pca_fonctio$log_curr_surf <- log(df_pca_fonctio$curr_surf_mea + 1)


df_pca_fonctio2 <- df_pca_fonctio %>% 
  dplyr::group_by(ID_unique,  Year, Quarter, depth,Ecoregion,
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
                   sal_surf = mean(sal_surf_mea, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE),
                   
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
                   curr_surf = mean(curr_surf, na.rm = TRUE),
                   sal_surf = mean(sal_surf, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE),)

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
df_pca_fonctio3$sal_surf_scale <- scale(df_pca_fonctio3$sal_surf)

df_pca_fonctio3$chloro_std_scale <- scale(df_pca_fonctio3$chloro_std)
df_pca_fonctio3$mlotst_std_scale <- scale(df_pca_fonctio3$mlotst_std)
df_pca_fonctio3$oxy_bottom_std_scale <- scale(df_pca_fonctio3$oxy_bottom_std)
df_pca_fonctio3$oxy_surf_std_scale <- scale(df_pca_fonctio3$oxy_surf_std)
df_pca_fonctio3$temp_bottom_std_scale <- scale(df_pca_fonctio3$temp_bottom_std)
df_pca_fonctio3$temp_surf_std_scale <- scale(df_pca_fonctio3$temp_surf_std)
df_pca_fonctio3$curr_surf_std_scale <- scale(df_pca_fonctio3$curr_surf_std)
df_pca_fonctio3$curr_bottom_std_scale <- scale(df_pca_fonctio3$curr_bottom_std)
df_pca_fonctio3$sal_surf_std_scale <- scale(df_pca_fonctio3$sal_surf_std)


df_pca_fonctio3 <- data.frame(df_pca_fonctio3)
df_pca_fonctio3$index <-  row_number(df_pca_fonctio3$x_my_spatial_id)

############ 

for(rep_dim in c('dim1_m','dim2_m')){
  print(rep_dim)

  temp <- df_pca_fonctio3 %>% 
    dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                    !is.na(temp_bottom) &  !is.na(temp_surf) &  
                    !is.na(depth) &  !is.na(curr_bottom))

  gls1 = nlme::gls(formula(paste0(rep_dim, " ~ depth_scale + depth_span_scale + 
                          temp_bottom_scale + temp_surf_scale +sal_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale  +
                      
                      temp_bottom_std_scale +   oxy_bottom_std_scale + chloro_std_scale ")),
                   data = temp,
                   correlation = corGaus(form =~ x_my_spatial_id + y_my_spatial_id))
  
  
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
  
  vect_obc <- dede[,rep_dim]
  dede$pred <- predict(object = mod, newdata = dede)
  SSE1 <- sum(c(vect_obc - dede$pred) ^ 2, na.rm = TRUE)
  SSTot1 <-  sum(c(vect_obc - mean(vect_obc)) ^ 2)
  my_rsquared1 <- round((1 - SSE1/SSTot1) * 100, 3)

  df_results$R2 <- my_rsquared1
  
  
  if( rep_dim == 'dim1_m'){
    df_results_last_spatial <- df_results } else {
    df_results_last_spatial <- rbind(df_results_last_spatial, df_results) }
  
}

names(df_results_last_spatial)[names(df_results_last_spatial) == "p-value"] <- 'p_val'
names(df_results_last_spatial)[names(df_results_last_spatial) == "Std.Error"] <- 'std_error'

df_results_spatial_functio <- df_results_last_spatial %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(dim)

df_results_spatial_functio$Ecoregion <- 'Spatial distribution'

df_all_func <- rbind(df_results_temporel_functio, 
                     df_results_spatial_functio )

write.csv2(file = 'data/intermediate/df_all_func.csv', 
           df_all_func, 
           row.numbers = FALSE)

################################################################################
################################### GLS effet env taxo #########################
################################################################################

head(df_PCA)
df_PCA$log_chloro <- log(df_PCA$chloro_mea + 1)
df_PCA$log_mlotst <- log(df_PCA$mlotst_mea + 1)

df_PCA$log_curr_bottom <- log(df_PCA$curr_bottom_mea + 1) 
df_PCA$log_curr_surf <- log(df_PCA$curr_surf_mea + 1)

df_PCA_sm <- df_PCA %>% 
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
                   sal_surf = mean(sal_surf_mea, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE),
                   
                   n = n_distinct(index_old))
df_pca_taxo_sm$chloro_scale <- scale(df_pca_taxo_sm$chloro)
df_pca_taxo_sm$mlotst_scale <- scale(df_pca_taxo_sm$mlotst)
df_pca_taxo_sm$oxy_bottom_scale <- scale(df_pca_taxo_sm$oxy_bottom)
df_pca_taxo_sm$oxy_surf_scale <- scale(df_pca_taxo_sm$oxy_surf)
df_pca_taxo_sm$temp_bottom_scale <- scale(df_pca_taxo_sm$temp_bottom)
df_pca_taxo_sm$temp_surf_scale <- scale(df_pca_taxo_sm$temp_surf)
df_pca_taxo_sm$curr_bottom_scale <- scale(df_pca_taxo_sm$curr_bottom)
df_pca_taxo_sm$curr_surf_scale <- scale(df_pca_taxo_sm$curr_surf)
df_pca_taxo_sm$sal_surf_scale <- scale(df_pca_taxo_sm$sal_surf)

df_pca_taxo_sm$chloro_std_scale <- scale(df_pca_taxo_sm$chloro_std)
df_pca_taxo_sm$mlotst_std_scale <- scale(df_pca_taxo_sm$mlotst_std)
df_pca_taxo_sm$oxy_bottom_std_scale <- scale(df_pca_taxo_sm$oxy_bottom_std)
df_pca_taxo_sm$oxy_surf_std_scale <- scale(df_pca_taxo_sm$oxy_surf_std)
df_pca_taxo_sm$temp_bottom_std_scale <- scale(df_pca_taxo_sm$temp_bottom_std)
df_pca_taxo_sm$temp_surf_std_scale <- scale(df_pca_taxo_sm$temp_surf_std)
df_pca_taxo_sm$curr_surf_std_scale <- scale(df_pca_taxo_sm$curr_surf_std)
df_pca_taxo_sm$curr_bottom_std_scale <- scale(df_pca_taxo_sm$curr_bottom_std)
df_pca_taxo_sm$sal_surf_std_scale <- scale(df_pca_taxo_sm$sal_surf_std)

df_pca_taxo_sm <- data.frame(df_pca_taxo_sm)
liste_ecoregion <- unique(df_pca_taxo_sm$Ecoregion )

for(rep_ecoR in liste_ecoregion){
  print(rep_ecoR)

  for(rep_dim in c('dim1_m','dim2_m')){
    print(rep_dim)
    
    temp <- df_pca_taxo_sm %>% 
      dplyr::filter(Ecoregion == rep_ecoR) %>% 
      dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                      !is.na(temp_bottom) &  !is.na(temp_surf) &  !is.na(curr_bottom))

    
      gls1 = nlme::gls(formula(paste0(rep_dim, " ~  temp_bottom_scale + temp_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + 
                      
                      temp_bottom_std_scale + temp_surf_std_scale + mlotst_std_scale +
                      oxy_bottom_std_scale + chloro_std_scale ")),
                       data = temp,
                       correlation = corGaus(form =~ Year)
      )

    
    
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
    
    vect_obc <- temp[,rep_dim]
    temp$pred <- predict(object = mod, newdata = temp)
    SSE1 <- sum(c(vect_obc - temp$pred) ^ 2, na.rm = TRUE)
    SSTot1 <-  sum(c(vect_obc - mean(vect_obc)) ^ 2)
    my_rsquared1 <- round((1 - SSE1/SSTot1) * 100, 3)

    df_results <- table_coeff
    df_results$Ecoregion <- rep_ecoR
    df_results$dim <- rep_dim
    df_results$R2 <- my_rsquared1
    

    
    if(rep_ecoR == liste_ecoregion[1] & rep_dim == 'dim1_m'){
      df_results_last_t <- df_results
    } else {
      df_results_last_t <- rbind(df_results_last_t, df_results)
    }
    
  }
  
}
names(df_results_last_t)[names(df_results_last_t) == "p-value"] <- 'p_val'
names(df_results_last_t)[names(df_results_last_t) == "Std.Error"] <- 'std_error'

df_results_temporel_taxo <- df_results_last_t %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(Ecoregion, dim)

############ spatial 
head(df_PCA)
df_PCA$log_chloro <- log(df_PCA$chloro_mea + 1)
df_PCA$log_mlotst <- log(df_PCA$mlotst_mea + 1)

df_PCA$log_curr_bottom <- log(df_PCA$curr_bottom_mea + 1) 
df_PCA$log_curr_surf <- log(df_PCA$curr_surf_mea + 1)

df_pca_taxo2 <- df_PCA %>% 
  dplyr::group_by(ID_unique,  Year, Quarter, depth,Ecoregion,
                  x_my_spatial_id, y_my_spatial_id, my_spatial_id) %>% 
  dplyr::summarise(pca1_m = mean(Dim.1),
                   pco2_m  =  mean(Dim.2),

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
                   sal_surf = mean(sal_surf_mea, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE),
                   
                   n = n_distinct(ID_unique))

df_pca_taxo3 <- df_pca_taxo2 %>% 
  dplyr::group_by(x_my_spatial_id, y_my_spatial_id) %>% 
  dplyr::summarise(dim1_m = mean(pca1_m),
                   dim2_m  =  mean(pco2_m),
                   
                   depth = mean(depth, na.rm = TRUE),
                   depth_span = mean(depth_span, na.rm = TRUE),
                   chloro = mean(chloro, na.rm = TRUE),
                   mlotst = mean(mlotst, na.rm = TRUE),
                   oxy_bottom = mean(oxy_bottom, na.rm = TRUE),
                   oxy_surf = mean(oxy_surf, na.rm = TRUE),
                   temp_bottom = mean(temp_bottom, na.rm = TRUE),
                   temp_surf = mean(temp_surf, na.rm = TRUE),
                   curr_bottom = mean(curr_bottom, na.rm = TRUE),
                   curr_surf = mean(curr_surf, na.rm = TRUE), 
                   sal_surf = mean(sal_surf, na.rm = TRUE),
                   
                   chloro_std = mean(chloro_std, na.rm = TRUE),
                   mlotst_std = mean(mlotst_std, na.rm = TRUE),
                   oxy_bottom_std = mean(oxy_bottom_std, na.rm = TRUE),
                   oxy_surf_std = mean(oxy_surf_std, na.rm = TRUE),
                   temp_bottom_std = mean(temp_bottom_std, na.rm = TRUE),
                   temp_surf_std = mean(temp_surf_std, na.rm = TRUE),
                   curr_surf_std = mean(curr_surf_std, na.rm = TRUE),
                   curr_bottom_std = mean(curr_bottom_std, na.rm = TRUE),
                   sal_surf_std = mean(sal_surf_std, na.rm = TRUE))

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

df_pca_taxo3$sal_surf_scale <- scale(df_pca_taxo3$sal_surf)

df_pca_taxo3$chloro_std_scale <- scale(df_pca_taxo3$chloro_std)
df_pca_taxo3$mlotst_std_scale <- scale(df_pca_taxo3$mlotst_std)
df_pca_taxo3$oxy_bottom_std_scale <- scale(df_pca_taxo3$oxy_bottom_std)
df_pca_taxo3$oxy_surf_std_scale <- scale(df_pca_taxo3$oxy_surf_std)
df_pca_taxo3$temp_bottom_std_scale <- scale(df_pca_taxo3$temp_bottom_std)
df_pca_taxo3$temp_surf_std_scale <- scale(df_pca_taxo3$temp_surf_std)
df_pca_taxo3$curr_surf_std_scale <- scale(df_pca_taxo3$curr_surf_std)
df_pca_taxo3$curr_bottom_std_scale <- scale(df_pca_taxo3$curr_bottom_std)
df_pca_taxo3$sal_surf_std_scale <- scale(df_pca_taxo3$sal_surf_std)


df_pca_taxo3 <- data.frame(df_pca_taxo3)
df_pca_taxo3$index <-  row_number( df_pca_taxo3$x_my_spatial_id)


for(rep_dim in c('dim1_m','dim2_m')){
  print(rep_dim)
  
  temp <- df_pca_taxo3 %>% 
    dplyr::filter(!is.na(chloro) & !is.na(mlotst) & !is.na(oxy_bottom) &
                    !is.na(temp_bottom) &  !is.na(temp_surf) &  
                    !is.na(depth) &  !is.na(curr_bottom))

  
  gls1 = nlme::gls(formula(paste0(rep_dim, " ~ depth_scale + depth_span_scale + 
                          temp_bottom_scale + temp_surf_scale + sal_surf_scale +
                       mlotst_scale + oxy_bottom_scale  + chloro_scale + 
                      
                      temp_bottom_std_scale +   oxy_bottom_std_scale + chloro_std_scale ")),
                   data = temp,
                   correlation = corGaus(form =~ x_my_spatial_id + y_my_spatial_id))
  
  options(na.action = "na.fail")
  dredge.glm1.BIC.p <- dredge(gls1,  rank ="AICc")
  options(na.action = "na.omit")
  aa = subset(dredge.glm1.BIC.p, delta < 2)
  
  mod <- get.models(aa, 1)[[1]]

  vect_obc <- dede[,rep_dim]
  dede$pred <- predict(object = mod, newdata = dede)
  SSE1 <- sum(c(vect_obc - dede$pred) ^ 2, na.rm = TRUE)
  SSTot1 <-  sum(c(vect_obc - mean(vect_obc)) ^ 2)
  my_rsquared1 <- round((1 - SSE1/SSTot1) * 100, 3)

  table_coeff <- as.data.frame(summary(mod)$tTable)
  table_coeff$variable_env <- rownames(table_coeff)
  table_coeff <- reshape2::melt(table_coeff)
  table_coeff <- table_coeff %>% dplyr::filter(variable %in% c("Value", "p-value", "Std.Error"))
  table_coeff <- reshape2::dcast(table_coeff, variable_env  ~ variable )
  
  df_results <- table_coeff
  df_results$dim <- rep_dim
  df_results$R2 <- my_rsquared1

  if( rep_dim == 'dim1_m'){
    df_results_last_spatial_t <- df_results
  } else {
    df_results_last_spatial_t <- rbind(df_results_last_spatial_t, df_results)
  }
  
}

names(df_results_last_spatial_t)[names(df_results_last_spatial_t) == "p-value"] <- 'p_val'
names(df_results_last_spatial_t)[names(df_results_last_spatial_t) == "Std.Error"] <- 'std_error'

df_results_spatial_taxo <- df_results_last_spatial_t %>% 
  dplyr::filter(variable_env != '(Intercept)' & p_val <= 0.05) %>% 
  dplyr::arrange(dim)

df_results_spatial_taxo$Ecoregion <- 'Spatial distribution'

df_all_taxo <- rbind(df_results_temporel_taxo, df_results_spatial_taxo )


write.csv2(file = 'data/intermediate/df_all_taxo.csv', 
           df_all_taxo, 
           row.numbers = FALSE)
