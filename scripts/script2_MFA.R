

################################################################################
#################################### MFA CWMs ##################################
################################################################################

########################

df_CMW_quanti <- df %>% 
  dplyr::select(ID_unique, genus_sp, log_abundance, 
                log_length.maturity, log_age.maturity, 
                log_growth.coefficient, log_tl)

df_CMW_quanti2 <- df_CMW_quanti %>% 
  dplyr::group_by(ID_unique) %>% 
  dplyr::summarise(length.maturity = weighted.mean(log_length.maturity, log_abundance),
                   age.maturity = weighted.mean(log_age.maturity, log_abundance),
                   growth.coefficient = weighted.mean(log_growth.coefficient, log_abundance),
                   tl = weighted.mean(log_tl, log_abundance))


#### spawning.type
df_CMW_quali <- df_last %>% 
  dplyr::select(ID_unique, genus_sp, log_abundance, 
                spawning.type, habitat, feeding.mode)

df_CMW_quali_cast_spaw <- reshape2::dcast(df_CMW_quali, 
                                          ID_unique + genus_sp + log_abundance ~ spawning.type,
                                          length,
                                          value.var = 'spawning.type')


df_CMW_quali_cast_spaw2 <- df_CMW_quali_cast_spaw %>% 
  dplyr::group_by(ID_unique) %>% 
  dplyr::summarise(bearer = weighted.mean(bearer, log_abundance)*100,
                   guarder = weighted.mean(guarder, log_abundance)*100,
                   non_guarder = weighted.mean(non_guarder, log_abundance)*100)

#### habitat
df_CMW_quali_cast_hab <- reshape2::dcast(df_CMW_quali,
                                         ID_unique + genus_sp + log_abundance ~ habitat,
                                         length,
                                         value.var = "habitat")

df_CMW_quali_cast_hab2 <- df_CMW_quali_cast_hab %>% 
  dplyr::group_by(ID_unique) %>% 
  dplyr::summarise(bathydemersal = weighted.mean(bathydemersal, log_abundance)*100,
                   bathypelagic = weighted.mean(bathypelagic, log_abundance)*100,
                   benthopelagic = weighted.mean(benthopelagic, log_abundance)*100,
                   demersal = weighted.mean(demersal, log_abundance)*100,
                   pelagic = weighted.mean(pelagic, log_abundance)*100,
                   reef_associated = weighted.mean(reef_associated, log_abundance)*100)


#### feeding.mode
df_CMW_quali_cast_feed <- reshape2::dcast(df_CMW_quali,
                                          ID_unique + genus_sp + log_abundance ~ feeding.mode,
                                          length,
                                          value.var = "feeding.mode")

df_CMW_quali_cast_feed2 <- df_CMW_quali_cast_feed %>% 
  dplyr::group_by(ID_unique) %>% 
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

result_mfa <- MFA(df_CMW_all, 
                  group = c(1, 1, 1, 1, 
                            3, 6, 5), 
                  type = c("c", "c", "c", "c", "c", "c", "c"),
                  name.group = c("length.maturity","age.maturity","growth.coefficient","tl", 
                                 "spawning",
                                 "habitat", 
                                 "diet"),
                  graph = FALSE)


df_pca_fonctio <-  data.frame(result_mfa$ind$coord)
df_pca_fonctio$ID_unique <- df_CMW_all$ID_unique

write.csv2(df_pca_fonctio,
           file = '/data/raw/df_MFA.csv')
