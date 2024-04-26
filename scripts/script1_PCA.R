
## this script gives the code tu run the PCA. The 'df' dataframe has three columns : the index ID of each trawl, the column with the species names, and the column with the log-transformed abundances

trawl_sp_matrice <- reshape2::acast(data = df,
                                    ID_unique ~ genus_sp,
                                    value.var = 'log_abundance',
                                    fun.aggregate = mean)
trawl_sp_matrice[is.na(trawl_sp_matrice)] <- 0
trawl_sp_matrice_hell <- decostand(trawl_sp_matrice, method = 'hellinger')


res_acp_taxo <- PCA(trawl_sp_matrice_hell, scale.unit = FALSE, ncp = 50,
                    graph = FALSE)

############ results 

df_pca_taxo <-  data.frame(res_acp_taxo$ind$coord)
df_pca_taxo$ID_unique <- rownames(df_pca_taxo)

write.csv2(df_pca_taxo, file = '/data/raw/df_PCA.csv')
