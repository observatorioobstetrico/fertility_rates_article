# Loading the necessary libraries
library(knitr)
library(dplyr)
library(factoextra)
library(clusterSim)
library(clValid)
library(htmltools)
library(gridExtra)
library(grid)
library(rstatix)
library(sf)
library(geobr)

# Reading and manipulating the data -------------------------------------------
df_indicadores <- read.csv("databases/data_fertility_rates_article.csv")

## Filtering and calculating the fertility rates for each period
df_indicadores_18_19 <- df_indicadores |>
  filter(ano %in% c(2018, 2019)) |>
  group_by(codmunres) |>
  summarise(
    tx_fecundidade_10_a_14 = round(sum(nvm_10_a_14) / sum(pop_feminina_10_a_14) * 1000, 1),
    tx_fecundidade_15_a_19 = round(sum(nvm_15_a_19) / sum(pop_feminina_15_a_19) * 1000, 1),
    cobertura_ab = round(sum(media_cobertura_ab) / sum(populacao_total) * 100, 1),
    porc_exclusivas_sus = round(sum(pop_fem_10_49_com_plano_saude) / sum(pop_feminina_10_a_49) * 100, 1),
    idhm = mean(idhm)
  ) |>
  ungroup()

df_indicadores_20_21 <- df_indicadores |>
  filter(ano %in% c(2020, 2021)) |>
  group_by(codmunres) |>
  summarise(
    tx_fecundidade_10_a_14 = round(sum(nvm_10_a_14) / sum(pop_feminina_10_a_14) * 1000, 1),
    tx_fecundidade_15_a_19 = round(sum(nvm_15_a_19) / sum(pop_feminina_15_a_19) * 1000, 1),
    cobertura_ab = round(sum(media_cobertura_ab) / sum(populacao_total) * 100, 1),
    porc_exclusivas_sus = round(sum(pop_fem_10_49_com_plano_saude) / sum(pop_feminina_10_a_49) * 100, 1),
    idhm = mean(idhm)
  ) |>
  ungroup()

## Removing auxiliary objects
rm(df_indicadores)

# For the cluster analysis ----------------------------------------------------
## Standardizing the data for each of the four cluster analysis
dados1_18_19 <- df_indicadores_18_19 |> select(tx_fecundidade_10_a_14) |> scale()
dados2_18_19 <- df_indicadores_18_19 |> select(tx_fecundidade_15_a_19) |> scale()

dados1_20_21 <- df_indicadores_20_21 |> select(tx_fecundidade_10_a_14) |> scale()
dados2_20_21 <- df_indicadores_20_21 |> select(tx_fecundidade_15_a_19) |> scale()

## Calculating the distance matrices
dados1_18_19_dist <- dist(dados1_18_19, method = "euclidean")
dados2_18_19_dist <- dist(dados2_18_19, method = "euclidean")

dados1_20_21_dist <- dist(dados1_20_21, method = "euclidean")
dados2_20_21_dist <- dist(dados2_20_21, method = "euclidean")


## For the 10 to 14 age group -------------------------------------------------
### K-means
#### Ploting the k-means knee-plots
fviz_nbclust(dados1_18_19, kmeans, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

fviz_nbclust(dados1_20_21, kmeans, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

#### Adjusting the k-means method with the chosen numbers of clusters 
set.seed(1504)
d1_kmeans3_18_19 <- kmeans(dados1_18_19, 3)

set.seed(1504)
d1_kmeans4_18_19 <- kmeans(dados1_18_19, 4)

set.seed(1504)
d1_kmeans3_20_21 <- kmeans(dados1_20_21, 3)

set.seed(1504)
d1_kmeans4_20_21 <- kmeans(dados1_20_21, 4)


### K-medians
#### Ploting the k-medians knee-plots
fviz_nbclust(dados1_18_19, cluster::pam, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

fviz_nbclust(dados1_20_21, cluster::pam, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

#### Adjusting the k-medians method with the chosen numbers of clusters 
set.seed(1504)
d1_pam3_18_19 <- cluster::pam(dados1_18_19, k = 3)

set.seed(1504)
d1_pam4_18_19 <- cluster::pam(dados1_18_19, k = 4)

set.seed(1504)
d1_pam3_20_21 <- cluster::pam(dados1_20_21, k = 3)

set.seed(1504)
d1_pam4_20_21 <- cluster::pam(dados1_20_21, k = 4)


### Ward's method
#### Ploting the drodrograms
set.seed(1504)
d1_ward_18_19 <- hclust(dados1_18_19_dist, method = "ward.D2")
plot(as.dendrogram(d1_ward_18_19), ylab = "Altura") 
rect.hclust(d1_ward_18_19, k = 2, border = 2:5) 
rect.hclust(d1_ward_18_19, k = 3, border = 2:5) # Candidates: K = 2 and K = 3.

set.seed(1504)
d1_ward_20_21 <- hclust(dados1_20_21_dist, method = "ward.D2")
plot(as.dendrogram(d1_ward_20_21), ylab = "Altura") 
rect.hclust(d1_ward_20_21, k = 2, border = 2:5) 
rect.hclust(d1_ward_20_21, k = 4, border = 2:5) # Candidates: K = 2 or K = 4.

#### Adjusting the Ward's method with the chosen numbers of clusters 
d1_ward2_18_19_class <- cutree(d1_ward_18_19, k = 2)
d1_ward3_18_19_class <- cutree(d1_ward_18_19, k = 3)

d1_ward2_20_21_class <- cutree(d1_ward_20_21, k = 2)
d1_ward4_20_21_class <- cutree(d1_ward_20_21, k = 4)


### Single linkage
#### Ploting the dendrograms
set.seed(1504)
d1_single_18_19 <- hclust(dados1_18_19_dist, method = "single")
plot(as.dendrogram(d1_single_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d1_single_20_21 <- hclust(dados1_20_21_dist, method = "single")
plot(as.dendrogram(d1_single_20_21), ylab = "Altura") # No candidates. Bad fit.


### Complete linkage
#### Ploting the dendrograms
set.seed(1504)
d1_complete_18_19 <- hclust(dados1_18_19_dist, method = "complete")
plot(as.dendrogram(d1_complete_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d1_complete_20_21 <- hclust(dados1_20_21_dist, method = "complete")
plot(as.dendrogram(d1_complete_20_21), ylab = "Altura") # No candidates. Bad fit.


### Average linkage
#### Ploting the dendrograms
set.seed(1504)
d1_average_18_19 <- hclust(dados1_18_19_dist, method = "average")
plot(as.dendrogram(d1_average_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d1_average_20_21 <- hclust(dados1_20_21_dist, method = "average")
plot(as.dendrogram(d1_average_20_21), ylab = "Altura") # No candidates. Bad fit. 


### Choosing the best clustering method
#### For the 2018-2019 period 
d1_kmeans_index_18_19 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("kmeans", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados1_18_19, get(paste0("d1_kmeans", i, "_18_19"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados1_18_19, get(paste0("d1_kmeans", i, "_18_19"))$cluster, d = dados1_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados1_18_19_dist, get(paste0("d1_kmeans", i, "_18_19"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados1_18_19_dist, get(paste0("d1_kmeans", i, "_18_19"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados1_18_19, get(paste0("d1_kmeans", i, "_18_19"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados1_18_19, get(paste0("d1_kmeans", i, "_18_19"))$cluster, d = dados1_18_19_dist, centrotypes = "medoids")))
)

d1_pam_index_18_19 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("pam", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados1_18_19, get(paste0("d1_pam", i, "_18_19"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados1_18_19, get(paste0("d1_pam", i, "_18_19"))$cluster, d = dados1_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados1_18_19_dist, get(paste0("d1_pam", i, "_18_19"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados1_18_19_dist, get(paste0("d1_pam", i, "_18_19"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados1_18_19, get(paste0("d1_pam", i, "_18_19"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados1_18_19, get(paste0("d1_pam", i, "_18_19"))$cluster, d = dados1_18_19_dist, centrotypes = "medoids")))
)

d1_ward_index_18_19 <- data.frame(
  metodo = unlist(lapply(2:3, function(i) paste0("ward", i))),
  db_index_cent = unlist(lapply(2:3, function(i) index.DB(dados1_18_19, get(paste0("d1_ward", i, "_18_19", "_class")))$DB)),
  db_index_med = unlist(lapply(2:3, function(i) index.DB(dados1_18_19, get(paste0("d1_ward", i, "_18_19", "_class")), d = dados1_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(2:3, function(i) dunn(distance = dados1_18_19_dist, get(paste0("d1_ward", i, "_18_19", "_class"))))),
  silh_index = unlist(lapply(2:3, function(i) index.S(dados1_18_19_dist, get(paste0("d1_ward", i, "_18_19", "_class"))))),
  ch_index_cent = unlist(lapply(2:3, function(i) index.G1(dados1_18_19, get(paste0("d1_ward", i, "_18_19", "_class"))))),
  ch_index_med = unlist(lapply(2:3, function(i) index.G1(dados1_18_19, get(paste0("d1_ward", i, "_18_19", "_class")), d = dados1_18_19_dist, centrotypes = "medoids")))
)

d1_avaliacao_18_19 <- rbind(
  d1_pam_index_18_19,
  d1_kmeans_index_18_19, 
  d1_ward_index_18_19
)
d1_avaliacao_18_19
#### Comments:
##### db_index_cent: the lower the better (kmeans4, followed by pam4)
##### db_index_med: the lower the better (kmeans3, followed by pam3 and kmeans4)
##### dunn_index: the higher the better (kmeans4, followed by ward3 and kmeans3)
##### silh_index: the higher the better (ward2, followed by pam4, kmeans4 and kmeans3)
##### ch_index_cent: the higher the better (kmeans4, followed by kmeans3)
##### ch_index_med: the higher the better (kmeans4, followed by kmeans3)

#### Choice: despite the better perfomance of k-means with K = 4, one of its clusters contained only a few municipalities (127).
d1_kmeans4_18_19$size
d1_kmeans3_18_19$size
#### Thus, the chosen method is the k-means with K = 3.

#### For the 2020-2021 period
d1_kmeans_index_20_21 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("kmeans", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados1_20_21, get(paste0("d1_kmeans", i, "_20_21"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados1_20_21, get(paste0("d1_kmeans", i, "_20_21"))$cluster, d = dados1_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados1_20_21_dist, get(paste0("d1_kmeans", i, "_20_21"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados1_20_21_dist, get(paste0("d1_kmeans", i, "_20_21"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados1_20_21, get(paste0("d1_kmeans", i, "_20_21"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados1_20_21, get(paste0("d1_kmeans", i, "_20_21"))$cluster, d = dados1_20_21_dist, centrotypes = "medoids")))
)

d1_pam_index_20_21 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("pam", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados1_20_21, get(paste0("d1_pam", i, "_20_21"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados1_20_21, get(paste0("d1_pam", i, "_20_21"))$cluster, d = dados1_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados1_20_21_dist, get(paste0("d1_pam", i, "_20_21"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados1_20_21_dist, get(paste0("d1_pam", i, "_20_21"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados1_20_21, get(paste0("d1_pam", i, "_20_21"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados1_20_21, get(paste0("d1_pam", i, "_20_21"))$cluster, d = dados1_20_21_dist, centrotypes = "medoids")))
)

d1_ward_index_20_21 <- data.frame(
  metodo = unlist(lapply(c(2, 4), function(i) paste0("ward", i))),
  db_index_cent = unlist(lapply(c(2, 4), function(i) index.DB(dados1_20_21, get(paste0("d1_ward", i, "_20_21", "_class")))$DB)),
  db_index_med = unlist(lapply(c(2, 4), function(i) index.DB(dados1_20_21, get(paste0("d1_ward", i, "_20_21", "_class")), d = dados1_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(c(2, 4), function(i) dunn(distance = dados1_20_21_dist, get(paste0("d1_ward", i, "_20_21", "_class"))))),
  silh_index = unlist(lapply(c(2, 4), function(i) index.S(dados1_20_21_dist, get(paste0("d1_ward", i, "_20_21", "_class"))))),
  ch_index_cent = unlist(lapply(c(2, 4), function(i) index.G1(dados1_20_21, get(paste0("d1_ward", i, "_20_21", "_class"))))),
  ch_index_med = unlist(lapply(c(2, 4), function(i) index.G1(dados1_20_21, get(paste0("d1_ward", i, "_20_21", "_class")), d = dados1_20_21_dist, centrotypes = "medoids")))
)

d1_avaliacao_20_21 <- rbind(
  d1_pam_index_20_21,
  d1_kmeans_index_20_21, 
  d1_ward_index_20_21
)
d1_avaliacao_20_21
#### Comments:
##### db_index_cent: the lower the better (ward4, followed by kmeans4 and kmeans3)
##### db_index_med: the lower the better (ward4, followed by kmeans4 and kmeans3)
##### dunn_index: the higher the better (ward4, followed by kmeans4 and kmeans3)
##### silh_index: the higher the better (ward2, followed by kmeans4 and kmeans3)
##### ch_index_cent: the higher the better (kmeans4, followed by ward4 and kmeans3)
##### ch_index_med: the higher the better (kmeans4, followed by ward4 and kmeans3)

#### Choice: again, despite the better perfomance of k-means with K = 4 and Ward's method with K = 4, one of their clusters contained 
#### only a few municipalities (90 and 70, respectively).
d1_kmeans4_20_21$size
table(d1_ward4_20_21_class)
d1_kmeans3_18_19$size
#### Thus, the chosen method is the k-means with K = 3.

### Adding the columns with the clusters informations to the original data.frames
#### For the 2018-2019 period
df_indicadores_18_19$cluster_10_a_14 <- d1_kmeans3_18_19$cluster

ordem_taxa1_18_19 <- df_indicadores_18_19 |>
  group_by(cluster_10_a_14) |>
  summarise(tx_fecundidade_10_a_14 = mean(tx_fecundidade_10_a_14)) |>
  ungroup() |>
  arrange(tx_fecundidade_10_a_14) |>
  pull(cluster_10_a_14)

df_indicadores_18_19$cluster_10_a_14 <- factor(case_when(
  df_indicadores_18_19$cluster_10_a_14 == ordem_taxa1_18_19[1] ~ "1: lowest fertility rates",
  df_indicadores_18_19$cluster_10_a_14 == ordem_taxa1_18_19[2] ~ "2: intermediary fertility rates",
  df_indicadores_18_19$cluster_10_a_14 == ordem_taxa1_18_19[3] ~ "3: highest fertility rates"
))

#### For the 2020-2021 period
df_indicadores_20_21$cluster_10_a_14 <- d1_kmeans3_20_21$cluster

ordem_taxa1_20_21 <- df_indicadores_20_21 |>
  group_by(cluster_10_a_14) |>
  summarise(tx_fecundidade_10_a_14 = mean(tx_fecundidade_10_a_14)) |>
  ungroup() |>
  arrange(tx_fecundidade_10_a_14) |>
  pull(cluster_10_a_14)

df_indicadores_20_21$cluster_10_a_14 <- factor(case_when(
  df_indicadores_20_21$cluster_10_a_14 == ordem_taxa1_20_21[1] ~ "1: lowest fertility rates",
  df_indicadores_20_21$cluster_10_a_14 == ordem_taxa1_20_21[2] ~ "2: intermediary fertility rates",
  df_indicadores_20_21$cluster_10_a_14 == ordem_taxa1_20_21[3] ~ "3: highest fertility rates"
))


## For the 15 to 19 age group -------------------------------------------------
### K-means
#### Ploting the k-means knee-plots
set.seed(1504)
fviz_nbclust(dados2_18_19, kmeans, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

set.seed(1504)
fviz_nbclust(dados2_20_21, kmeans, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

#### Adjusting the k-means method with the chosen numbers of clusters 
set.seed(1504)
d2_kmeans3_18_19 <- kmeans(dados2_18_19, 3)

set.seed(1504)
d2_kmeans4_18_19 <- kmeans(dados2_18_19, 4)

set.seed(1504)
d2_kmeans3_20_21 <- kmeans(dados2_20_21, 3)

set.seed(1504)
d2_kmeans4_20_21 <- kmeans(dados2_20_21, 4)


### K-medians
#### Ploting the k-medians knee-plots
set.seed(1504)
fviz_nbclust(dados2_18_19, cluster::pam, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3.

set.seed(1504)
fviz_nbclust(dados2_20_21, cluster::pam, method = "wss") +
  labs(
    x = "Número de clusters", y = "Variância total intragrupo", title = ""
  ) # Candidates: K = 3 and K = 4.

#### Adjusting the k-medians method with the chosen numbers of clusters 
set.seed(1504)
d2_pam3_18_19 <- cluster::pam(dados2_18_19, k = 3)

set.seed(1504)
d2_pam3_20_21 <- cluster::pam(dados2_20_21, k = 3)

set.seed(1504)
d2_pam4_20_21 <- cluster::pam(dados2_20_21, k = 4)


### Ward's method
#### Ploting the drodrograms
set.seed(1504)
d2_ward_18_19 <- hclust(dados2_18_19_dist, method = "ward.D2")
plot(as.dendrogram(d2_ward_18_19), ylab = "Altura") 
rect.hclust(d2_ward_18_19, k = 3, border = 2:5) # Candidates: K = 3.

set.seed(1504)
d2_ward_20_21 <- hclust(dados2_20_21_dist, method = "ward.D2")
plot(as.dendrogram(d2_ward_20_21), ylab = "Altura") 
rect.hclust(d2_ward_20_21, k = 3, border = 2:5) # Candidates: K = 3.

#### Adjusting the Ward's method with the chosen numbers of clusters 
d2_ward3_18_19_class <- cutree(d2_ward_18_19, k = 3)
d2_ward3_20_21_class <- cutree(d2_ward_20_21, k = 3)


### Single linkage
#### Ploting the drodrograms
set.seed(1504)
d2_single_18_19 <- hclust(dados2_18_19_dist, method = "single")
plot(as.dendrogram(d2_single_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d2_single_20_21 <- hclust(dados2_20_21_dist, method = "single")
plot(as.dendrogram(d2_single_20_21), ylab = "Altura") # No candidates. Bad fit.


### Complete linkage
#### Ploting the drodrograms
set.seed(1504)
d2_complete_18_19 <- hclust(dados2_18_19_dist, method = "complete")
plot(as.dendrogram(d2_complete_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d2_complete_20_21 <- hclust(dados2_20_21_dist, method = "complete")
plot(as.dendrogram(d2_complete_20_21), ylab = "Altura") # No candidates. Bad fit.


### Average linkage
#### Ploting the drodrograms
set.seed(1504)
d2_average_18_19 <- hclust(dados2_18_19_dist, method = "average")
plot(as.dendrogram(d2_average_18_19), ylab = "Altura") # No candidates. Bad fit.

set.seed(1504)
d2_average_20_21 <- hclust(dados2_20_21_dist, method = "average")
plot(as.dendrogram(d2_average_20_21), ylab = "Altura") # No candidates. Bad fit.


### Chosing the best clustering method
#### For the 2018-2019 period
d2_kmeans_index_18_19 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("kmeans", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados2_18_19, get(paste0("d2_kmeans", i, "_18_19"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados2_18_19, get(paste0("d2_kmeans", i, "_18_19"))$cluster, d = dados2_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados2_18_19_dist, get(paste0("d2_kmeans", i, "_18_19"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados2_18_19_dist, get(paste0("d2_kmeans", i, "_18_19"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados2_18_19, get(paste0("d2_kmeans", i, "_18_19"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados2_18_19, get(paste0("d2_kmeans", i, "_18_19"))$cluster, d = dados2_18_19_dist, centrotypes = "medoids")))
)

d2_pam_index_18_19 <- data.frame(
  metodo = unlist(lapply(3, function(i) paste0("pam", i))),
  db_index_cent = unlist(lapply(3, function(i) index.DB(dados2_18_19, get(paste0("d2_pam", i, "_18_19"))$cluster)$DB)),
  db_index_med = unlist(lapply(3, function(i) index.DB(dados2_18_19, get(paste0("d2_pam", i, "_18_19"))$cluster, d = dados2_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3, function(i) dunn(distance = dados2_18_19_dist, get(paste0("d2_pam", i, "_18_19"))$cluster))),
  silh_index = unlist(lapply(3, function(i) index.S(dados2_18_19_dist, get(paste0("d2_pam", i, "_18_19"))$cluster))),
  ch_index_cent = unlist(lapply(3, function(i) index.G1(dados2_18_19, get(paste0("d2_pam", i, "_18_19"))$cluster))),
  ch_index_med = unlist(lapply(3, function(i) index.G1(dados2_18_19, get(paste0("d2_pam", i, "_18_19"))$cluster, d = dados2_18_19_dist, centrotypes = "medoids")))
)

d2_ward_index_18_19 <- data.frame(
  metodo = unlist(lapply(3, function(i) paste0("ward", i))),
  db_index_cent = unlist(lapply(3, function(i) index.DB(dados2_18_19, get(paste0("d2_ward", i, "_18_19", "_class")))$DB)),
  db_index_med = unlist(lapply(3, function(i) index.DB(dados2_18_19, get(paste0("d2_ward", i, "_18_19", "_class")), d = dados2_18_19_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3, function(i) dunn(distance = dados2_18_19_dist, get(paste0("d2_ward", i, "_18_19", "_class"))))),
  silh_index = unlist(lapply(3, function(i) index.S(dados2_18_19_dist, get(paste0("d2_ward", i, "_18_19", "_class"))))),
  ch_index_cent = unlist(lapply(3, function(i) index.G1(dados2_18_19, get(paste0("d2_ward", i, "_18_19", "_class"))))),
  ch_index_med = unlist(lapply(3, function(i) index.G1(dados2_18_19, get(paste0("d2_ward", i, "_18_19", "_class")), d = dados2_18_19_dist, centrotypes = "medoids")))
)

d2_avaliacao_18_19 <- rbind(
  d2_pam_index_18_19,
  d2_kmeans_index_18_19, 
  d2_ward_index_18_19
)
d2_avaliacao_18_19
#### Comments:
##### db_index_cent: the lower the better (kmeans3, followed by ward3 and kmeans4)
##### db_index_med: the lower the better (kmeans3, followed by ward3)
##### dunn_index: the higher the better (kmeans3, followed by ward3)
##### silh_index: the higher the better (kmeans3, followed by kmeans4)
##### ch_index_cent: the higher the better (kmeans4, followed by kmeans3)
##### ch_index_med: the higher the better (kmeans4, followed by kmeans3)

#### Choice: k-means with K = 3.

#### For the 2020-2021 period
d2_kmeans_index_20_21 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("kmeans", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados2_20_21, get(paste0("d2_kmeans", i, "_20_21"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados2_20_21, get(paste0("d2_kmeans", i, "_20_21"))$cluster, d = dados2_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados2_20_21_dist, get(paste0("d2_kmeans", i, "_20_21"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados2_20_21_dist, get(paste0("d2_kmeans", i, "_20_21"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados2_20_21, get(paste0("d2_kmeans", i, "_20_21"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados2_20_21, get(paste0("d2_kmeans", i, "_20_21"))$cluster, d = dados2_20_21_dist, centrotypes = "medoids")))
)

d2_pam_index_20_21 <- data.frame(
  metodo = unlist(lapply(3:4, function(i) paste0("pam", i))),
  db_index_cent = unlist(lapply(3:4, function(i) index.DB(dados2_20_21, get(paste0("d2_pam", i, "_20_21"))$cluster)$DB)),
  db_index_med = unlist(lapply(3:4, function(i) index.DB(dados2_20_21, get(paste0("d2_pam", i, "_20_21"))$cluster, d = dados2_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3:4, function(i) dunn(distance = dados2_20_21_dist, get(paste0("d2_pam", i, "_20_21"))$cluster))),
  silh_index = unlist(lapply(3:4, function(i) index.S(dados2_20_21_dist, get(paste0("d2_pam", i, "_20_21"))$cluster))),
  ch_index_cent = unlist(lapply(3:4, function(i) index.G1(dados2_20_21, get(paste0("d2_pam", i, "_20_21"))$cluster))),
  ch_index_med = unlist(lapply(3:4, function(i) index.G1(dados2_20_21, get(paste0("d2_pam", i, "_20_21"))$cluster, d = dados2_20_21_dist, centrotypes = "medoids")))
)

d2_ward_index_20_21 <- data.frame(
  metodo = unlist(lapply(3, function(i) paste0("ward", i))),
  db_index_cent = unlist(lapply(3, function(i) index.DB(dados2_20_21, get(paste0("d2_ward", i, "_20_21", "_class")))$DB)),
  db_index_med = unlist(lapply(3, function(i) index.DB(dados2_20_21, get(paste0("d2_ward", i, "_20_21", "_class")), d = dados2_20_21_dist, centrotypes = "medoids")$DB)),
  dunn_index = unlist(lapply(3, function(i) dunn(distance = dados2_20_21_dist, get(paste0("d2_ward", i, "_20_21", "_class"))))),
  silh_index = unlist(lapply(3, function(i) index.S(dados2_20_21_dist, get(paste0("d2_ward", i, "_20_21", "_class"))))),
  ch_index_cent = unlist(lapply(3, function(i) index.G1(dados2_20_21, get(paste0("d2_ward", i, "_20_21", "_class"))))),
  ch_index_med = unlist(lapply(3, function(i) index.G1(dados2_20_21, get(paste0("d2_ward", i, "_20_21", "_class")), d = dados2_20_21_dist, centrotypes = "medoids")))
)

d2_avaliacao_20_21 <- rbind(
  d2_pam_index_20_21,
  d2_kmeans_index_20_21, 
  d2_ward_index_20_21
)
d2_avaliacao_20_21
#### Comments:
##### db_index_cent: the lower the better (kmeans4, followed by ward3 and kmeans3)
##### db_index_med: the lower the better (kmeans4, followed by kmeans3 and ward3)
##### dunn_index: the higher the better (kmeans4, followed by kmeans3)
##### silh_index: the higher the better (kmeans4, followed by kmeans3)
##### ch_index_cent: the higher the better (kmeans4, followed by kmeans3)
##### ch_index_med: the higher the better (kmeans4, followed by kmeans3)

#### Choice: even though, this time around, the k-means with K = 4 didn't have problems with too little municipalities in a cluster,
#### I'm chosing the k-means with K = 3 to make it "comparable" with the other chosen clustering methods.
d2_kmeans4_20_21$size
#### Thus, the chosen method is the k-means with K = 3.

### Adding the columns with the clusters informations to the original data.frames
#### For the 2018-2019 period
df_indicadores_18_19$cluster_15_a_19 <- d2_kmeans3_18_19$cluster

ordem_taxa2_18_19 <- df_indicadores_18_19 |>
  group_by(cluster_15_a_19) |>
  summarise(tx_fecundidade_15_a_19 = mean(tx_fecundidade_15_a_19)) |>
  ungroup() |>
  arrange(tx_fecundidade_15_a_19) |>
  pull(cluster_15_a_19)

df_indicadores_18_19$cluster_15_a_19 <- factor(case_when(
  df_indicadores_18_19$cluster_15_a_19 == ordem_taxa2_18_19[1] ~ "1: lowest fertility rates",
  df_indicadores_18_19$cluster_15_a_19 == ordem_taxa2_18_19[2] ~ "2: intermediary fertility rates",
  df_indicadores_18_19$cluster_15_a_19 == ordem_taxa2_18_19[3] ~ "3: highest fertility rates"
))

#### For the 2020-2021 period
df_indicadores_20_21$cluster_15_a_19 <- d2_kmeans3_20_21$cluster

ordem_taxa2_20_21 <- df_indicadores_20_21 |>
  group_by(cluster_15_a_19) |>
  summarise(tx_fecundidade_15_a_19 = mean(tx_fecundidade_15_a_19)) |>
  ungroup() |>
  arrange(tx_fecundidade_15_a_19) |>
  pull(cluster_15_a_19)

df_indicadores_20_21$cluster_15_a_19 <- factor(case_when(
  df_indicadores_20_21$cluster_15_a_19 == ordem_taxa2_20_21[1] ~ "1: lowest fertility rates",
  df_indicadores_20_21$cluster_15_a_19 == ordem_taxa2_20_21[2] ~ "2: intermediary fertility rates",
  df_indicadores_20_21$cluster_15_a_19 == ordem_taxa2_20_21[3] ~ "3: highest fertility rates"
))


## Visualizing the clusters in a map ---------------------------------------
### Downloading the geometry data
df_muni_sf <- read_municipality(year = 2019, showProgress = FALSE) |>
  mutate(codmunres = as.numeric(substr(code_muni, 1, 6)))

df_ufs_sf <- read_state(year = 2019, showProgress = FALSE)

### Joining the geometry and the clustering data for each period
df_dados_mapa_18_19 <- left_join(df_indicadores_18_19, df_muni_sf) |>
  mutate(periodo = "2018-2019") 

df_dados_mapa_20_21 <- left_join(df_indicadores_20_21, df_muni_sf) |>
  mutate(periodo = "2020-2021") 

df_dados_mapa_completo <- full_join(df_dados_mapa_18_19, df_dados_mapa_20_21) |>
  st_as_sf()

### Plotting the maps
#### For the 10 to 14 age group
plot_clusters_10_a_14 <- ggplot() +
  geom_sf(data = df_dados_mapa_completo, aes(fill = cluster_10_a_14), color = NA) +
  facet_wrap(vars(periodo)) +
  scale_fill_viridis_d(name = "Groups", end = 0.8, alpha = 0.6, direction = -1) +
  geom_sf(data = df_ufs_sf, fill = NA, linewidth = 0.08, color = "black") +
  theme_bw() +
  theme(legend.position = "bottom") 

##### Exporting the plot
ggsave(
  "figures/Fig2.tiff", 
  plot_clusters_10_a_14,
  width = 7.5, height = 6, units = "in", 
  dpi = 600, compression = "lzw"
)

#### For the 15 to 19 age group
plot_clusters_15_a_19 <- ggplot() +
  geom_sf(data = df_dados_mapa_completo, aes(fill = cluster_15_a_19), color = NA) +
  facet_wrap(vars(periodo)) +
  scale_fill_viridis_d(name = "Groups", end = 0.8, alpha = 0.6, direction = -1) +
  geom_sf(data = df_ufs_sf, fill = NA, linewidth = 0.08, color = "black") +
  theme_bw() +
  theme(legend.position = "bottom")

##### Exporting the plot
ggsave(
  "figures/Fig3.tiff", 
  plot_clusters_15_a_19,
  width = 7.5, height = 6, units = "in", 
  dpi = 600, compression = "lzw"
)


## Comparing the groups -------------------------------------------------------
### Creating a function that creates a table with the summary measures 
cria_tabelas <- function(variaveis, df, var_grupos, label_tabela, faixa_etaria, periodo) {
  grupos <- levels(df[[var_grupos]])
  
  for (variavel in variaveis) {
    tabela <- kable(
      data.frame(
        grupo = grupos,
        n = unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> nrow())),
        minimo = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> min(na.rm = TRUE))), 2),
        primeiro_qt = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> quantile(0.25, na.rm = TRUE))), 2),
        media = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> mean(na.rm = TRUE))), 2),
        mediana = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> median(na.rm = TRUE))), 2),
        dp = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> sd(na.rm = TRUE))), 2),
        terceiro_qt = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> quantile(0.75, na.rm = TRUE))), 2),
        maximo = round(unlist(lapply(grupos, function(grupo) df[df[[var_grupos]] == grupo, ] |> pull(variavel) |> max(na.rm = TRUE))), 2)
      ),
      align = "cccccccc",
      col.names = c("Grupo (cluster)", "n", "Mín.", "1º Quartil", "Média", "Mediana", "D.P.", "3º Quartil", "Máx."),
      caption = HTML(ifelse(variavel == variaveis[1], paste0(glue::glue("Medidas resumo das variáveis de interesse para os grupos de municípios agrupados pela taxa de fecundidade de mulheres {faixa_etaria} (por mil) no período de {periodo}. <br><br>"), variavel), variavel)),
      label = ifelse(variavel == variaveis[1], paste0("tabela", label_tabela), NA)
    )
    
    print(tabela)
  }
}

### Creating a function that creates a table with the multiple comparisions
cria_tabelas_testes <- function(variaveis, df, var_grupos, label_tabela, faixa_etaria, periodo) {
  for (variavel in variaveis) {
    
    formula <- as.formula(glue::glue("{variavel} ~ {var_grupos}"))
    resultados_dunn <- dunn_test(data = df, formula = formula, p.adjust.method = "bonferroni")
    
    # Calcula a medida de efeito d de Cohen
    efeito_d <- cohens_d(data = df, formula = formula)
    
    monta_linhas <- function(num_grupo) {
      paste(
        -1 * round(resultados_dunn$statistic[which(resultados_dunn$group1 == unique(resultados_dunn$group1)[num_grupo])], 3),
        "(Z) <br>",
        ifelse(
          round(resultados_dunn$p.adj[which(resultados_dunn$group1 == unique(resultados_dunn$group1)[num_grupo])], 3) < 0.05,
          ifelse(
            round(resultados_dunn$p.adj[which(resultados_dunn$group1 == unique(resultados_dunn$group1)[num_grupo])], 3) == 0,
            "< 0.001*",
            paste0(round(resultados_dunn$p.adj[which(resultados_dunn$group1 == unique(resultados_dunn$group1)[num_grupo])], 3), "*")
          ),
          round(resultados_dunn$p.adj[which(resultados_dunn$group1 == unique(resultados_dunn$group1)[num_grupo])], 3)
        ),
        "(valor-p) <br>",
        round(efeito_d$effsize[efeito_d$group1 == unique(efeito_d$group1)[num_grupo]], 3),
        "(D de Cohen)"
      )
    }
    
    df_teste <- data.frame(
      grupo1 = monta_linhas(1),
      grupo2 = c(rep("", 1), monta_linhas(2)),
      row.names = unique(resultados_dunn$group2)
    )
    
    colnames(df_teste) <- unique(resultados_dunn$group1)
    
    tabela <- kable(
      df_teste,
      align = c("cc"),
      caption = HTML(ifelse(variavel == variaveis[1], paste0(glue::glue("Resultados dos testes de Dunn (com correção de Bonferroni) para as comparações múltiplas entre os pares de grupos de municípios agrupados pela taxa de fecundidade de mulheres {faixa_etaria} (por mil) no período de {periodo}. <br><br>"), variavel), variavel)),
      label = ifelse(variavel == variaveis[1], paste0("tabela-comparacoes-", label_tabela), NA)
    )
    
    print(tabela)
  }
}
nomes_variaveis1 <- c("idhm", "cobertura_ab", "tx_fecundidade_10_a_14")
nomes_variaveis2 <- c("idhm", "cobertura_ab", "tx_fecundidade_15_a_19")

### Creating the summary tables
#### For the 10 to 14 age group
cria_tabelas(
  df = df_indicadores_18_19, var_grupos = "cluster_10_a_14", label_tabela = "1",
  faixa_etaria = "de 10 a 14 anos", variaveis = nomes_variaveis1, periodo = "2018 a 2019"
)

cria_tabelas(
  df = df_indicadores_20_21, var_grupos = "cluster_10_a_14", label_tabela = "2", 
  faixa_etaria = "de 10 a 14 anos", variaveis = nomes_variaveis1, periodo = "2020 a 2021"
)

#### For the 15 to 19 age group
cria_tabelas(
  df = df_indicadores_18_19, var_grupos = "cluster_15_a_19", label_tabela = "3", 
  faixa_etaria = "de 15 a 19 anos", variaveis = nomes_variaveis2, periodo = "2018 a 2019"
)

cria_tabelas(
  df = df_indicadores_20_21, var_grupos = "cluster_15_a_19", label_tabela = "4", 
  faixa_etaria = "de 15 a 19 anos", variaveis = nomes_variaveis2, periodo = "2020 a 2021"
)


### Creating the multiple comparisions tables
#### For the 10 to 14 age group
cria_tabelas_testes(
  df = df_indicadores_18_19, var_grupos = "cluster_10_a_14", label_tabela = "1",
  faixa_etaria = "de 10 a 14 anos", variaveis = nomes_variaveis1, periodo = "2018 a 2019"
)

cria_tabelas_testes(
  df = df_indicadores_20_21, var_grupos = "cluster_10_a_14", label_tabela = "2", 
  faixa_etaria = "de 10 a 14 anos", variaveis = nomes_variaveis1, periodo = "2020 a 2021"
)

#### For the 15 to 19 age group
cria_tabelas_testes(
  df = df_indicadores_18_19, var_grupos = "cluster_15_a_19", label_tabela = "3", 
  faixa_etaria = "de 15 a 19 anos", variaveis = nomes_variaveis2, periodo = "2018 a 2019"
)

cria_tabelas_testes(
  df = df_indicadores_20_21, var_grupos = "cluster_15_a_19", label_tabela = "4", 
  faixa_etaria = "de 15 a 19 anos", variaveis = nomes_variaveis2, periodo = "2020 a 2021"
)
