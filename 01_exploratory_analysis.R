# Loading the necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(trend)

# Reading the data ------------------------------------------------------------
## Reading a table with auxiliary data for the municipalities
df_aux_municios <- read.csv("databases/dados_aux_municipios.csv") |>
  dplyr::select(codmunres, municipio, uf, regiao, idhm) |>
  mutate(regiao = factor(regiao, levels = c("Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul")))

## Reading a table with the necessary data for calculating the fertility rates
df_fecundidade_aux <- read.csv("databases/dados_fecundidade_outras_faixas.csv") 

## Joining both data.frames
df_fecundidade <- left_join(df_fecundidade_aux, df_aux_municios)

## Aggregating the data for the hole country
df_fecundidade_br <- df_fecundidade |>
  group_by(ano) |>
  summarise(
    nvm_10_a_14 = sum(nvm_10_a_14),
    pop_feminina_10_a_14 = sum(pop_feminina_10_a_14),
    tx_fecundidade_10_a_14 = round(nvm_10_a_14 / pop_feminina_10_a_14 * 1000, 1),
    
    nvm_15_a_19 = sum(nvm_15_a_19),
    pop_feminina_15_a_19 = sum(pop_feminina_15_a_19),
    tx_fecundidade_15_a_19 = round(nvm_15_a_19 / pop_feminina_15_a_19 * 1000, 1)
  ) |>
  ungroup() 

## Aggregating the data for each state
df_fecundidade_ufs <- df_fecundidade |>
  group_by(ano, uf) |>
  summarise(
    nvm_10_a_14 = sum(nvm_10_a_14),
    pop_feminina_10_a_14 = sum(pop_feminina_10_a_14),
    tx_fecundidade_10_a_14 = round(nvm_10_a_14 / pop_feminina_10_a_14 * 1000, 1),
    
    nvm_15_a_19 = sum(nvm_15_a_19),
    pop_feminina_15_a_19 = sum(pop_feminina_15_a_19),
    tx_fecundidade_15_a_19 = round(nvm_15_a_19 / pop_feminina_15_a_19 * 1000, 1)
  ) |>
  ungroup() |>
  mutate(
    uf = factor(uf, levels = df_aux_municios |> dplyr::select(uf, regiao) |> arrange(regiao, uf) |> pull(uf) |> unique())
  ) |>
  arrange(uf)


# Ploting Brazil's time series for both indicators ----------------------------
## Transforming the data.frame to the long format
df_fecundidade_br_unico <- df_fecundidade_br |>
  pivot_longer(
    cols = starts_with("tx"),
    names_to = "faixa_etaria",
    values_to = "tx_fecundidade"
  ) |>
  mutate(
    faixa_etaria = case_when(
      faixa_etaria == "tx_fecundidade_10_a_14" ~ "10 to 14 years",
      faixa_etaria == "tx_fecundidade_15_a_19" ~ "15 to 19 years"
    )
  )

## Ploting the time series
plot_br_time_series <- ggplot(data = df_fecundidade_br_unico, mapping = aes(x = ano, y = tx_fecundidade, color = faixa_etaria)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = faixa_etaria)) +
  scale_x_continuous(breaks = unique(df_fecundidade_br_unico$ano), guide = guide_axis(angle = 45)) +
  theme_bw(base_size = 15) +
  labs(
    color = "Age groups",
    shape = "Age groups",
    x = "Year",
    y = "Fertility rate"
  ) +
  geom_text(label = df_fecundidade_br_unico$tx_fecundidade, nudge_y = 2, show.legend = FALSE) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("Salmon", "DodgerBlue"))

## Exporting the plot
ggsave(
  "figures/Fig1.tiff", plot_br_time_series,
  width = 7.5, height = 6, units = "in", 
  dpi = 600, compression = "lzw"
)


# Mann-Kendall tests for trend analysis ---------------------------------------
## Joining the country and the states data
df_fecundidade_tendencia <- full_join(
  df_fecundidade_br |> mutate(local = "Brasil", .before = "ano"),
  df_fecundidade_ufs |> rename(local = uf)
)

## Selecting the desired variables
variaveis <- c("tx_fecundidade_10_a_14", "tx_fecundidade_15_a_19")

## Creating a data.frame with the Mann-Kendall tests results
results_table_completa <- data.frame()

for (i in 1:length(unique(df_fecundidade_tendencia$local))) {
  localidade <- unique(df_fecundidade_tendencia$local)[i]
  
  mann_kendall_results <- lapply(df_fecundidade_tendencia |> filter(local == localidade) |> dplyr::select(all_of(variaveis)), function(x) {
    mann_kendall_test <- mk.test(x)
    return(c(local = localidade, p_value = mann_kendall_test$p.value, z_value = mann_kendall_test$statistic))
  })
  
  p_value <- round(as.numeric(unlist(lapply(mann_kendall_results, `[[`, "p_value"))), 5)
  
  results_table <- data.frame(
    local = as.vector(unlist(lapply(mann_kendall_results, `[[`, "local"))),
    Variavel = names(mann_kendall_results), 
    valor_2012 = lapply(variaveis, function(variavel) df_fecundidade_tendencia |> filter(local == localidade, ano == 2012) |> pull(variavel)) |> as.numeric(),
    valor_2019 = lapply(variaveis, function(variavel) df_fecundidade_tendencia |> filter(local == localidade, ano == 2019) |> pull(variavel)) |> as.numeric(),
    valor_2020 = lapply(variaveis, function(variavel) df_fecundidade_tendencia |> filter(local == localidade, ano == 2020) |> pull(variavel)) |> as.numeric(),
    valor_2021 = lapply(variaveis, function(variavel) df_fecundidade_tendencia |> filter(local == localidade, ano == 2021) |> pull(variavel)) |> as.numeric(),
    Mann_Kendall_z = round(as.numeric(unlist(lapply(mann_kendall_results, `[[`, "z_value.z"))), 3),
    Mann_Kendall_p = p_value
  )
  
  results_table_completa <- bind_rows(results_table_completa, results_table)
}

## Organizing the resulting table
results_table_completa_organizada <- results_table_completa |>
  filter(startsWith(Variavel, "tx")) |>
  pivot_wider(
    names_from = Variavel,
    values_from = c(valor_2012, valor_2019, valor_2020, valor_2021, Mann_Kendall_z, Mann_Kendall_p)
  ) |>
  select_at(vars("local", ends_with("10_a_14"), ends_with("15_a_19")))

## Exportando the final table
write.csv2(results_table_completa_organizada, "databases/tabela_mann_kendall_completa.csv", row.names = FALSE)
