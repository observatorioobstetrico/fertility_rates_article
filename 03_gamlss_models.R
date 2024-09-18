library(dplyr)
library(tidyr)
library(janitor)
library(gamlss)
library(ggplot2)

# Reading and manipulating the data -------------------------------------------
## Reading a table with auxiliary data for the municipalities
df_aux_municios <- read.csv("databases/dados_aux_municipios.csv") |>
  dplyr::select(codmunres, municipio, uf, regiao, idhm) |>
  mutate(regiao = factor(regiao, levels = c("Norte", "Nordeste", "Centro-Oeste", "Sudeste", "Sul")))

## Reading a table with the necessary data for calculating the fertility rate
df_fecundidade_muni <- read.csv("databases/dados_fecundidade_menores_20.csv") |>
  dplyr::select(codmunres, ano, nvm_menor_que_20, pop_feminina_10_a_19) |>
  filter(ano >= 2018 & ano <= 2021) |>
  mutate(
    codmunres = as.character(codmunres),
    tx_fecundidade_menor_20 = round(nvm_menor_que_20 / pop_feminina_10_a_19 * 1000, 1)
  ) 

## Reading a table with the HDI-M and the necessary data for calculating the primary care coverage 
df_primary_care <- read.csv("databases/dados_indicadores_auxiliares.csv") |>
  dplyr::select(codmunres, ano, media_cobertura_esf, populacao_total, pop_fem_10_49_com_plano_saude, populacao_feminina_10_a_49) |>
  filter(ano >= 2018 & ano <= 2021) |>
  mutate(
    codmunres = as.character(codmunres),
    cobertura_ab = round(media_cobertura_esf / populacao_total * 100, 1),
    .keep = "unused"
  ) |>
  left_join(df_aux_municios |> mutate(codmunres = as.character(codmunres)) |> dplyr::select(!uf))

## Joining the fertility rate and the primary care coverage data
df_indicadores <- left_join(df_fecundidade_muni, df_primary_care) |>
  mutate(
    pandemia = factor(ifelse(ano < 2020, "before", "after"), levels = c("before", "after")),
    ano = as.factor(ano)
  ) |>
  drop_na()  # There are 5 municipalities without the HDI-M information

## Removing auxiliary data
rm(df_aux_municios, df_primary_care, df_fecundidade_muni)
gc()


# Adjusting the models --------------------------------------------------------
## For the ZAGA distribution --------------------------------------------------
### Adjusting the full model 
fit1_zaga <- gamlss(
  sqrt(tx_fecundidade_menor_20) ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  data = df_indicadores,
  family = ZAGA(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
summary(fit1_zaga)
plot(fit1_zaga)
wp(fit1_zaga, ylim.all = 0.6)

### Utilizing the two strategies described in the book "Flexible Regression and Smoothing: Using GAMLSS in R", page 397, to get the best models for each parameter
#### Adjusting a model with only the intercept for mu and considering sigma and nu as constants 
fit2_zaga <- gamlss(
  sqrt(tx_fecundidade_menor_20) ~ 1,
  data = df_indicadores,
  family = ZAGA,
  control = gamlss.control(n.cyc = 500)
)

#### Strategy 1:
#### Starting from the simplest model and utilizing the forward stepwise method to find the model with the lowest AIC for each parameter
#### Obs.: in this strategy, different variables can be selected for each parameter's model
fit3_zaga <- stepGAICAll.A(
  fit2_zaga,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
summary(fit3_zaga)

#### Strategy 2:
#### This strategy forces all os the models to select the same variables subsets 
fit4_zaga <- stepGAICAll.B(
  fit2_zaga,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
summary(fit4_zaga)

#### Using the LR test to verify if the simpler models are better than the full model 
LR.test(fit3_zaga, fit1_zaga)  # The model selected by the first strategy ISN'T better than the full model 
LR.test(fit4_zaga, fit1_zaga)  # The model selected by the second strategy ISN'T better than the full model 
GAIC(fit1_zaga, fit3_zaga, fit4_zaga)

#### Residual analysis of the selected model
plot(fit1_zaga)
wp(fit1_zaga, ylim.all = 0.6)


## For the ST4 distribution ---------------------------------------------------
### Adjusting the full model
fit1_st4 <- gamlss(
  sqrt(tx_fecundidade_menor_20) ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  data = df_indicadores,
  family = ST4(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  tau.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
summary(fit1_st4)
plot(fit1_st4)
wp(fit1_st4, ylim.all = 0.6)

### Utilizing the two strategies described in the book "Flexible Regression and Smoothing: Using GAMLSS in R", page 397, to get the best models for each parameter
#### Adjusting a model with only the intercept for mu and considering sigma and nu as constants 
fit2_st4 <- gamlss(
  sqrt(tx_fecundidade_menor_20) ~ 1,
  data = df_indicadores,
  family = ST4,
  control = gamlss.control(n.cyc = 500)
)
summary(fit2_st4)

#### Strategy 1:
#### Starting from the simplest model and utilizing the forward stepwise method to find the model with the lowest AIC for each parameter
#### Obs.: in this strategy, different variables can be selected for each parameter's model
fit3_st4 <- stepGAICAll.A(
  fit2_st4,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
summary(fit3_st4)

#### Strategy 2:
#### This strategy forces all os the models to select the same variables subsets 
fit4_st4 <- stepGAICAll.B(
  fit2_st4,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
summary(fit4_st4)  

#### Using the LR test to verify if the simpler models are better than the full model 
LR.test(fit3_st4, fit1_st4)  # The model selected by the first strategy ISN'T better than the full model 
LR.test(fit4_st4, fit1_st4)  # The model selected by the second strategy ISN'T better than the full model 
GAIC(fit1_st4, fit3_st4, fit4_st4)

#### Residual analysis of the selected model
plot(fit1_st4)
wp(fit1_st4, ylim.all = 0.6)

#### Final choice: ST4's full model.


# Interpreting the selected model ---------------------------------------------
## Estimated coefficients each of the parameter's model
summary(fit1_st4)
getSmo(fit1_st4)$coef

## Trying to better understand the impact of each variable on the fertility rate
### Creating a dataframe with the desired observations to predict
idhms <- c(as.numeric(quantile(df_indicadores$idhm, 0.25)), median(df_indicadores$idhm), as.numeric(quantile(df_indicadores$idhm, 0.75)))
anos <- 2018:2021
coberturas_ab <- seq(0, 100, by = 10)

df_newdata <- expand.grid(idhms, anos, coberturas_ab) |>
  arrange(Var1) |>
  mutate(pandemia = factor(ifelse(Var2 < 2020, "before", "after"), levels = c("before", "after")))

colnames(df_newdata) <- c("idhm", "ano", "cobertura_ab", "pandemia")

### Making the predictions 
predicoes <- predict(fit1_st4, newdata = df_newdata, type = "response")

df_newdata_completo <- df_newdata |>
  mutate(
    idhm = as.factor(idhm),
    predicoes = predicoes^2
  ) |>
  group_by(idhm, cobertura_ab, pandemia) |>
  summarise(predicoes = mean(predicoes))

### Creating a dataframe with the mean variations of the predicted fertility rates
### for each combination of HDI-M and the pandemic indicator when the primary care 
### coverage is increased by 10
df_variacao <- df_newdata_completo |>
  group_by(idhm, pandemia) |>
  mutate(variacao = (predicoes - lag(predicoes))) |>
  summarise(variacao_media = round(mean(variacao, na.rm = T), 3)) |>
  arrange(idhm, pandemia)
df_variacao

### Plotting a similar information
plot_variacao <- ggplot(
  data = df_newdata_completo,
  mapping = aes(x = cobertura_ab, y = predicoes, colour = idhm, linetype = pandemia)
) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_continuous(breaks = unique(df_newdata_completo$cobertura_ab)) +
  labs(
    x = "Primary care coverage",
    y = "Predicted fertility rate per thousand women aged 10 to 19",
    colour = "HDI-M",
    linetype = "COVID-19 pandemic"
  ) +
  theme_bw() +
  theme(legend.position = "top")
plot_variacao

### Exporting the plot
ggsave(
  "figures/Fig4.tiff", 
  plot_variacao,
  width = 7.5, height = 6, units = "in", 
  dpi = 600, compression = "lzw"
)
