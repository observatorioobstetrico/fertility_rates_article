library(dplyr)
library(tidyr)
library(janitor)
library(gamlss)
library(ggplot2)

# Reading and manipulating the data -------------------------------------------
df_indicadores <- read.csv("databases/data_fertility_rates_article.csv") |>
  filter(ano >= 2018 & ano <= 2021) |>
  mutate(
    codmunres = as.character(codmunres),
    tx_fecundidade_menor_20 = round(nvm_10_a_19 / pop_feminina_10_a_19 * 1000, 1),
    cobertura_ab = round(media_cobertura_ab / populacao_total * 100, 1),
    pandemia = factor(ifelse(ano < 2020, "before", "during"), levels = c("before", "during")),
    ano = as.factor(ano)
  ) |>
  drop_na() # There are 5 municipalities without the M-HDI information


# Adjusting the models --------------------------------------------------------
## For the ZAGA distribution --------------------------------------------------
### Adjusting the full model 
fit1_zaga <- gamlss(
  tx_fecundidade_menor_20 ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  data = df_indicadores,
  family = ZAGA(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
saveRDS(fit1_zaga, "r_objects/fit1_zaga.RDS")

summary(fit1_zaga)
saveRDS(summary(fit1_zaga, save = TRUE), "r_objects/sumario_fit1_zaga.RDS")

plot(fit1_zaga)
wp(fit1_zaga, ylim.all = 0.6)

### Utilizing the two strategies described in the book "Flexible Regression and Smoothing: Using GAMLSS in R", page 397, to get the best models for each parameter
#### Adjusting a model with only the intercept for mu and considering sigma and nu as constants 
fit2_zaga <- gamlss(
  tx_fecundidade_menor_20 ~ 1,
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
saveRDS(fit3_zaga, "r_objects/fit3_zaga.RDS")

summary(fit3_zaga)
saveRDS(summary(fit3_zaga, save = TRUE), "r_objects/sumario_fit3_zaga.RDS")

plot(fit3_zaga)
wp(fit3_zaga, ylim.all = 0.6)

#### Strategy 2:
#### This strategy forces all os the models to select the same variables subsets 
fit4_zaga <- stepGAICAll.B(
  fit2_zaga,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
saveRDS(fit4_zaga, "r_objects/fit4_zaga.RDS")

summary(fit4_zaga)
saveRDS(summary(fit4_zaga, save = TRUE), "r_objects/sumario_fit4_zaga.RDS")

plot(fit4_zaga)
wp(fit4_zaga, ylim.all = 0.6)

### Adjusting a third model with only the significant terms for mu
fit5_zaga <- gamlss(
  tx_fecundidade_menor_20 ~ random(ano) + pandemia + cobertura_ab + idhm + pandemia*idhm + cobertura_ab*idhm,
  data = df_indicadores,
  family = ZAGA(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  tau.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
saveRDS(fit5_zaga, "r_objects/fit5_zaga.RDS")

summary(fit5_zaga)
saveRDS(summary(fit5_zaga, save = TRUE), "r_objects/sumario_fit5_zaga.RDS")

plot(fit5_zaga)
wp(fit5_zaga, ylim.all = 0.6)

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
  tx_fecundidade_menor_20 ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  data = df_indicadores,
  family = ST4(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  tau.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
saveRDS(fit1_st4, "r_objects/fit1_st4.RDS")

summary(fit1_st4)
saveRDS(summary(fit1_st4, save = TRUE), "r_objects/sumario_fit1_st4.RDS")

plot(fit1_st4)
wp(fit1_st4, ylim.all = 0.6)

### Utilizing the two strategies described in the book "Flexible Regression and Smoothing: Using GAMLSS in R", page 397, to get the best models for each parameter
#### Adjusting a model with only the intercept for mu and considering sigma and nu as constants 
fit2_st4 <- gamlss(
  tx_fecundidade_menor_20 ~ 1,
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
saveRDS(fit3_st4, "r_objects/fit3_st4.RDS")

summary(fit3_st4)
saveRDS(summary(fit3_st4, save = TRUE), "r_objects/sumario_fit3_st4.RDS")

plot(fit3_st4)
wp(fit3_st4, ylim.all = 0.6)

#### Strategy 2:
#### This strategy forces all os the models to select the same variables subsets 
fit4_st4 <- stepGAICAll.B(
  fit2_st4,
  scope = list(lower =~ 1, upper =  ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm))
)
saveRDS(fit4_st4, "r_objects/fit4_st4.RDS")

summary(fit4_st4)  
saveRDS(summary(fit4_st4, save = TRUE), "r_objects/sumario_fit4_st4.RDS")

plot(fit4_st4)
wp(fit4_st4, ylim.all = 0.6)

### Adjusting a third model with only the significant terms for mu
fit5_st4 <- gamlss(
  tx_fecundidade_menor_20 ~ random(ano) + pandemia + cobertura_ab + idhm + pandemia*cobertura_ab + cobertura_ab*idhm,
  data = df_indicadores,
  family = ST4(),
  sigma.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  nu.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  tau.formula = ~ random(ano) + (pandemia + cobertura_ab + idhm) * (pandemia + cobertura_ab + idhm),
  control = gamlss.control(n.cyc = 500)
)
saveRDS(fit5_st4, "r_objects/fit5_st4.RDS")

summary(fit5_st4)
saveRDS(summary(fit5_st4, save = TRUE), "r_objects/sumario_fit5_st4.RDS")

plot(fit5_st4)
wp(fit5_st4, ylim.all = 0.6)

#### Using the LR test to verify if the simpler models are better than the full model 
LR.test(fit3_st4, fit1_st4)  # The model selected by the first strategy ISN'T better than the full model 
LR.test(fit4_st4, fit1_st4)  # The model selected by the second strategy IS better than the full model 
LR.test(fit5_st4, fit1_st4)  # The model selected by the third strategy IS better than the full model 
GAIC(fit1_st4, fit3_st4, fit4_st4, fit5_st4)

#### Residual analysis of the selected model
fit_final <- fit5_st4

plot(fit_final)
wp(fit_final, ylim.all = 0.6)

#### Final choice: ST4's final reduced model.


# Interpreting the selected model ---------------------------------------------
## Estimated coefficients each of the parameter's model
summary(fit_final)
getSmo(fit_final)$coef

## Trying to better understand the impact of each variable on the fertility rate
### Creating a dataframe with the desired observations to predict
idhms <- c(as.numeric(quantile(df_indicadores$idhm, 0.25)), median(df_indicadores$idhm), as.numeric(quantile(df_indicadores$idhm, 0.75)))
anos <- 2018:2021
coberturas_ab <- seq(0, 100, by = 10)

df_newdata <- expand.grid(idhms, anos, coberturas_ab) |>
  arrange(Var1) |>
  mutate(pandemia = factor(ifelse(Var2 < 2020, "before", "during"), levels = c("before", "during")))

colnames(df_newdata) <- c("idhm", "ano", "cobertura_ab", "pandemia")

### Making the predictions 
predicoes <- predict(fit_final, newdata = df_newdata, type = "response")

df_newdata_completo <- df_newdata |>
  mutate(
    idhm = as.factor(idhm),
    predicoes = predicoes
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

### Creating a dataframe with the mean variations of the predicted fertility rates
### for each value of HDI-M when the primary care coverage is increased by 10
df_variacoes_sem_pandemia <- df_newdata_completo |>
  group_by(idhm, pandemia) |>
  mutate(variacao = (predicoes - lag(predicoes))) |>
  ungroup() |>
  group_by(idhm) |>
  summarise(variacao_media = round(mean(variacao, na.rm = T), 3)) |>
  arrange(idhm)

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
  theme(legend.position = "bottom")
plot_variacao

### Exporting the plot
ggsave(
  "figures/Fig4.tiff", 
  plot_variacao,
  width = 7.5, height = 6, units = "in", 
  dpi = 600, compression = "lzw"
)