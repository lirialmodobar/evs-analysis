# Set working directory
setwd ("C:/Users/Belle/Documents/Belle - Nanosight/dados_raw")

nanosight_intersect_ev_pequena <- read.csv("nanosight_intersect_ev_pequena.csv") 

library(dplyr)

nanosight_intersect_ev_pequena_long <- nanosight_intersect_ev_pequena %>%
  mutate(
    wave = factor(wave, levels = c("t1", "t2")),
    dcmadep = factor(dcmadep, levels = c(FALSE, TRUE),
                     labels = c("Ausente", "Presente")),
    dcanyhk = factor(dcanyhk, levels = c(FALSE, TRUE),
                     labels = c("Ausente", "Presente")),
    subjectid = factor(subjectid)
  )

#Analise exploratoria
nanosight_intersect_ev_pequena_long %>%
  group_by(wave, dcmadep) %>%
  summarise(
    media = mean(concentracao_real, na.rm = TRUE),
    sd = sd(concentracao_real, na.rm = TRUE),
    n = n()
  )
#Grafico
library(ggplot2)
ggplot(nanosight_intersect_ev_pequena_long, aes(x = wave, y = concentracao_real,
                  color = dcmadep, group = dcmadep)) +
  stat_summary(fun = mean, geom = "line") +
  stat_summary(fun = mean, geom = "point") +
  theme_minimal()
#Modelo Linear simples
##Depressao
modelo_dep_lm <- lm(
  concentracao_real ~ wave * dcmadep,
  data = nanosight_intersect_ev_pequena_long
)
anova(modelo_dep_lm)
summary(modelo_dep_lm)
### A análise por modelo linear revelou efeito significativo da depressão sobre a concentração de vesículas extracelulares (F(1,177)=4.22, p=0.041), bem como interação significativa entre tempo e depressão (F(1,177)=7.26, p=0.008), indicando que a variação temporal da concentração difere de acordo com a presença de depressão. Enquanto indivíduos sem depressão apresentaram tendência à redução da concentração entre t1 e t2, indivíduos com depressão apresentaram aumento significativo ao longo do tempo.

##TDAH
modelo_tdah_lm <- lm(
  concentracao_real ~ wave * dcanyhk,
  data = nanosight_intersect_ev_pequena_long
)
anova(modelo_tdah_lm)
summary(modelo_tdah_lm)

## Ajuste para os dois transtornos juntos
modelo_ajustado <- lm(
  concentracao_real ~ wave * dcmadep + wave * dcanyhk,
  data = nanosight_intersect_ev_pequena_long
)
anova(modelo_ajustado)

## Verificando a robustez
lm(log(concentracao_real) ~ wave * dcmadep + wave * dcanyhk,
   data = nanosight_intersect_ev_pequena_long)
## Mais robusto
anova(lm(log(concentracao_real) ~ wave * dcmadep + wave * dcanyhk,
         data = nanosight_intersect_ev_pequena_long))
## Modelo misto, efeito aleatorio por subjectid
lmer(log(concentracao_real) ~ wave * dcmadep + wave * dcanyhk + (1 | subjectid),
     data = nanosight_intersect_ev_pequena_long)
##### Grafico
library(lme4)
library(lmerTest)
modelo_misto <- lmer(
  log(concentracao_real) ~ wave * dcmadep + wave * dcanyhk + (1 | subjectid),
  data = nanosight_intersect_ev_pequena_long
)
library(emmeans)
emm_dep <- emmeans(
  modelo_misto,
  ~ wave | dcmadep,
  type = "response"   # back-transforma do log
)
emm_dep
emm_dep_df <- as.data.frame(emm_dep)
library(ggplot2)
ggplot(emm_dep_df,
       aes(x = wave,
           y = response,
           color = dcmadep,
           group = dcmadep)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.1,
    linewidth = 0.8
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    x = "Tempo",
    y = "Concentração real de EVs pequenas",
    color = "Depressão",
    title = "Variação longitudinal da concentração de EVs pequenas",
    subtitle = "Estimativas do modelo linear misto ajustado para TDAH"
  ) +
  theme_minimal(base_size = 14)





