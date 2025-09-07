# modelo de análise para investigar, tempo (W1/W2), grupo (4 categorias) e interação tempo grupo e setar os contrastes após a análise com correção. Dá para fazer isso com modelo misto também.

#Carregar pacotes
library(lme4)
library(lmerTest)   # para valores-p
library(emmeans)    # para contrastes pós-hoc
library(ggplot2)    # para visualização
library(writexl)    # para salvar em excel
library(patchwork)
library(glmmTMB)
library(car)
library(openxlsx)

###Analise post hoc com contraste por tempo e grupo

##Size
# Modelo linear misto
mod_lin_mist_size <- lmer(tamanho_mean_average ~ wave * Trajetoria + (1 | subjectid), data = nanosight_intersect)
# ANOVA para efeitos fixos (Type III SS)
anova_mod_size <- anova(mod_lin_mist_size, type = 3)
df_anova_mod_size <- as.data.frame(anova_mod_size)
df_anova_mod_size$Effect <- rownames(df_anova_mod_size) #Adiciona a coluna com o nome dos efeitos
df_anova_mod_size <- df_anova_mod_size[, c("Effect", setdiff(names(df_anova_mod_size), "Effect"))] # Reordena para deixar "Effect" como primeira coluna
write_xlsx(df_anova_mod_size, path = "anova_modelo_misto_size.xlsx") # Salvar no Excel
# Resumo dos efeitos fixos do modelo (coeficientes)
summary_mod_size <- summary(mod_lin_mist_size)
df_fixed_effects_size <- as.data.frame(coef(summary_mod_size))
df_fixed_effects_size$Coefficient <- rownames(df_fixed_effects_size)
# Reorganizar colunas (Coeficiente primeiro)
df_fixed_effects_size <- df_fixed_effects_size[, c("Coefficient", setdiff(names(df_fixed_effects_size), "Coefficient"))]
#Médias marginais por grupo em cada tempo:
emmeans(mod_lin_mist_size, ~ Trajetoria | wave)
#Médias marginais por tempo em cada grupo:
emmeans(mod_lin_mist_size, ~ wave | Trajetoria)
#Comparações dentro de cada tempo:
pairs(emmeans(mod_lin_mist_size, ~ Trajetoria | wave), adjust = "tukey")
#Comparações dentro de cada grupo:
pairs(emmeans(mod_lin_mist_size, ~ wave | Trajetoria), adjust = "tukey")
#####Para salvar em excel
# 1. Calcular médias marginais (EMMeans)
emm_group_by_wave_size <- emmeans(mod_lin_mist_size, ~ Trajetoria | wave)
emm_wave_by_group_size <- emmeans(mod_lin_mist_size, ~ wave | Trajetoria)
# 2. Comparações (contrastes) entre níveis
pairs_group_by_wave_size <- pairs(emm_group_by_wave_size, adjust = "tukey")
pairs_wave_by_group_size <- pairs(emm_wave_by_group_size, adjust = "tukey")
# 3. Converter todos os resultados para data.frames
df_emm_group_by_wave_size <- as.data.frame(emm_group_by_wave_size)
df_emm_wave_by_group_size <- as.data.frame(emm_wave_by_group_size)
df_pairs_group_by_wave_size <- as.data.frame(pairs_group_by_wave_size)
df_pairs_wave_by_group_size <- as.data.frame(pairs_wave_by_group_size)
# 4. Escrever tudo em um arquivo Excel com múltiplas abas
write_xlsx(
  list(
    "Fixed_effects(coefficients)" = df_fixed_effects_size,
    "EMMeans_Trajectory_by_Wave" = df_emm_group_by_wave_size,
    "EMMeans_Wave_by_Trajectory" = df_emm_wave_by_group_size,
    "Contrasts_Trajectory_by_Wave" = df_pairs_group_by_wave_size,
    "Contrasts_Wave_by_Trajectory" = df_pairs_wave_by_group_size
  ),
  path = "results_emmeans_size.xlsx"
)

##concentration
mod_lin_mist_conc <- lmer(concentracao_real ~ wave * Trajetoria + (1 | subjectid), data = nanosight_intersect)
summary(mod_lin_mist_conc)
# ANOVA para efeitos fixos (Type III SS)
anova_mod_conc <- anova(mod_lin_mist_conc, type = 3)
df_anova_mod_conc <- as.data.frame(anova_mod_conc)
df_anova_mod_conc$Effect <- rownames(df_anova_mod_conc) #Adiciona a coluna com o nome dos efeitos
df_anova_mod_conc <- df_anova_mod_conc[, c("Effect", setdiff(names(df_anova_mod_conc), "Effect"))] # Reordena para deixar "Effect" como primeira coluna
write_xlsx(df_anova_mod_conc, path = "anova_modelo_misto_conc.xlsx") # Salvar no Excel
# Resumo dos efeitos fixos do modelo (coeficientes)
summary_mod_conc <- summary(mod_lin_mist_conc)
df_fixed_effects_conc <- as.data.frame(coef(summary_mod_conc))
df_fixed_effects_conc$Coefficient <- rownames(df_fixed_effects_conc)
# Reorganizar colunas (Coeficiente primeiro)
df_fixed_effects_conc <- df_fixed_effects_conc[, c("Coefficient", setdiff(names(df_fixed_effects_conc), "Coefficient"))]
#Médias marginais por grupo em cada tempo:
emmeans(mod_lin_mist_conc, ~ Trajetoria | wave)
#Médias marginais por tempo em cada grupo:
emmeans(mod_lin_mist_conc, ~ wave | Trajetoria)
#Comparações dentro de cada tempo:
pairs(emmeans(mod_lin_mist_conc, ~ Trajetoria | wave), adjust = "tukey")
#Comparações dentro de cada grupo:
pairs(emmeans(mod_lin_mist_conc, ~ wave | Trajetoria), adjust = "tukey")
#####Para salvar em excel
# 1. Calcular médias marginais (EMMeans)
emm_group_by_wave <- emmeans(mod_lin_mist_conc, ~ Trajetoria | wave)
emm_wave_by_group <- emmeans(mod_lin_mist_conc, ~ wave | Trajetoria)
# 2. Comparações (contrastes) entre níveis
pairs_group_by_wave <- pairs(emm_group_by_wave, adjust = "tukey")
pairs_wave_by_group <- pairs(emm_wave_by_group, adjust = "tukey")
# 3. Converter todos os resultados para data.frames
df_emm_group_by_wave <- as.data.frame(emm_group_by_wave)
df_emm_wave_by_group <- as.data.frame(emm_wave_by_group)
df_pairs_group_by_wave <- as.data.frame(pairs_group_by_wave)
df_pairs_wave_by_group <- as.data.frame(pairs_wave_by_group)
# 4. Escrever tudo em um arquivo Excel com múltiplas abas
write_xlsx(
  list(
    "Fixed_effects(coefficients)" = df_fixed_effects_conc,
    "EMMeans_Trajectory_by_Wave" = df_emm_group_by_wave,
    "EMMeans_Wave_by_Trajectory" = df_emm_wave_by_group,
    "Contrasts_Trajectory_by_Wave" = df_pairs_group_by_wave,
    "Contrasts_Wave_by_Trajectory" = df_pairs_wave_by_group
  ),
  path = "results_emmeans_concentration.xlsx"
)

#concentration sem interação com fator aleatorio
mod_lin_mist_conc_sem_interacao <- lmer(concentracao_real ~ wave + Trajetoria + (1 | subjectid), data = nanosight_intersect)
summary(mod_lin_mist_conc_sem_interacao)
# ANOVA para efeitos fixos (Type III SS)
anova_mod_conc_sem_interacao <- anova(mod_lin_mist_conc_sem_interacao, type = 3)
df_anova_mod_conc_sem_interacao <- as.data.frame(anova_mod_conc_sem_interacao)
df_anova_mod_conc_sem_interacao$Effect <- rownames(df_anova_mod_conc_sem_interacao) #Adiciona a coluna com o nome dos efeitos
df_anova_mod_conc_sem_interacao <- df_anova_mod_conc_sem_interacao[, c("Effect", setdiff(names(df_anova_mod_conc_sem_interacao), "Effect"))] # Reordena para deixar "Effect" como primeira coluna
write_xlsx(df_anova_mod_conc_sem_interacao, path = "anova_modelo_misto_conc_sem_interacao.xlsx") # Salvar no Excel
# Resumo dos efeitos fixos do modelo (coeficientes)
summary_mod_conc_sem_interacao <- summary(mod_lin_mist_conc_sem_interacao)
df_fixed_effects_conc_sem_interacao <- as.data.frame(coef(summary_mod_conc_sem_interacao))
df_fixed_effects_conc_sem_interacao$Coefficient <- rownames(df_fixed_effects_conc_sem_interacao)
# Reorganizar colunas (Coeficiente primeiro)
df_fixed_effects_conc_sem_interacao <- df_fixed_effects_conc_sem_interacao[, c("Coefficient", setdiff(names(df_fixed_effects_conc_sem_interacao), "Coefficient"))]
#Médias marginais por grupo em cada tempo:
emmeans(mod_lin_mist_conc_sem_interacao, ~ Trajetoria | wave)
#Médias marginais por tempo em cada grupo:
emmeans(mod_lin_mist_conc_sem_interacao, ~ wave | Trajetoria)
#Comparações dentro de cada tempo:
pairs(emmeans(mod_lin_mist_conc_sem_interacao, ~ Trajetoria | wave), adjust = "tukey")
#Comparações dentro de cada grupo:
pairs(emmeans(mod_lin_mist_conc_sem_interacao, ~ wave | Trajetoria), adjust = "tukey")
# Função para calcular o ICC de um modelo misto
calcular_ICC <- function(modelo) {
  variancias <- as.data.frame(VarCorr(modelo))
  var_aleatoria <- variancias$vcov[variancias$grp != "Residual"]
  var_residual  <- variancias$vcov[variancias$grp == "Residual"]
  ICC <- var_aleatoria / (var_aleatoria + var_residual)
  return(round(ICC, 3))
}
#Aplicação ao modelo de concentracao sem interação
icc_conc_sem_interacao <- calcular_ICC(mod_lin_mist_conc_sem_interacao)
cat("ICC (concentracao sem interacao):", icc_conc_sem_interacao, "\n")
###Calculo do AIC e BIC
AIC(mod_lin_mist_conc_sem_interacao)
BIC(mod_lin_mist_conc_sem_interacao)


#concentration sem interação e sem fator aleatorio
# Ajustar o modelo linear sem efeito aleatório
mod_lin_conc_sem_interacao_sem_random <- lm(concentracao_real ~ wave + Trajetoria, data = nanosight_intersect)
# Resumo do modelo (coeficientes)
summary(mod_lin_conc_sem_interacao_sem_random)
# ANOVA para efeitos fixos (Type III SS) - usando a função anova do pacote 'car'
# O Type III SS é mais adequado para modelos desbalanceados.
anova_mod_conc_sem_interacao_sem_random <- Anova(mod_lin_conc_sem_interacao_sem_random, type = "III")
df_anova_mod_conc_sem_interacao_sem_random <- as.data.frame(anova_mod_conc_sem_interacao_sem_random)
df_anova_mod_conc_sem_interacao_sem_random$Effect <- rownames(df_anova_mod_conc_sem_interacao_sem_random) # Adiciona a coluna com o nome dos efeitos
df_anova_mod_conc_sem_interacao_sem_random <- df_anova_mod_conc_sem_interacao_sem_random[, c("Effect", setdiff(names(df_anova_mod_conc_sem_interacao_sem_random), "Effect"))] # Reordena para deixar "Effect" como primeira coluna
# salvando anova no excel
write_xlsx(df_anova_mod_conc_sem_interacao, path = "anova_modelo_lin_conc_sem_interacao_sem_random.xlsx")
# Resumo dos efeitos fixos do modelo (coeficientes)
summary_mod_conc_sem_interacao_sem_random <- summary(mod_lin_conc_sem_interacao_sem_random)
df_fixed_effects_conc_sem_interacao_sem_random <- as.data.frame(summary_mod_conc_sem_interacao_sem_random$coefficients)
df_fixed_effects_conc_sem_interacao_sem_random$Coefficient <- rownames(df_fixed_effects_conc_sem_interacao_sem_random)
df_fixed_effects_conc_sem_interacao_sem_random <- df_fixed_effects_conc_sem_interacao_sem_random[, c("Coefficient", setdiff(names(df_fixed_effects_conc_sem_interacao_sem_random), "Coefficient"))]
# Médias marginais por grupo em cada tempo:
emmeans(mod_lin_conc_sem_interacao_sem_random, ~ Trajetoria | wave)
# Médias marginais por tempo em cada grupo:
emmeans(mod_lin_conc_sem_interacao_sem_random, ~ wave | Trajetoria)
# Comparações dentro de cada tempo:
pairs(emmeans(mod_lin_conc_sem_interacao_sem_random, ~ Trajetoria | wave), adjust = "tukey")
# Comparações dentro de cada grupo:
pairs(emmeans(mod_lin_conc_sem_interacao_sem_random, ~ wave | Trajetoria), adjust = "tukey")
### Cálculo do AIC e BIC
AIC(mod_lin_conc_sem_interacao_sem_random)
BIC(mod_lin_conc_sem_interacao_sem_random)

##percentage
mod_lin_mist_perc <- lmer(EV_pequenas_porcentagem ~ wave * Trajetoria + (1 | subjectid), data = nanosight_intersect)
summary(mod_lin_mist_perc)
# ANOVA para efeitos fixos (Type III SS)
anova_mod_perc <- anova(mod_lin_mist_perc, type = 3)
df_anova_mod_perc <- as.data.frame(anova_mod_perc)
df_anova_mod_perc$Effect <- rownames(df_anova_mod_perc) #Adiciona a coluna com o nome dos efeitos
df_anova_mod_perc <- df_anova_mod_perc[, c("Effect", setdiff(names(df_anova_mod_perc), "Effect"))] # Reordena para deixar "Effect" como primeira coluna
write_xlsx(df_anova_mod_perc, path = "anova_modelo_misto_perc.xlsx") # Salvar no Excel
# Resumo dos efeitos fixos do modelo (coeficientes)
summary_mod_perc <- summary(mod_lin_mist_perc)
df_fixed_effects_perc <- as.data.frame(coef(summary_mod_perc))
df_fixed_effects_perc$Coefficient <- rownames(df_fixed_effects_perc)
# Reorganizar colunas (Coeficiente primeiro)
df_fixed_effects_perc <- df_fixed_effects_perc[, c("Coefficient", setdiff(names(df_fixed_effects_perc), "Coefficient"))]
#Médias marginais por grupo em cada tempo:
emmeans(mod_lin_mist_perc, ~ Trajetoria | wave)
#Médias marginais por tempo em cada grupo:
emmeans(mod_lin_mist_perc, ~ wave | Trajetoria)
#Comparações dentro de cada tempo:
pairs(emmeans(mod_lin_mist_perc, ~ Trajetoria | wave), adjust = "tukey")
#Comparações dentro de cada grupo:
pairs(emmeans(mod_lin_mist_perc, ~ wave | Trajetoria), adjust = "tukey")
#####Para salvar em excel
# 1. Calcular médias marginais (EMMeans)
emm_group_by_wave_perc <- emmeans(mod_lin_mist_perc, ~ Trajetoria | wave)
emm_wave_by_group_perc <- emmeans(mod_lin_mist_perc, ~ wave | Trajetoria)
# 2. Comparações (contrastes) entre níveis
pairs_group_by_wave_perc <- pairs(emm_group_by_wave_perc, adjust = "tukey")
pairs_wave_by_group_perc <- pairs(emm_wave_by_group_perc, adjust = "tukey")
# 3. Converter todos os resultados para data.frames
df_emm_group_by_wave_perc <- as.data.frame(emm_group_by_wave_perc)
df_emm_wave_by_group_perc <- as.data.frame(emm_wave_by_group_perc)
df_pairs_group_by_wave_perc <- as.data.frame(pairs_group_by_wave_perc)
df_pairs_wave_by_group_perc <- as.data.frame(pairs_wave_by_group_perc)
# 4. Escrever tudo em um arquivo Excel com múltiplas abas
write_xlsx(
  list(
    "Fixed_effects(coefficients)" = df_fixed_effects_perc,
    "EMMeans_Trajectory_by_Wave" = df_emm_group_by_wave_perc,
    "EMMeans_Wave_by_Trajectory" = df_emm_wave_by_group_perc,
    "Contrasts_Trajectory_by_Wave" = df_pairs_group_by_wave_perc,
    "Contrasts_Wave_by_Trajectory" = df_pairs_wave_by_group_perc
  ),
  path = "results_emmeans_percentage.xlsx"
)

###Calculo do ICC
# Função para calcular o ICC de um modelo misto
calcular_ICC <- function(modelo) {
  variancias <- as.data.frame(VarCorr(modelo))
  var_aleatoria <- variancias$vcov[variancias$grp != "Residual"]
  var_residual  <- variancias$vcov[variancias$grp == "Residual"]
  ICC <- var_aleatoria / (var_aleatoria + var_residual)
  return(round(ICC, 3))
}
#Aplicação ao modelo de tamanho
icc_size <- calcular_ICC(mod_lin_mist_size)
cat("ICC (tamanho):", icc_size, "\n")
#Aplicação ao modelo de concentracao
icc_conc <- calcular_ICC(mod_lin_mist_conc)
cat("ICC (concentracao):", icc_conc, "\n")
#Aplicação ao modelo de concentracao
icc_perc <- calcular_ICC(mod_lin_mist_perc)
cat("ICC (porcentagem):", icc_perc, "\n")

###Calculo do AIC e BIC
AIC(mod_lin_mist_size)
BIC(mod_lin_mist_size)
AIC(mod_lin_mist_conc)
BIC(mod_lin_mist_conc)
AIC(mod_lin_mist_perc)
BIC(mod_lin_mist_perc)

###Modelagem tamanho médio das VEs
# Ajustar o GLMM com a família Gamma.
modelo_gamma_bobyqa_size <- glmer(tamanho_mean_average ~ wave * Trajetoria + (1 | subjectid),
                             data = nanosight_intersect,
                             family = Gamma(link = "log"),
                             control = glmerControl(optimizer = "bobyqa"))
# Tabela principal com o teste F para a interação (Wald Chi-Square ou F-test dependendo do pacote/versão)
anova(modelo_gamma_bobyqa_size)
# Obter os contrastes (comparações par a par)
emmeans(modelo_gamma_bobyqa_size, specs = pairwise ~ Trajetoria | wave, type = "response") # 'type = "response"' dá os resultados na escala original (nm)
# Calculo do ICC
library(performance)
# Calcular o ICC para o modelo de tamanho (glmer)
icc_tamanho_glmer <- performance::icc(modelo_gamma_bobyqa_size)
summary(icc_tamanho_glmer)

###Modelagem concentração
# Distribuição Binomial Negativa Tipo 1
modelo_nbinom1_conc <- glmmTMB(concentracao_real ~ wave + Trajetoria,
                               data = nanosight_intersect,
                               family = nbinom1(link = "log"))
# Para glmmTMB, a função anova() do pacote 'car' é recomendada
Anova(modelo_nbinom1_conc, type = "II") # Teste Wald Chi-Square
# Obter os contrastes
emmeans(modelo_nbinom1_conc, specs = pairwise ~ Trajetoria | wave, type = "response") # 'type = "response"' dá os resultados na escala original (partículas/mL)
#Salvar em excel
anova_nbinom1_conc <- Anova(modelo_nbinom1_conc, type = "II")
df_anova_nbinom1_conc <- as.data.frame(anova_nbinom1_conc)
contrastes_nbinom1_conc <- emmeans(modelo_nbinom1_conc, specs = pairwise ~ Trajetoria | wave, type = "response")
df_contrastes_nbinom1_conc <- summary(contrastes_nbinom1_conc)
nbinom1_conc <- createWorkbook()
addWorksheet(nbinom1_conc, "ANOVA")
writeData(nbinom1_conc, "ANOVA", df_anova_nbinom1_conc, rowNames = TRUE)
addWorksheet(nbinom1_conc, "Contrastes")
writeData(nbinom1_conc, "Contrastes", df_contrastes_nbinom1_conc)
saveWorkbook(nbinom1_conc, "resultados_nbinom1_conc.xlsx", overwrite = TRUE)


###Modelagem Porcentagem
# Pré-processamento: Converter percentagem (0-100) para proporção (0-1)
nanosight_intersect$percentage_prop <- nanosight_intersect$EV_pequenas_porcentagem / 100
# Ajustar o GLMM com a família Beta
modelo_perc_beta <- glmmTMB(percentage_prop ~ wave * Trajetoria + (1 | subjectid),
                            data = nanosight_intersect,
                            family = beta_family(link = "logit"))
# Obter os resultados
Anova(modelo_perc_beta, type = "II")
# Obter os contrastes
emmeans(modelo_perc_beta, specs = pairwise ~ Trajetoria | wave, type = "response") # 'type = "response"' dá os resultados na escala original (proporções 0-1)
# Salvar em excer
anova_perc_beta <- Anova(modelo_perc_beta, type = "II")
df_anova_perc_beta <- as.data.frame(anova_perc_beta)
emmeans_results_perc_beta <- emmeans(modelo_perc_beta,
                           specs = pairwise ~ Trajetoria | wave,
                           type = "response")
df_emmeans_perc_beta <- as.data.frame(emmeans_results_perc_beta$emmeans)
df_contrastes_perc_beta <- as.data.frame(emmeans_results_perc_beta$contrasts)
resultados_perc_beta_workbook <- createWorkbook()
addWorksheet(resultados_perc_beta_workbook, "ANOVA")
writeData(resultados_perc_beta_workbook, "ANOVA", df_anova_perc_beta, rowNames = TRUE)
addWorksheet(resultados_perc_beta_workbook, "Medias Estimadas")
writeData(resultados_perc_beta_workbook, "Medias Estimadas", df_emmeans_perc_beta)
addWorksheet(resultados_perc_beta_workbook, "Contrastes")
writeData(resultados_perc_beta_workbook, "Contrastes", df_contrastes_perc_beta)
saveWorkbook(resultados_perc_beta_workbook, "resultados_porcentagem_beta.xlsx", overwrite = TRUE)


###Calculo do AIC e BIC
AIC(modelo_gamma_bobyqa_size)
BIC(modelo_gamma_bobyqa_size)
AIC(modelo_nbinom1_conc)
BIC(modelo_nbinom1_conc)
AIC(modelo_perc_beta)
BIC(modelo_perc_beta)

#ICC GLMM
#size
icc_size_glmm <- calcular_ICC(modelo_gamma_bobyqa_size)
cat("ICC (tamanho):", icc_size_glmm, "\n")
#concentracao
icc_conc_glmm <- calcular_ICC(modelo_nbinom1_conc)
cat("ICC (concentracao):", icc_conc_glmm, "\n")
#Aplicação ao modelo de concentracao
icc_perc_glmm <- calcular_ICC(modelo_perc_beta)
cat("ICC (porcentagem):", icc_perc_glmm, "\n")

###Excel resultados GLMM
pacman::p_load(openxlsx)
# 1. Obter e preparar os resultados dos modelos
# Modelo de Tamanho
# Obter a tabela ANOVA
tabela_anova_tamanho <- as.data.frame(anova(modelo_gamma_bobyqa_size))
# Obter os resultados de emmeans
tabela_emmeans_tamanho <- as.data.frame(emmeans(modelo_gamma_bobyqa_size, specs = pairwise ~ Trajetoria | wave, type = "response")$contrasts)
# Modelo de Concentração
# Obter a tabela ANOVA
tabela_anova_conc <- as.data.frame(Anova(modelo_nbinom1_conc, type = "II"))
# Obter os resultados de emmeans
tabela_emmeans_conc <- as.data.frame(emmeans(modelo_nbinom1_conc, specs = pairwise ~ Trajetoria | wave, type = "response")$contrasts)
# Modelo de Porcentagem
# Obter a tabela ANOVA
tabela_anova_perc <- as.data.frame(Anova(modelo_perc_beta, type = "II"))
# Obter os resultados de emmeans
tabela_emmeans_perc <- as.data.frame(emmeans(modelo_perc_beta, specs = pairwise ~ Trajetoria | wave, type = "response")$contrasts)
# 2. Salvar os resultados em um arquivo Excel com múltiplas abas
# Cria um objeto Workbook
wb <- createWorkbook()
# Adicionar as tabelas como abas no Workbook
addWorksheet(wb, "Tamanho_ANOVA")
writeData(wb, "Tamanho_ANOVA", tabela_anova_tamanho, rowNames = TRUE)
addWorksheet(wb, "Tamanho_Contrastes")
writeData(wb, "Tamanho_Contrastes", tabela_emmeans_tamanho, rowNames = FALSE)
addWorksheet(wb, "Conc_ANOVA")
writeData(wb, "Conc_ANOVA", tabela_anova_conc, rowNames = TRUE)
addWorksheet(wb, "Conc_Contrastes")
writeData(wb, "Conc_Contrastes", tabela_emmeans_conc, rowNames = FALSE)
addWorksheet(wb, "Perc_ANOVA")
writeData(wb, "Perc_ANOVA", tabela_anova_perc, rowNames = TRUE)
addWorksheet(wb, "Perc_Contrastes")
writeData(wb, "Perc_Contrastes", tabela_emmeans_perc, rowNames = FALSE)
# Salva o arquivo Excel
saveWorkbook(wb, "resultados_modelos.xlsx", overwrite = TRUE)
# Mensagem de confirmação
cat("Resultados salvos com sucesso em 'resultados_modelos.xlsx'\n")


###Graficos visuais modelos lineares mistos
# Gráfico A) Size
gA <- emmip(mod_lin_mist_size, Trajetoria ~ wave, CIs = TRUE) +
  labs(title = "A) Interaction Time x Trajectory", y = "Mean Size (nm)", x = "Time") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  guides(linetype = guide_legend(title = "Trajectory:"),
         color = guide_legend(title = "Trajectory:"))
# Gráfico B) Concentration
gB <- emmip(mod_lin_mist_conc, Trajetoria ~ wave, CIs = TRUE) +
  labs(title = "B) Interaction Time x Trajectory", y = "Concentration (particles/mL)", x = "Time") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  guides(linetype = guide_legend(title = "Trajectory:"),
         color = guide_legend(title = "Trajectory:"))
# Gráfico C) Percentage
gC <- emmip(mod_lin_mist_perc, Trajetoria ~ wave, CIs = TRUE) +
  labs(title = "C) Interaction Time x Trajectory", y = "Small EVs (%)", x = "Time") +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash")) +
  guides(linetype = guide_legend(title = "Trajectory:"),
         color = guide_legend(title = "Trajectory:"))
# Combinar os gráficos
painel_emmip <- gA / gB / gC + plot_layout(guides = "collect") & theme(legend.position = "bottom")
# Salvar
ggsave("painel_emmip.png", painel_emmip, width = 8, height = 14, dpi = 300)


###Graficos descritivos
##Ggplot size
ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = tamanho_mean_average, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 2.2, alpha = 0.8) +
  labs(title = "EV's mean size for Trajectory and Time", 
       x = "Trajectory", 
       y = "Size (nm)", 
       fill = "Time point", 
       color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
ggsave("size_trajectory_time_readable.pdf", width = 9, height = 6)
##Ggplot concentration
ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = concentracao_real, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 2.2, alpha = 0.8) +
  labs(title = "EV's concentration size for Trajectory and Time", 
       x = "Trajectory", 
       y = "Concentration (particles/mL)", 
       fill = "Time point", 
       color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
ggsave("concentration_trajectory_time_readable.pdf", width = 9, height = 6)
##Ggplot percentage of small EVs
ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = EV_pequenas_porcentagem, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 2.2, alpha = 0.8) +
  labs(title = "Small EV's percentage for Trajectory and Time", 
       x = "Trajectory", 
       y = "Percentage (%)", 
       fill = "Time point", 
       color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
ggsave("percentage_trajectory_time_readable.pdf", width = 9, height = 6)

###Gráficos descritivos para salvar tudo junto em um pdf
# Gráfico A - Size
plot_a <- ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = tamanho_mean_average, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "A) EV's Mean Size", x = "Trajectory", y = "Size (nm)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Gráfico B - Concentration
plot_b <- ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = concentracao_real, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "B) EV's Concentration", x = "Trajectory", y = "Concentration (particles/mL)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Gráfico C - Small EV %
plot_c <- ggplot(nanosight_plus_sampleinfo, aes(x = Trajetoria, y = EV_pequenas_porcentagem, fill = wave)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA, alpha = 0.7, width = 0.6) +
  geom_jitter(aes(color = wave), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 2.2, alpha = 0.8) +
  labs(title = "C) Small EV's Percentage", x = "Trajectory", y = "Percentage (%)", fill = "Time point", color = "Time point") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40")
# Combinar todos em um painel
painel_final <- (plot_a / plot_b / plot_c) + plot_layout(guides = "collect")
# Salvar em um único arquivo PDF
ggsave("painel_EV_boxplots.pdf", painel_final, width = 11, height = 14)
# Salvar em png
ggsave("painel_EV_boxplots.png", painel_final, width = 11, height = 14, dpi = 300)

###Testes ANOVA sem efeito aleatorio
mod_lm_size <- lm(tamanho_mean_average ~ wave * Trajetoria, data = nanosight_intersect)
anova(mod_lm_size)
mod_lm_conc <- lm(concentracao_real ~ wave * Trajetoria, data = nanosight_intersect)
anova(mod_lm_conc)
mod_lm_perc <- lm(EV_pequenas_porcentagem ~ wave * Trajetoria, data = nanosight_intersect)
anova(mod_lm_perc)
