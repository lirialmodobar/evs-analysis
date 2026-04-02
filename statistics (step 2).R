# SET WORKING DIRECTORY
project_directory <- "C:/Users/lirie/belle"
tryCatch({ setwd(project_directory); cat("Working directory successfully set to:\n", getwd(), "\n\n")},
         error = function(e) { stop("ERROR: The specified directory was not found. Please check the path.") })

# ENVIRONMENT AND PACKAGES
cat("--- SECTION 0: SETTING UP THE ENVIRONMENT ---\n")
pacotes_necessarios <- c("lme4", "lmerTest", "glmmTMB", "car", "emmeans",
                         "performance", "DHARMa", "gtsummary", "dplyr", "tidyr", "correlation", "see")
for (pacote in pacotes_necessarios) {
  if (!require(pacote, character.only = TRUE)) {
    install.packages(pacote, dependencies = TRUE); library(pacote, character.only = TRUE)
  }
}
if (!dir.exists("analysis_results")) { dir.create("analysis_results") }
cat("Packages loaded and results directory is ready.\n\n")

options(warn = 1)

transtorno_tabela <- list(
  list(arquivo = "nanosight_intersect_ev_pequena_cross.csv", coluna_trajetoria = "Trajetoria",tag_transtorno = "cross", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"),
  list(arquivo = "nanosight_intersect_p90_cross.csv", coluna_trajetoria = "Trajetoria", tag_transtorno = "cross", coluna_perc = "p90_porcentagem", tag_perc = "P90"),
  list(arquivo = "nanosight_intersect_p95_cross.csv", coluna_trajetoria = "Trajetoria", tag_transtorno = "cross", coluna_perc = "p95_porcentagem", tag_perc = "P95"),
  list(arquivo = "nanosight_intersect_ev_pequena_dep.csv", coluna_trajetoria = "grupo_analise_dep", tag_transtorno = "dep", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"),
  list(arquivo = "nanosight_intersect_p90_dep.csv", coluna_trajetoria = "grupo_analise_dep", tag_transtorno = "dep", coluna_perc = "p90_porcentagem", tag_perc = "P90"),
  list(arquivo = "nanosight_intersect_p95_dep.csv", coluna_trajetoria = "grupo_analise_dep", tag_transtorno = "dep", coluna_perc = "p95_porcentagem", tag_perc = "P95"),
  list(arquivo = "nanosight_intersect_ev_pequena_anx_inc_control.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p90_anx_inc_control.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control", coluna_perc = "p90_porcentagem", tag_perc = "P90"),#n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p95_anx_inc_control.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control", coluna_perc = "p95_porcentagem", tag_perc = "P95"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_ev_pequena_adhd_inc_control.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p90_adhd_inc_control.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control", coluna_perc = "p90_porcentagem", tag_perc = "P90"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p95_adhd_inc_control.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control", coluna_perc = "p95_porcentagem", tag_perc = "P95") ,#n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_ev_pequena_anx_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control_remit", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p90_anx_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control_remit", coluna_perc = "p90_porcentagem", tag_perc = "P90"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p95_anx_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_gena", tag_transtorno = "anx_inc_control_remit", coluna_perc = "p95_porcentagem", tag_perc = "P95"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_ev_pequena_adhd_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control_remit", coluna_perc = "EV_pequenas_porcentagem", tag_perc = "EV_PEQUENA"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p90_adhd_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control_remit", coluna_perc = "p90_porcentagem", tag_perc = "P90"), #n baixo, so incidente vs controle
  list(arquivo = "nanosight_intersect_p95_adhd_inc_control_remit.csv", coluna_trajetoria = "grupo_analise_anyhk", tag_transtorno = "adhd_inc_control_remit", coluna_perc = "p95_porcentagem", tag_perc = "P95") #n baixo, so incidente vs controle
 ) ##Para testar transtornos individualmente e cross tb ver se as maiores EVs e ver se n tem influencia apenas do dado ser muito em um extremo

# Tabela global para acumular resultados de LRT
tabela_global_lrt <- data.frame()
# HELPER FUNCTIONS
# CREATE DHARMa diagonosis table
criar_tabela_diagnostico_dharma_qq <- function(simulationOutput) {
  alfa <- 0.05
  teste_uniformidade <- testUniformity(simulationOutput, plot = FALSE); teste_dispersao <- testDispersion(simulationOutput, plot = FALSE); teste_outliers <- testOutliers(simulationOutput, plot = FALSE)
  tabela <- data.frame(Diagnostic_Test = c("Overall Uniformity (KS)", "Dispersion", "Outliers"), p_value = c(round(teste_uniformidade$p.value, 3), round(teste_dispersao$p.value, 3), round(teste_outliers$p.value, 3)))
  tabela$Significance <- ifelse(tabela$p_value < alfa, "Significant Violation", "OK")
  return(tabela)
}

vis <- c("coluna_traj", "wave", "sex", "bage", "site")

for (transtorno in transtorno_tabela) {
  
  arquivo_csv <- transtorno$arquivo
  tag_perc <- transtorno$tag_perc
  tag_transtorno <- transtorno$tag_transtorno
  coluna_perc <- transtorno$coluna_perc
  coluna_traj <- transtorno$coluna_trajetoria
  if (!dir.exists(paste0("analysis_results/", tag_perc, "_", tag_transtorno))) { dir.create(paste0("analysis_results/", tag_perc,  "_", tag_transtorno)) }
  
  # LOG
  log_file_name <- paste0("analysis_results/", tag_perc,"_", tag_transtorno, "/log_analise_", tag_perc, "_",  tag_transtorno, "_", Sys.Date(), ".txt")
  ## ABRINDO CONEXÃO DE LOG
  con <- file(log_file_name, open = "wt")
  sink(con, split = TRUE)            # Saída padrão (cat, print) vai para Console + Log
  sink(con, type = "message")        # Mensagens de Erro/Aviso vão APENAS para o Log (limpa o console)
  
  cat("### INÍCIO DO LOG DE ANÁLISE - GRUPO:", tag_transtorno, " ###\n")
  cat("Data e Hora:", as.character(Sys.time()), "\n")
  cat("Diretório de Trabalho:", getwd(), "\n\n")
  cat("Arquivo de Origem:", arquivo_csv, "\n")
  cat("Coluna de porcentagem Analisada:", coluna_perc, "\n\n")
  cat("Transtorno analisado:", tag_transtorno, "\n\n")
  
  
  tryCatch({
    
    # LOADING AND PREPARING DATA
    cat("--- SECTION 1: LOADING AND PREPARING DATA ---\n")
    tryCatch({
      nanosight_intersect <- read.csv(arquivo_csv)
      cat("File ", arquivo_csv, " loaded successfully.\n")
      nanosight_intersect[[coluna_traj]] <- as.factor(nanosight_intersect[[coluna_traj]])
      nanosight_intersect$wave <- as.factor(nanosight_intersect$wave)
      nanosight_intersect$subjectid <- as.factor(nanosight_intersect$subjectid)
      nanosight_intersect$sex <- as.factor(nanosight_intersect$sex)
      nanosight_intersect$site <- as.factor(nanosight_intersect$site)
      nanosight_intersect$coluna_traj <- nanosight_intersect[[coluna_traj]]
      nanosight_intersect$bage <- as.numeric(nanosight_intersect$bage)
      tabela_desenho <- nanosight_intersect %>%
        mutate(grupo = interaction(wave, coluna_traj, drop = TRUE)) %>%
        select(grupo, sex, site, bage) %>%
        tbl_summary(
          by = grupo,
          statistic = all_continuous() ~ "{mean} ± {sd}",
          digits = all_continuous() ~ 2
        ) %>%
        add_ci() %>%
        add_n()
      tabela_desenho_print <- tabela_desenho
      tabela_desenho_print[is.na(tabela_desenho_print)] <- ""
      print(as.data.frame(as_tibble(tabela_desenho_print)))
      gtsummary::as_flex_table(tabela_desenho) |>
        flextable::save_as_docx(
          path = paste0(
            "analysis_results/", tag_perc, "_", tag_transtorno,
            "/tabela_desenho_", tag_transtorno, "_", tag_perc, ".docx"
          )
        )
     # write.csv(tabela_desenho, paste0("analysis_results/", tag_transtorno, "_", tag_perc, "/tabela_desenho_", tag_transtorno, "_", tag_perc), row.names = F, quote = F)
      nanosight_intersect$bage <- as.numeric(scale(nanosight_intersect$bage))
    }, error = function(e) {
      message("ERRO: ", e$message)
      stop(e)
    })
    cat("Verifying and transforming the percentage variable...\n")
    nanosight_intersect$coluna_porcentagem_em_uso <- nanosight_intersect[[coluna_perc]]
    proporcao_pura <- nanosight_intersect$coluna_porcentagem_em_uso / 100
    if (any(nanosight_intersect$coluna_porcentagem_em_uso == 0, na.rm = TRUE) || any(nanosight_intersect$coluna_porcentagem_em_uso == 100, na.rm = TRUE)) {
      cat("Values of 0 or 100 detected. Applying transformation for the Beta model: (valor/100)*((n-1) + 0.5)/n.\n")
      n <- nrow(nanosight_intersect)
      nanosight_intersect$percentage_prop <- (proporcao_pura * (n - 1) + 0.5) / n
    } else {
      cat("No 0 or 100 values detected. Using simple division by 100.\n")
      nanosight_intersect$percentage_prop <- proporcao_pura
    }
    cat("Data preparation complete.\n\n")
    
    resultados_finais_aic_bic <- data.frame()
    alfa <- 0.05
    
  
    # MODEL 1: EV's Sizes
    cat("\n\n=========================================================\n"); cat("  ANALYSIS 1: EV MEAN SIZE (tamanho_mode_average)\n"); cat("=========================================================\n\n")
    cat("--- 1.1 Gaussian model ---\n")
    if ((tag_perc == "P95" && tag_transtorno %in% c("adhd_inc_control_remit", "dep")) || (tag_perc == "P90" && tag_transtorno %in% c("dep", "adhd_inc_control_remit"))) {
      equacao_gaussiano <- "glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj)"
      cat("Gaussiano: ", equacao_gaussiano,"\n"  )
      modelo_size_random <- tryCatch({glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj)}, warning = function(w) {cat("WARNING\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
      modelo_size_vif <- glmmTMB(tamanho_mode_average ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj)
      modelo_size_no_random <- glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj)
      if(tag_perc %in% c("P90", "P95") && tag_transtorno == "dep") {
        equacao_gaussiano <- "glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj)"
        cat("Gaussiano final: ", equacao_gaussiano,"\n"  )
        modelo_size <- modelo_size_no_random
        size_no_random_sem_disp <- glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = gaussian())
        lrt_size_fixo <- anova(modelo_size, size_no_random_sem_disp)
        cat("\n LRT FIXO \n"); print(lrt_size_fixo)
        tabela_global_lrt <- rbind(tabela_global_lrt, data.frame(
          Disorder = tag_transtorno, Percentage_threshold = tag_perc, Dependent_Variable = "Size FIXO", 
          LRT_p = round(lrt_size_fixo$`Pr(>Chisq)`[2], 4), 
          Decision = ifelse(lrt_size_fixo$`Pr(>Chisq)`[2] < 0.05, "Disp", "No disp")
        ))
        cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
        residuos_nodisp <- simulateResiduals(fittedModel = size_no_random_sem_disp); tabela_diag_nodisp <- criar_tabela_diagnostico_dharma_qq(residuos_nodisp)
        png(paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_size_FIXO_NO_DISP_GAUSSIAN_", tag_perc, "_",  tag_transtorno, ".png"), width = 3000, height = 1200, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 0.8); plotQQunif(residuos_nodisp, main = "A) QQ Plot - Size FIXO SEM DISP", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_nodisp, main = ""); title("B) Residuals vs. Predicted – Size FIXO SEM DISP", line = 2); dev.off()
        cat("AIC"); print(AIC(size_no_random_sem_disp))
        cat("BIC"); print(BIC(size_no_random_sem_disp))
        } else {
        modelo_size <- modelo_size_random
        cat("\nVarrCorr:\n"); print(VarCorr(modelo_size)) #disp nao vai ICC
      }
      } else if (tag_perc %in% c("P95", "P90") && tag_transtorno == "adhd_inc_control"){
        equacao_gaussiano <- "glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj * bage)"
        modelo_size_random <- tryCatch({glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj * bage)}, warning = function(w) {cat("WARNING\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
        modelo_size_vif <- glmmTMB(tamanho_mode_average ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj * bage)
        modelo_size_no_random <- glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = gaussian(), dispformula = ~ coluna_traj + bage)
        modelo_size  <- modelo_size_random
        cat("\nVarrCorr:\n"); print(VarCorr(modelo_size))
        } else {
        equacao_gaussiano <- "glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian())"
        cat("Gaussiano: ", equacao_gaussiano,"\n"  )
        modelo_size_random <- tryCatch({glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian())}, warning = function(w) {cat("WARNING\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
        modelo_size_vif <- glmmTMB(tamanho_mode_average ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = gaussian())
        modelo_size_no_random <- glmmTMB(tamanho_mode_average ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = gaussian())
        modelo_size <- modelo_size_random
        cat("\nVarrCorr:\n"); print(VarCorr(modelo_size))
        icc_size <- icc(modelo_size)
        cat("\nICC:\n"); print(icc_size)
        }
    cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
    residuos_gaussiano <- simulateResiduals(fittedModel = modelo_size); tabela_diag_gaussiano <- criar_tabela_diagnostico_dharma_qq(residuos_gaussiano)
    png(paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_size_GAUSSIAN_", tag_perc, "_",  tag_transtorno, ".png"), width = 3000, height = 1200, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 0.8); plotQQunif(residuos_gaussiano, main = "A) QQ Plot - Size (Gaussian)", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_gaussiano, main = ""); title("B) Residuals vs. Predicted – Size (Gaussian)", line = 2); dev.off()
    # Extract only the rows actually used in fitting
    df_model <- model.frame(modelo_size)
    vars <- c("sex", "site", "wave", "bage", coluna_traj)
    
    for(v in vars){
      # Define filename for each plot
      filename <-paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_size_GAUSSIAN_", v, "_", tag_perc, "_",  tag_transtorno, ".png")
      
      # Open PNG device
      png(filename)
      
      # Plot residuals
      plotResiduals(residuos_gaussiano, df_model[[v]], main = paste("Residuals vs", v))
      
      # Close the device to save the file
      dev.off()
    }
    cat("\nDiagnostics for Gaussian Model:\n"); print(tabela_diag_gaussiano); write.csv(tabela_diag_gaussiano, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_diagnostics_dharma_size_GAUSSIAN_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Size", Model = "Gaussian", AIC = AIC(modelo_size), BIC = BIC(modelo_size)))
    cat("\n--- 1.3. ANOVA Table (Type III) for Size ---\n")
    anova_size <- Anova(modelo_size, type = "III"); df_anova_size <- as.data.frame(anova_size); df_anova_size$Significance <- ifelse(df_anova_size$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_size); write.csv(df_anova_size, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_anova_size_", tag_perc, "_", tag_transtorno, ".csv"))
    cat("\n--- 1.4. Conditional Pairwise Comparisons (Size) ---\n")
    p_interacao_size <- df_anova_size["wave:coluna_traj", "Pr(>Chisq)"]
    if (p_interacao_size < 0.05) {
      cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
      emm_size_p1 <- emmeans(modelo_size, specs = pairwise ~ coluna_traj | wave, adjust = "tukey"); df_emm_size_p1 <- as.data.frame(emm_size_p1$contrasts); df_emm_size_p1$Significance <- ifelse(df_emm_size_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p1); write.csv(df_emm_size_p1, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
      emm_size_p2 <- emmeans(modelo_size, specs = pairwise ~ wave | coluna_traj, adjust = "tukey"); df_emm_size_p2 <- as.data.frame(emm_size_p2$contrasts); df_emm_size_p2$Significance <- ifelse(df_emm_size_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_p2); write.csv(df_emm_size_p2, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    } else {
     cat("Interaction is not significant. Proceeding to analyze main effects with a Type II ANOVA.\n")
      anova_size_type2 <- Anova(modelo_size, type = "II"); print(anova_size_type2)
      if (anova_size_type2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_size_wave <- emmeans(modelo_size, specs = pairwise ~ wave, adjust = "tukey"); df_emm_size_wave <- as.data.frame(emm_size_wave$contrasts); df_emm_size_wave$Significance <- ifelse(df_emm_size_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_wave); write.csv(df_emm_size_wave, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
      if (anova_size_type2["coluna_traj", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'coluna_traj' is significant. Performing contrasts.\n\n"); emm_size <- emmeans(modelo_size, specs = pairwise ~ coluna_traj, adjust = "tukey"); df_emm_size <- as.data.frame(emm_size$contrasts); df_emm_size$Significance <- ifelse(df_emm_size$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size); write.csv(df_emm_size, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_size", tag_transtorno,".csv"), row.names = FALSE)}
      if (anova_size_type2["sex", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'sex' is significant. Performing contrasts.\n\n"); emm_size_sex <- emmeans(modelo_size, specs = pairwise ~ sex, adjust = "tukey"); df_emm_size_sex <- as.data.frame(emm_size_sex$contrasts); df_emm_size_sex$Significance <- ifelse(df_emm_size_sex$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_sex); write.csv(df_emm_size_sex, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_sex_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
      if (anova_size_type2["site", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'site' is significant. Performing contrasts.\n\n"); emm_size_site <- emmeans(modelo_size, specs = pairwise ~ site, adjust = "tukey"); df_emm_size_site <- as.data.frame(emm_size_site$contrasts); df_emm_size_site$Significance <- ifelse(df_emm_size_site$p.value < alfa, "Significant", "Not Significant"); print(df_emm_size_site); write.csv(df_emm_size_site, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_site_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    }
    cat("\n--- 1.5. Model Diagnostics (Size) ---\n")
    cat("Assessing VIF:\n")
    vif_df_size <- as.data.frame(check_collinearity(modelo_size_vif)); vif_df_size$Interpretation <- ifelse(vif_df_size$VIF >= 10, "High", ifelse(vif_df_size$VIF >= 5, "Moderate", "Low"))
    cat("VIF:\n"); print(vif_df_size); write.csv(vif_df_size, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_vif_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    #Random ou nao
    lrt_size <- anova(modelo_size_random, modelo_size_no_random)
    cat("\n LRT \n"); print(lrt_size)
    tabela_global_lrt <- rbind(tabela_global_lrt, data.frame(
      Disorder = tag_transtorno, Percentage_threshold = tag_perc, Dependent_Variable = "Size", 
      LRT_p = round(lrt_size$`Pr(>Chisq)`[2], 4), 
      Decision = ifelse(lrt_size$`Pr(>Chisq)`[2] < 0.05, "Random factor", "No random factor")
    ))
    
    # MODEL 2: EV's concentration
    cat("\n\n====================================================================\n"); cat("  ANALYSIS 2: EV CONCENTRATION (concentracao_real)\n"); cat("====================================================================\n\n")
   if (tag_transtorno == "anx_inc_control") {
      equacao_conc <- "glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = log), dispformula = ~ bage)"
      cat("Gamma GLMM...\n", "equacao: ", equacao_conc, "\n", "\n")
      modelo_conc_random <- tryCatch({glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ bage)}, warning = function(w) {cat("CONVERGENCE WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
      modelo_conc_vif <- glmmTMB(concentracao_real ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ bage)
      modelo_conc_no_random <- glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ bage)
      modelo_conc <- modelo_conc_random
      cat("\nVarrCorr:\n"); print(VarCorr(modelo_conc))
      } else if (tag_transtorno == "adhd_inc_control_remit" && tag_perc %in% c("P95", "P90")) {
      equacao_conc <- "glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = log), dispformula = ~ coluna_traj)"
      cat("Gamma GLMM...\n", "equacao: ", equacao_conc, "\n", "\n")
      modelo_conc_random <- tryCatch({glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ coluna_traj)}, warning = function(w) {cat("CONVERGENCE WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
      modelo_conc_vif <- glmmTMB(concentracao_real ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ coluna_traj)
      modelo_conc_no_random <- glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = Gamma(link = "log"), dispformula = ~ coluna_traj)
      modelo_conc <- modelo_conc_random
      cat("\nVarrCorr:\n"); print(VarCorr(modelo_conc))
       } else {
        equacao_conc <- "glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = log))"
        cat("Gamma GLMM...\n", "equacao: ", equacao_conc, "\n", "\n")
        modelo_conc_random <- tryCatch({glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log"))}, warning = function(w) {cat("CONVERGENCE WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("CONVERGENCE ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
        modelo_conc_vif <- glmmTMB(concentracao_real ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = Gamma(link = "log")) 
        modelo_conc_no_random <- glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = Gamma(link = "log")) 
        if (tag_transtorno %in% c("adhd_inc_control", "anx_inc_control_remit")) {
          modelo_conc <- modelo_conc_no_random
          equacao_conc <- "glmmTMB(concentracao_real ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = Gamma(link = log))"
          cat("Gamma GLMM Conc...\n", "equacao final: ", equacao_conc, "\n", "\n")
        } else {
          modelo_conc <- modelo_conc_random
          cat("\nVarrCorr:\n"); print(VarCorr(modelo_conc))
          cat("\nICC:\n"); print(icc(modelo_conc))
        }
        }
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Concentration", Model = "Gamma", AIC =AIC(modelo_conc), BIC = BIC(modelo_conc)))
    cat("\n--- 1.3. ANOVA Table (Type III) for Concentration ---\n")
    anova_conc <- Anova(modelo_conc, type = "III"); df_anova_conc <- as.data.frame(anova_conc); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_anova_conc_", tag_perc, "_", tag_transtorno, ".csv"))
    cat("\n--- 1.4. Conditional Pairwise Comparisons (Conc) ---\n")
    p_interacao_conc <- df_anova_conc["wave:coluna_traj", "Pr(>Chisq)"]
    if (p_interacao_conc < 0.05) {
      cat("Interaction is significant. Performing contrasts for BOTH interaction perspectives.\n")
      emm_conc_p1 <- emmeans(modelo_conc, specs = pairwise ~ coluna_traj | wave, adjust = "tukey"); df_emm_conc_p1 <- as.data.frame(emm_conc_p1$contrasts); df_emm_conc_p1$Significance <- ifelse(df_emm_conc_p1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p1); write.csv(df_emm_conc_p1, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_conc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
      emm_conc_p2 <- emmeans(modelo_conc, specs = pairwise ~ wave | coluna_traj, adjust = "tukey"); df_emm_conc_p2 <- as.data.frame(emm_conc_p2$contrasts); df_emm_conc_p2$Significance <- ifelse(df_emm_conc_p2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_p2); write.csv(df_emm_conc_p2, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_conc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    } else {
      cat("Interaction is not significant. Proceeding to analyze main effects with a Type II ANOVA.\n")
      anova_conc_type2 <- Anova(modelo_conc, type = "II"); df_anova_conc <- as.data.frame(anova_conc_type2); df_anova_conc$Significance <- ifelse(df_anova_conc$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_conc); write.csv(df_anova_conc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_anova_concentration_", tag_perc, "_", tag_transtorno, ".csv"))
      if (anova_conc_type2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_conc_wave <- emmeans(modelo_conc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_conc_wave <- as.data.frame(emm_conc_wave$contrasts); df_emm_conc_wave$Significance <- ifelse(df_emm_conc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_wave); write.csv(df_emm_conc_wave, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_size_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
      if (anova_conc_type2["coluna_traj", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'coluna_traj' is significant. Performing contrasts.\n\n"); emm_conc <- emmeans(modelo_conc, specs = pairwise ~ coluna_traj, adjust = "tukey"); df_emm_conc <- as.data.frame(emm_conc$contrasts); df_emm_conc$Significance <- ifelse(df_emm_conc$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc); write.csv(df_emm_conc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_conc", tag_transtorno,".csv"), row.names = FALSE)}
      if (anova_conc_type2["sex", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'sex' is significant. Performing contrasts.\n\n"); emm_conc_sex <- emmeans(modelo_conc, specs = pairwise ~ sex, adjust = "tukey"); df_emm_conc_sex <- as.data.frame(emm_conc_sex$contrasts); df_emm_conc_sex$Significance <- ifelse(df_emm_conc_sex$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_sex); write.csv(df_emm_conc_sex, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_sex_conc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
      if (anova_conc_type2["site", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'site' is significant. Performing contrasts.\n\n"); emm_conc_site <- emmeans(modelo_conc, specs = pairwise ~ site, adjust = "tukey"); df_emm_conc_site <- as.data.frame(emm_conc_site$contrasts); df_emm_conc_site$Significance <- ifelse(df_emm_conc_site$p.value < alfa, "Significant", "Not Significant"); print(df_emm_conc_site); write.csv(df_emm_conc_site, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_site_conc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    }  
    cat("\n--- 2.4. Model Diagnostics (Concentration) ---\n")
    vif_df_conc <- as.data.frame(check_collinearity(modelo_conc_vif)); vif_df_conc$Interpretation <- ifelse(vif_df_conc$VIF >= 10, "High", ifelse(vif_df_conc$VIF >= 5, "Moderate", "Low"));
    cat("VIF:\n"); print(vif_df_conc); write.csv(vif_df_conc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_vif_concentration_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE);
    #Random ou nao
    lrt_conc <- anova(modelo_conc_random, modelo_conc_no_random)
    cat("\n LRT \n"); print(lrt_conc)
    tabela_global_lrt <- rbind(tabela_global_lrt, data.frame(
      Disorder = tag_transtorno, Percentage_threshold = tag_perc, Dependent_Variable = "Concentration", 
      LRT_p = round(lrt_conc$`Pr(>Chisq)`[2], 4), 
      Decision = ifelse(lrt_conc$`Pr(>Chisq)`[2] < 0.05, "Random factor", "No random factor")
    ))
    #Dharma
    cat("\n--- 2.5. DHARMa Residual Diagnostics (Concentration) ---\n")
    residuos_conc <- simulateResiduals(fittedModel = modelo_conc); tabela_diag_conc <- criar_tabela_diagnostico_dharma_qq(residuos_conc)
    png(paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_concentration_", tag_perc, "_",  tag_transtorno, ".png"), width = 3000, height = 1200, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 0.8); plotQQunif(residuos_conc, main = "A) QQ Plot - Concentration", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_conc, main = ""); title("B) Residuals vs. Predicted - Concentration", line = 2); dev.off()
    df_model <- model.frame(modelo_conc)
    vars <- c("sex", "site", "wave", "bage", coluna_traj)
    
    for(v in vars){
      # Define filename for each plot
      filename <-paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_conc_", v, "_",  tag_perc, "_",  tag_transtorno, ".png")
      
      # Open PNG device
      png(filename)
      
      # Plot residuals
      plotResiduals(residuos_conc, df_model[[v]], main = paste("Residuals vs", v))
      
      # Close the device to save the file
      dev.off()
    }
    print(tabela_diag_conc); write.csv(tabela_diag_conc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_diagnostics_dharma_concentration_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    
    # MODEL 3: EVs percentages
    cat("\n\n======================================================================\n"); cat("  ANALYSIS 3: percentage OF SMALL EVs \n"); cat("======================================================================\n\n")
   if (tag_perc == "P90" && tag_transtorno == "dep")  {
     equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit), dispformula = ~ coluna_traj + site"
     cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao: ", equacao_perc, "\n")
     modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site)}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
     modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site)
     modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site)
     modelo_perc <- modelo_perc_random
     cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
     } else if (tag_perc == "EV_PEQUENA" && tag_transtorno == "adhd_inc_control_remit") {
     equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit), dispformula = ~ coluna_traj + site + bage"
     cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao: ", equacao_perc, "\n")
     modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site + bage)}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
     modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site + bage)
     modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + site + bage)
     modelo_perc <- modelo_perc_random
     cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
     } else if (tag_perc == "EV_PEQUENA" && tag_transtorno == "adhd_inc_control"){
     equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit), dispformula = ~ coluna_traj +  bage"
     cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao: ", equacao_perc, "\n")
     modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + bage)}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
     modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + bage)
     modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj + bage)
     modelo_perc <- modelo_perc_random
     cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
    } else if (tag_perc == "P95" && tag_transtorno == "cross") {
       equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit), dispformula = ~ coluna_traj)"
       cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao: ", equacao_perc, "\n")
       modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj)}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
       modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj)
       modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ coluna_traj)
       modelo_perc <- modelo_perc_random
       cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
    } else if  (tag_perc %in% c("P90", "P95") && tag_transtorno == "adhd_inc_control"){
       equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = logit), dispformula = ~ site)"
       cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao final: ", equacao_perc, "\n")
       modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ site)}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
       modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ site)
       modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"), dispformula = ~ site)
       modelo_perc <- modelo_perc_no_random
       perc_no_random_sem_disp <-  glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"))
       lrt_perc_fixo <- anova(modelo_perc, perc_no_random_sem_disp)
       cat("\n LRT FIXO \n"); print(lrt_perc_fixo)
       tabela_global_lrt <- rbind(tabela_global_lrt, data.frame(
         Disorder = tag_transtorno, Percentage_threshold = tag_perc, Dependent_Variable = "Perc FIXO", 
         LRT_p = round(lrt_perc_fixo$`Pr(>Chisq)`[2], 4), 
         Decision = ifelse(lrt_perc_fixo$`Pr(>Chisq)`[2] < 0.05, "Disp", "No disp")
       ))
       cat("\n--- 1.2. DHARMa Residual Diagnostics ---\n")
       residuos_nodisp <- simulateResiduals(fittedModel = perc_no_random_sem_disp); tabela_diag_nodisp <- criar_tabela_diagnostico_dharma_qq(residuos_nodisp)
       png(paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_perc_FIXO_NO_DISP_", tag_perc, "_",  tag_transtorno, ".png"), width = 3000, height = 1200, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 0.8); plotQQunif(residuos_nodisp, main = "A) QQ Plot - PERC FIXO SEM DISP", testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_nodisp, main = ""); title("B) Residuals vs. Predicted – Perc (Fixo sem disp)", line = 2); dev.off()
       cat("AIC"); print(AIC(perc_no_random_sem_disp))
       cat("BIC"); print(BIC(perc_no_random_sem_disp))
       } else {
       equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = logit))"
       cat("--- 3.1. Fitting the Beta GLMM for percentage ---\n", "equacao: ", equacao_perc, "\n")
       modelo_perc_random <- tryCatch({glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit"))}, warning = function(w) {cat("WARNING:\n", conditionMessage(w), "\n"); return(NULL)}, error = function(e) {cat("ERROR:\n", conditionMessage(e), "\n"); return(NULL)})
       modelo_perc_vif <- glmmTMB(percentage_prop ~ wave + coluna_traj + sex + site + bage + (1 | subjectid), data = nanosight_intersect, family = beta_family(link = "logit")) 
       modelo_perc_no_random <- glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = "logit"))
       if((tag_transtorno == "adhd_inc_control_remit" && tag_perc %in% c("P90", "P95")) || (tag_perc != "P90" && tag_transtorno == "anx_inc_control") || (tag_perc != "P95" && tag_transtorno == "anx_inc_control_remit")) {
       modelo_perc <- modelo_perc_no_random
       equacao_perc <- "glmmTMB(percentage_prop ~ wave * coluna_traj + sex + site + bage, data = nanosight_intersect, family = beta_family(link = logit))"
       cat("Beta GLMM for percentage ---\n", "equacao final: ", equacao_perc, "\n")
       } else if ((tag_transtorno == "anx_inc_control" && tag_perc == "P90") || (tag_transtorno == "anx_inc_control_remit" && tag_perc == "P95")) { #dharma showed heterocedasticity, so going back to the random
       modelo_perc <- modelo_perc_random
       cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
       icc_perc <- icc(modelo_perc)
       cat("\nICC:\n"); print(icc_perc)
       } else {
         modelo_perc <- modelo_perc_random
         cat("\nVarrCorr:\n"); print(VarCorr(modelo_perc))
         icc_perc <- icc(modelo_perc)
         cat("\nICC:\n"); print(icc_perc)
         
       }
      }
    resultados_finais_aic_bic <- rbind(resultados_finais_aic_bic, data.frame(Variable = "Percentage", Model = "Beta", AIC =AIC(modelo_perc), BIC = BIC(modelo_perc)))
    cat("\n--- 3.2. ANOVA Table (Type III) for percentage ---\n")
    anova_perc_type3 <- Anova(modelo_perc, type = "III"); df_anova_perc3 <- as.data.frame(anova_perc_type3); df_anova_perc3$Significance <- ifelse(df_anova_perc3$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc3); write.csv(df_anova_perc3, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_anova_percentage_", tag_perc, "_", tag_transtorno, ".csv"))
    cat("\n--- 3.3. Conditional Pairwise Comparisons (percentage) ---\n")
    p_interacao_perc <- df_anova_perc3["wave:coluna_traj", "Pr(>Chisq)"]
    if (p_interacao_perc < 0.05) {
      cat("Interaction is significant. Analyzing both perspectives.\n")
      emm_perc1 <- emmeans(modelo_perc, specs = pairwise ~ coluna_traj | wave, adjust = "tukey"); df_emm_perc1 <- as.data.frame(emm_perc1$contrasts); df_emm_perc1$Significance <- ifelse(df_emm_perc1$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc1); write.csv(df_emm_perc1, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_perc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
     emm_perc2 <- emmeans(modelo_perc, specs = pairwise ~ wave | coluna_traj, adjust = "tukey"); df_emm_perc2 <- as.data.frame(emm_perc2$contrasts); df_emm_perc2$Significance <- ifelse(df_emm_perc2$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc2); write.csv(df_emm_perc2, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_perc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    } else {
    cat("Interaction is not significant. Analyzing main effects with a Type II ANOVA.\n")
    anova_perc_type2 <- Anova(modelo_perc, type = "II"); df_anova_perc2 <- as.data.frame(anova_perc_type2); df_anova_perc2$Significance <- ifelse(df_anova_perc2$`Pr(>Chisq)` < alfa, "Significant", "Not Significant"); print(df_anova_perc2)
    if (df_anova_perc2["wave", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'wave' is significant. Performing contrasts.\n\n"); emm_perc_wave <- emmeans(modelo_perc, specs = pairwise ~ wave, adjust = "tukey"); df_emm_perc_wave <- as.data.frame(emm_perc_wave$contrasts); df_emm_perc_wave$Significance <- ifelse(df_emm_perc_wave$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_wave); write.csv(df_emm_perc_wave, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_wave_percentage_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    if (df_anova_perc2["coluna_traj", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'coluna_traj' is significant. Performing contrasts.\n\n");emm_perc <- emmeans(modelo_perc, specs = pairwise ~ coluna_traj, adjust = "tukey"); df_emm_perc <- as.data.frame(emm_perc$contrasts); df_emm_perc$Significance <- ifelse(df_emm_perc$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc); write.csv(df_emm_perc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_coluna_traj_percentage_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    if (anova_perc_type2["sex", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'sex' is significant. Performing contrasts.\n\n"); emm_perc_sex <- emmeans(modelo_perc, specs = pairwise ~ sex, adjust = "tukey"); df_emm_perc_sex <- as.data.frame(emm_perc_sex$contrasts); df_emm_perc_sex$Significance <- ifelse(df_emm_perc_sex$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_sex); write.csv(df_emm_perc_sex, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_sex_perc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    if (anova_perc_type2["site", "Pr(>Chisq)"] < 0.05) {cat("\n Main effect of 'site' is significant. Performing contrasts.\n\n"); emm_perc_site <- emmeans(modelo_perc, specs = pairwise ~ site, adjust = "tukey"); df_emm_perc_site <- as.data.frame(emm_perc_site$contrasts); df_emm_perc_site$Significance <- ifelse(df_emm_perc_site$p.value < alfa, "Significant", "Not Significant"); print(df_emm_perc_site); write.csv(df_emm_perc_site, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_contrasts_site_perc_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE) }
    }
    cat("\n--- 3.4. Model Diagnostics (percentage) ---\n")
    cat("Assessing VIF on a model without the interaction term to prevent artificial inflation:\n")
    vif_df_perc <- as.data.frame(check_collinearity(modelo_perc_vif)); vif_df_perc$Interpretation <- ifelse(vif_df_perc$VIF >= 10, "High", ifelse(vif_df_perc$VIF >= 5, "Moderate", "Low"));
    cat("VIF:\n"); print(vif_df_perc); write.csv(vif_df_perc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_vif_percentage_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    #Random ou nao
    lrt_perc <- anova(modelo_perc_random, modelo_perc_no_random)
    cat("\n LRT \n"); print(lrt_perc)
    tabela_global_lrt <- rbind(tabela_global_lrt, data.frame(
      Disorder = tag_transtorno, Percentage_threshold = tag_perc, Dependent_Variable = "Percentage", 
      LRT_p = round(lrt_perc$`Pr(>Chisq)`[2], 4), 
      Decision = ifelse(lrt_perc$`Pr(>Chisq)`[2] < 0.05, "Random factor", "No random factor")
    ))
    #Dharma
    cat("\n--- 3.5. DHARMa Residual Diagnostics (percentage) ---\n")
    residuos_perc <- simulateResiduals(fittedModel = modelo_perc); tabela_diag_perc <- criar_tabela_diagnostico_dharma_qq(residuos_perc)
    png(paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_dharma_percentage_", tag_perc, "_",  tag_transtorno, ".png"), width = 3000, height = 1200, res = 300); par(mfrow = c(1, 2), mar = c(4, 4, 3, 1) + 0.1, oma = c(0, 0, 1, 0), cex = 0.8); plotQQunif(residuos_perc, main = paste0("A) QQ Plot - percentage"), testUniformity = FALSE, testDispersion = FALSE, testOutliers = FALSE); plotResiduals(residuos_perc, main = ""); title("B) Residuals vs. Predicted - percentage", line = 2); dev.off()
    df_model <- model.frame(modelo_perc)
    vars <- c("sex", "site", "wave", "bage", coluna_traj)
    
    for(v in vars){
      # Define filename for each plot
      filename <-paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/plot_diagnostics_perc_", v, "_", tag_perc, "_",  tag_transtorno, ".png")
      
      # Open PNG device
      png(filename)
      
      # Plot residuals
      plotResiduals(residuos_perc, df_model[[v]], main = paste("Residuals vs", v))
      
      # Close the device to save the file
      dev.off()
    }
    print(tabela_diag_perc); write.csv(tabela_diag_perc, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_diagnostics_dharma_percentage_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    #AIC/BIC
    cat("\n\n=========================================================\n"); cat("  SUMMARY TABLE:AIC/BIC OF FINAL MODELS\n"); cat("=========================================================\n\n")
    print(resultados_finais_aic_bic)
    write.csv2(resultados_finais_aic_bic, paste0("analysis_results/", tag_perc, "_", tag_transtorno, "/table_summary_aic_bic_final_", tag_perc, "_", tag_transtorno, ".csv"), row.names = FALSE)
    #LRT
    write.csv2(tabela_global_lrt, "analysis_results/tabela_global_lrt.csv", row.names = F, quote = F)
    }, finally = {
    # LOG
    cat("\n\n### FIM DO LOG DE ANÁLISE ###\n")
    cat("Todos os resultados foram salvos no diretório 'analysis_results'.\n")
    cat("O log completo desta sessão foi salvo em:", log_file_name, "\n")
    sink(type = "message")
    while(sink.number() > 0) sink()
    if(exists("con")) close(con)
  }
  )}
