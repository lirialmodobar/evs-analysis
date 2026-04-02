#Importar bibliotecas e scripts necessarios
library(tidyr)
library(dplyr)
library(readr)
library(broom)

#Definir pasta de trabalho 
setwd("C:/Users/lirie/belle/dados_raw")
save_dir <- ("C:/Users/lirie/belle/")

#Trazer tabelas de interesse para o R
isolar_exp_sum <- list.files(pattern="*ExperimentSummary.csv") 
importar_exp_sum <- lapply(isolar_exp_sum, read.csv, header = FALSE, col.names = c(paste0("V", 1:6))) #importa todas as expsum corretamente ao nomear as colunas

#Funcoes

##Para ter o dado do ID
extract_standardized_id <- function(ids) {
  standardized_ids <- character()
  
  for (i in seq_along(ids)) {
    id <- ids[i]
    id_number <- sub('.*?(\\d{5}).*', '\\1', id, perl = TRUE)
    
    if (grepl("(?i)\\d{5}$", id)) {
      if (grepl("(?i)^[wW][12]", id)) {
        id_number <- paste0(id_number, "_", tolower(sub("^[wW]([12]).*", "\\1", id, perl = TRUE)))
        standardized_ids <- c(standardized_ids, id_number)
      } else {
        warning(paste("Sample at position", i, "(", id, ") does not contain 'w1' or 'w2' in its name. It has not been outputted."))
      }
    } else if (grepl("(?i)^[wW][12]", id)) {
      id_number <- paste0(id_number, "_", tolower(sub("^[wW]([12]).*", "\\1", id, perl = TRUE)))
      standardized_ids <- c(standardized_ids, id_number)
    } else if (grepl("(?i)([wW][12])", id)) {
      id_number <- paste0(id_number, "_", tolower(sub(".*?([wW][12]).*", "\\1", id, perl = TRUE)))
      standardized_ids <- c(standardized_ids, id_number)
    } else {
      warning(paste("Sample at position", i, "(", id, ") does not match the expected format. It has not been outputted."))
    }
  }
  
  return(standardized_ids)
}

index_lapply_sapply <- function(...) {
  ...[]
} 

isolar_bins <- function(df) {
  inicio <- which(df[, 1] == "Bin centre (nm)")[1] + 1
  fim    <- which(df[, 1] == "Percentile")[1] - 1
  tabela_dados <- df[inicio:fim, ]
  tabela_dados$V1 <- as.numeric(tabela_dados$V1)
  tabela_dados$V5 <- as.numeric(tabela_dados$V5)
  colnames(tabela_dados) <- c("diametros", "V2", "V3", "V4", "concentracoes", "V6")
  return(tabela_dados)
}

processar_exp <- function(df, corte) {
  corte <- corte
  diametros <- df$diametros
  concentracoes <- df$concentracoes
  pequenas <- sum(concentracoes[diametros <= corte], na.rm = TRUE)
  grandes  <- sum(concentracoes[diametros > corte], na.rm = TRUE)
  return(c(pequenas = pequenas, grandes = grandes))
}
#Extrair dados de interesse de ExpSum
##Para ter o dado do ID
isolar_exp_sum2 <- gsub (" W1", "_w1", isolar_exp_sum)
isolar_exp_sum2 <- gsub (" W2", "_w2", isolar_exp_sum2)
id_sample <- separate(as.data.frame(isolar_exp_sum2), 1, into=c("ids", "resto"), sep=" ") [,1]
id_sample <- extract_standardized_id(id_sample) #padroniza as ids para como o inpd usa
duplicados_indices <- duplicated(id_sample, fromLast = TRUE)
id_sample <- id_sample[!duplicados_indices]
importar_exp_sum <- importar_exp_sum[!duplicados_indices]
##Para ter o dado do quanto diluiu (assumindo que está no campo diluent, e não dilution factor, como a Jess fez)
diluicao <- as.data.frame(sapply(importar_exp_sum, index_lapply_sapply, 11, 2)) #se estiver no campo diluent (dados Jessica)
diluicao[diluicao == "" | is.na(diluicao)] <- deparse(quote(1/100)) #se nao estiver em campo nenhum e n tiver precisado diluir mais (dados padrao Belle)
diluicao[diluicao == "PBS"] <- "1/100"
diluicao[diluicao == "10X diluido" | diluicao == "10x diluida" | diluicao == "10x" | diluicao == "900ul PBS + 100ul amostra"] <- "1/1000"
##Para ter a concentração da amostra antes de diluir (concentração real)
concentracao_average <- as.numeric(sapply(importar_exp_sum, index_lapply_sapply, 44, 5)) #obter média da concentração
separar_fracao_diluicao <- separate(diluicao, 1, into=c("numerador_diluicao", "denominador_diluicao"), sep="/") #gera uma tabela com o numerador em uma coluna e denominador em outra para que seja possível ter um dado do tipo numeric a partir da fração de diluição anotada (usando só as.numeric, retorna NA, pq a barra confunde o R)
separar_fracao_diluicao <- separar_fracao_diluicao %>%
  mutate(numerador_diluicao = ifelse(row_number() == 175, "1", numerador_diluicao)) %>%
  mutate(denominador_diluicao = ifelse(row_number() == 175, "10", denominador_diluicao)) #uma amostra que estava desconfigurada
numerador <- as.numeric(separar_fracao_diluicao[,1]) #isolar o numerador da tabela gerada para poder transformar em numeric
denominador <- as.numeric(separar_fracao_diluicao[,2]) #isolar o denominador da tabela gerada para poder transformar em numeric
concentracao_real <- concentracao_average/(numerador/denominador) #obtém de fato a concentração real
##Para obter o restante dos dados
tamanho_mean_mode <- sapply(importar_exp_sum, function(df) {
  indice_tamanho <- which(df[, 1] == "Mean")[1] 
  indice_mode <- which(df[, 1] == "Mode")[1] 
  tamanho_mean_average <- as.numeric(df[indice_tamanho,5])
  tamanho_mode_average <- as.numeric(df[indice_mode,5])
  df <- data.frame(tamanho_mean_average, tamanho_mode_average)
  return(df)
})
tamanho_mean_mode <- t(tamanho_mean_mode)
##Porcentagem de vesículas pequenas e grandes 
###Adequar dados (info em linhas diferentes dependendo da tabela)
bins_isolados <- lapply(importar_exp_sum, isolar_bins)
evs_bin_size <- sapply(bins_isolados, processar_exp, 128.5)
concentration_EV_pequenas <- evs_bin_size["pequenas", ]
concentration_EV_grandes  <- evs_bin_size["grandes", ]
EV_pequenas_porcentagem <- (concentration_EV_pequenas/(concentration_EV_grandes + concentration_EV_pequenas)) * 100
### Calculando os valores para cada amostra
valores_percentis <- sapply(bins_isolados, function(df) {
  diametros <- df$diametros
  concentracoes <- df$concentracoes
  
  ### Calcula o acumulado e identifica os tamanhos
  total_conc <- sum(concentracoes)
  acumulado  <- cumsum(concentracoes)
  
  ### Acha o tamanho onde bate 90% e 95%
  p90_amostra <- diametros[which(acumulado >= total_conc * 0.90)[1]]
  p95_amostra <- diametros[which(acumulado >= total_conc * 0.95)[1]]
  
  return(c(p90 = p90_amostra, p95 = p95_amostra))
})
corte_p90 <- mean(valores_percentis["p90", ], na.rm = TRUE)
print(corte_p90)
corte_p95 <- mean(valores_percentis["p95", ], na.rm = TRUE)
print(corte_p95)
p90 <- sapply(bins_isolados, processar_exp, corte_p90)
concentration_EV_pequenas <- p90["pequenas", ]
concentration_EV_grandes  <- p90["grandes", ]
p90_porcentagem <- (concentration_EV_grandes/(concentration_EV_pequenas + concentration_EV_grandes)) * 100
p95 <- sapply(bins_isolados, processar_exp, corte_p95)
concentration_EV_pequenas <- p95["pequenas", ]
concentration_EV_grandes  <- p95["grandes", ]
p95_porcentagem <- (concentration_EV_grandes/(concentration_EV_pequenas + concentration_EV_grandes)) * 100
##Batch Nanosight
batch_nan <- sapply(importar_exp_sum, index_lapply_sapply, 8, 2)
batch_nan <- as.data.frame(batch_nan)
batch_nan <- sapply(strsplit(batch_nan$batch_nan, " "), tail, 1)
batch_nan <- as.data.frame(batch_nan)
datas_unicas <- duplicated(batch_nan$batch_nan)
datas_unicas <- subset(batch_nan, !datas_unicas) #para ver quantos grupos de datas (batch)
datas_unicas <- datas_unicas[,1]
numeros_batch <- seq(1, 22)
replacement <- setNames(numeros_batch, datas_unicas)
batch_nan$batch_nan <- replacement[batch_nan$batch_nan]
batch_nan <- batch_nan %>%
  mutate(batch_nan = case_when(
    batch_nan %in% c(1, 13) ~ "B",  # Agrupar 1 e 13 para B
    batch_nan == 2 ~ "J",           # Agrupar 2 para J
    batch_nan %in% c(4, 5, 9, 14, 20, 21) ~ "I1",
    batch_nan %in% c(3, 6, 7, 8, 10, 11, 12, 15, 16, 17, 18, 19, 22) ~ "I2",
    TRUE ~ as.character(batch_nan)  # Manter os outros números como estão
  )) #Jess, Belle e ICESP
#Juntando todos os dados obtidos do nanosight em uma tabela
tabela_tudo_nanosight <- data.frame(id_sample, concentracao_average, diluicao, concentracao_real, tamanho_mean_mode, EV_pequenas_porcentagem, p90_porcentagem, p95_porcentagem, batch_nan)
colnames(tabela_tudo_nanosight)[3] <- "diluicao"
duplicatas <- duplicated(tabela_tudo_nanosight[, c("id_sample")], fromLast = TRUE) #a 2a leitura é a correta
nanosight_sem_duplicatas <- subset(tabela_tudo_nanosight, !duplicatas)
#Juntar dados do nanosight com informações da amostra
sample_information <- read.csv("sampleinformation_liriel320_atualizada.csv") #chamando/lendo/importando infos da amostra
sample_information$id_sample <- paste(sample_information$subjectid, sample_information$wave, sep = "_")
nanosight_plus_sampleinfo <- inner_join(nanosight_sem_duplicatas, sample_information, by="id_sample")
nanosight_plus_sampleinfo$Grupo <- paste(nanosight_plus_sampleinfo$Trajetoria, nanosight_plus_sampleinfo$wave, sep = "_")
nanosight_plus_sampleinfo <- unnest(
  nanosight_plus_sampleinfo,
  cols = c("tamanho_mean_average", "tamanho_mode_average")
)

##Para analise de Levene as colunas de transtornos não podem ser numeros, portanto (F=False T=True)
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(
    dcany = ifelse(dcany == 0, "F", "T"),
    dcmadep = ifelse(dcmadep == 0, "F", "T"),
    dcanyanx = ifelse(dcanyanx == 0, "F", "T"),
    dcgena = ifelse(dcgena == 0, "F", "T"),
    dcanyhk = ifelse(dcanyhk == 0, "F", "T"),
    dcpsych = ifelse(dcpsych == 0, "F", "T"),
    dcmania = ifelse(dcmania == 0, "F", "T"),
    dcptsd = ifelse(dcptsd == 0, "F", "T")
  )
##Nomeando as trajetórias para deixar menos confuso
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(Trajetoria = case_when(
    Trajetoria == "A" ~ "Control",
    Trajetoria == "B" ~ "Incident",
    Trajetoria == "C" ~ "Remitted",
    Trajetoria == "D" ~ "Persistent"
  ))
##Alterando wave para time point
nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  mutate(wave = case_when(
    wave == "w1" ~ "t1",
    wave == "w2" ~ "t2"
  ))

nanosight_plus_sampleinfo <- nanosight_plus_sampleinfo %>%
  group_by(subjectid) %>%
  mutate(
    dep_t1 = first(dcmadep[wave == "t1"], default = NA),
    dep_t2 = first(dcmadep[wave == "t2"], default = NA),
    grupo_analise_dep = case_when(
      is.na(dep_t1) | is.na(dep_t2) ~ "None",
      Trajetoria == "Control" ~ "Control",
      dep_t1 == "F" & dep_t2 == "T" ~ "Incident",
      dep_t1 == "T" & dep_t2 == "T" ~ "Persistent",
      dep_t1 == "T" & dep_t2 == "F" ~ "Remitted",
      TRUE ~ NA_character_),
    gena_t1 = first(dcgena[wave == "t1"], default = NA),
    gena_t2 = first(dcgena[wave == "t2"], default = NA),
    grupo_analise_gena = case_when(
      is.na(gena_t1) | is.na(gena_t2) ~ "None",
      Trajetoria == "Control" ~ "Control",
      gena_t1 == "F" & gena_t2 == "T" ~ "Incident",
      gena_t1 == "T" & gena_t2 == "T" ~ "Persistent",
      gena_t1 == "T" & gena_t2 == "F" ~ "Remitted",
      TRUE ~ NA_character_),
    anyhk_t1 = first(dcanyhk[wave == "t1"], default = NA),
    anyhk_t2 = first(dcanyhk[wave == "t2"], default = NA),
    grupo_analise_anyhk = case_when(
      is.na(anyhk_t1) | is.na(anyhk_t2) ~ "None",
      Trajetoria == "Control" ~ "Control",
      anyhk_t1 == "F" & anyhk_t2 == "T" ~ "Incident",
      anyhk_t1 == "T" & anyhk_t2 == "T" ~ "Persistent",
      anyhk_t1 == "T" & anyhk_t2 == "F" ~ "Remitted",
      TRUE ~ NA_character_)
    ) %>% ungroup()
  
  
write.csv(nanosight_plus_sampleinfo, paste0(save_dir, "nanosight_plus_sampleinfo_all.csv"), row.names = FALSE)

transtorno_traj <- list(
  list(coluna_trajetoria = "Trajetoria", transtorno = "cross"),
  list(coluna_trajetoria = "grupo_analise_dep", transtorno = "dep"),
  list(coluna_trajetoria = "grupo_analise_gena", transtorno = "anx_inc_control"), #n baixo, so incidente vs controle
  list(coluna_trajetoria = "grupo_analise_gena", transtorno = "anx_inc_control_remit"), #n baixo, so incidente vs controle
  list(coluna_trajetoria = "grupo_analise_anyhk", transtorno = "adhd_inc_control"), #n baixo, so incidente vs controle
  list(coluna_trajetoria = "grupo_analise_anyhk", transtorno = "adhd_inc_control_remit") #n baixo, so incidente vs controle
)

tabela_decisao_modelagem <- data.frame()

for (transtornos in transtorno_traj) {

transtorno <- transtornos$transtorno
coluna_trajetoria <- transtornos$coluna_trajetoria

nanosight_plus_sampleinfo$coluna_traj <- nanosight_plus_sampleinfo[[coluna_trajetoria]]

nanosight_plus_sampleinfo_transt <- subset(nanosight_plus_sampleinfo, !is.na(coluna_traj)) 
nanosight_plus_sampleinfo_transt <- subset(nanosight_plus_sampleinfo_transt, coluna_traj != "None") 

if(transtorno == "adhd_inc_control" || transtorno == "anx_inc_control") {
nanosight_plus_sampleinfo_transt <- subset(nanosight_plus_sampleinfo_transt, coluna_traj %in% c("Control","Incident")) 
  
} else if (transtorno == "adhd_inc_control_remit" || transtorno == "anx_inc_control_remit") {
  nanosight_plus_sampleinfo_transt <- subset(nanosight_plus_sampleinfo_transt, coluna_traj != "Persistent") 
}


#Retirando outliers
nanosight_plus_sampleinfo_transt$zscore_mode <- scale(nanosight_plus_sampleinfo_transt$tamanho_mode_average)
nanosight_plus_sampleinfo_transt$zscore_porcentagem <- scale(nanosight_plus_sampleinfo_transt$EV_pequenas_porcentagem)
nanosight_plus_sampleinfo_transt$zscore_concentracao <- scale(nanosight_plus_sampleinfo_transt$concentracao_real)
nanosight_plus_sampleinfo_transt$zscore_p90 <- scale(nanosight_plus_sampleinfo_transt$p90_porcentagem)
nanosight_plus_sampleinfo_transt$zscore_p95 <- scale(nanosight_plus_sampleinfo_transt$p95_porcentagem)


nanosight_plus_sampleinfo_sem_outliers_mode <- nanosight_plus_sampleinfo_transt[nanosight_plus_sampleinfo_transt$zscore_mode >= -3.0 & nanosight_plus_sampleinfo_transt$zscore_mode <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_porcentagem <- nanosight_plus_sampleinfo_transt[nanosight_plus_sampleinfo_transt$zscore_porcentagem >= -3.0 & nanosight_plus_sampleinfo_transt$zscore_porcentagem <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_concentracao <- nanosight_plus_sampleinfo_transt[nanosight_plus_sampleinfo_transt$zscore_concentracao >= -3.0 & nanosight_plus_sampleinfo_transt$zscore_concentracao <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_p90 <- nanosight_plus_sampleinfo_transt[nanosight_plus_sampleinfo_transt$zscore_p90 >= -3.0 & nanosight_plus_sampleinfo_transt$zscore_p90 <= 3.0, ]
nanosight_plus_sampleinfo_sem_outliers_p95 <- nanosight_plus_sampleinfo_transt[nanosight_plus_sampleinfo_transt$zscore_p95 >= -3.0 & nanosight_plus_sampleinfo_transt$zscore_p95 <= 3.0, ]


nanosight_plus_sampleinfo_sem_outliers_mode_concentracao <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mode, nanosight_plus_sampleinfo_sem_outliers_concentracao, by = "id_sample")

nanosight_intersect_ev_pequena <- semi_join(nanosight_plus_sampleinfo_sem_outliers_mode_concentracao, nanosight_plus_sampleinfo_sem_outliers_porcentagem, by = ("id_sample"))
nanosight_intersect_p90 <-  semi_join(nanosight_plus_sampleinfo_sem_outliers_mode_concentracao, nanosight_plus_sampleinfo_sem_outliers_p90, by = ("id_sample"))
nanosight_intersect_p95 <-  semi_join(nanosight_plus_sampleinfo_sem_outliers_mode_concentracao, nanosight_plus_sampleinfo_sem_outliers_p95, by = ("id_sample"))

write.csv(nanosight_intersect_ev_pequena, paste0(save_dir, "nanosight_intersect_ev_pequena_", transtorno, ".csv"), row.names = FALSE, quote = FALSE)
write.csv(nanosight_intersect_p90, paste0(save_dir, "nanosight_intersect_p90_", transtorno, ".csv"), row.names = FALSE, quote = FALSE)
write.csv(nanosight_intersect_p95, paste0(save_dir, "nanosight_intersect_p95_", transtorno, ".csv"), row.names = FALSE, quote = FALSE)

# --- Definir quais modelos serao rodados e quais sao redundantes ---
# 1. Coleta as IDs de cada grupo
ids_90  <- nanosight_intersect_p90$id_sample
ids_95  <- nanosight_intersect_p95$id_sample
ids_peq <- nanosight_intersect_ev_pequena$id_sample

# 2. Une todas as IDs únicas que aparecem em qualquer um dos 3 arquivos
todos_ids_do_transtorno <- unique(c(ids_90, ids_95, ids_peq))

# 3. Identifica quem NÃO está presente nos três simultaneamente
ids_conflito <- todos_ids_do_transtorno[!(todos_ids_do_transtorno %in% ids_90 & 
                                            todos_ids_do_transtorno %in% ids_95 & 
                                            todos_ids_do_transtorno %in% ids_peq)]

# 4. Cria uma descrição detalhada para cada indivíduo problemático
if (length(ids_conflito) > 0) {
  detalhes <- sapply(ids_conflito, function(id) {
    faltas <- c()
    if (!(id %in% ids_90))  faltas <- c(faltas, "P90")
    if (!(id %in% ids_95))  faltas <- c(faltas, "P95")
    if (!(id %in% ids_peq)) faltas <- c(faltas, "Pequena")
    paste0(id, " (Falta em: ", paste(faltas, collapse = ", "), ")")
  })
  quem_e_diferente <- paste(detalhes, collapse = " | ")
} else {
  quem_e_diferente <- "Nenhum (Grupos Idênticos)"
}

# --- LÓGICA DE DECISÃO (IGUAL ANTERIOR) ---
comp_90_95  <- isTRUE(all.equal(ids_90, ids_95))
comp_90_peq <- isTRUE(all.equal(ids_90, ids_peq))
comp_95_peq <- isTRUE(all.equal(ids_95, ids_peq))

if (comp_90_95 & comp_90_peq) {
  decisao <- "ÚNICO"; arquivos_reportar <- "P90 (vale para todos)"
} else if (comp_90_95 & !comp_90_peq) {
  decisao <- "DUPLO"; arquivos_reportar <- "P90 (vale para P95); Pequena"
} else if (!comp_90_95 & comp_95_peq) {
  decisao <- "DUPLO"; arquivos_reportar <- "P90; P95 (vale para Pequena)"
} else if (!comp_90_95 & !comp_90_peq & !comp_95_peq) {
  decisao <- "TRIPLO"; arquivos_reportar <- "P90; P95; Pequena"
} else {
  decisao <- "DUPLO (P90=Peq)"; arquivos_reportar <- "P90 (vale para Pequena); P95"
}

# --- ARMAZENANDO NA TABELA ---
linha_decisao <- data.frame(
  Transtorno = transtorno,
  N_P90 = length(ids_90),
  N_P95 = length(ids_95),
  N_Pequena = length(ids_peq),
  Decisao_Relatorio_Size_Conc = decisao,
  Arquivos_para_Modelar = arquivos_reportar,
  Individuos_Diferentes = quem_e_diferente  # NOVA COLUNA
)
tabela_decisao_modelagem <- rbind(tabela_decisao_modelagem, linha_decisao)
write.csv(tabela_decisao_modelagem, "tabela_decisao_modelagem.csv", row.names = F, quote = F)
#Infos adicionais de numeros de amostras, transtornos, etc

##Outliers
nanosight_outliers_ev_pequena <- anti_join(nanosight_plus_sampleinfo, nanosight_intersect_ev_pequena, by = "id_sample")
nanosight_outliers_p95 <- anti_join(nanosight_plus_sampleinfo, nanosight_intersect_p90, by = "id_sample")
nanosight_outliers_p90 <- anti_join(nanosight_plus_sampleinfo, nanosight_intersect_p95, by = "id_sample")
print(transtorno)
print("outliers: pequena, p90, p95")
print(nanosight_outliers_ev_pequena)
print(nanosight_outliers_p90)
print(nanosight_outliers_p95)
##Amostras sem pares
linhas_sem_pares_intersect_ev_pequena <- !duplicated(nanosight_intersect_ev_pequena$subjectid) & !duplicated(nanosight_intersect_ev_pequena$subjectid, fromLast = TRUE)
linhas_sem_pares_intersect_p90 <- !duplicated(nanosight_intersect_p90$subjectid) & !duplicated(nanosight_intersect_p90$subjectid, fromLast = TRUE)
linhas_sem_pares_intersect_p95 <- !duplicated(nanosight_intersect_p95$subjectid) & !duplicated(nanosight_intersect_p95$subjectid, fromLast = TRUE)
nanosight_sem_pares_intersect_ev_pequena <- nanosight_intersect_ev_pequena[linhas_sem_pares_intersect_ev_pequena, ] 
nanosight_sem_pares_intersect_p90 <- nanosight_intersect_p90[linhas_sem_pares_intersect_p90, ] 
nanosight_sem_pares_intersect_p95 <- nanosight_intersect_p95[linhas_sem_pares_intersect_p95, ] 
print(transtorno)
print("sem pares: pequena, p90, p95")
print(nanosight_sem_pares_intersect_ev_pequena)
print(nanosight_sem_pares_intersect_p90)
print(nanosight_sem_pares_intersect_p95)

##Amostras com pares
nanosight_intersect_pares_ev_pequena <- anti_join (nanosight_intersect_ev_pequena, nanosight_sem_pares_intersect_ev_pequena, by = "id_sample")
nanosight_intersect_pares_p90 <- anti_join (nanosight_intersect_p90, nanosight_sem_pares_intersect_p90, by = "id_sample")
nanosight_intersect_pares_p95 <- anti_join (nanosight_intersect_p95, nanosight_sem_pares_intersect_p95, by = "id_sample")
print(transtorno)
print("pares: pequena, p90, p95")
print(nanosight_intersect_pares_ev_pequena)
print(nanosight_intersect_pares_p90)
print(nanosight_intersect_pares_p95)
}

#Infos adicionais gerais (com outliers) de numeros de amostras, transtornos, etc

##Tabela só com w1 e só w2
nanosight_w1 <- subset(nanosight_plus_sampleinfo, wave == "t1")
nanosight_w2 <- subset(nanosight_plus_sampleinfo, wave == "t2")

##Amostras sem análise no Nanosight ou sem correspondencia no sample info
sample_info_sem_nanosight <- anti_join(sample_information, nanosight_sem_duplicatas, by="id_sample")
sem_correspondencia <- anti_join(nanosight_sem_duplicatas, sample_information, by = "id_sample")


##Amostras sem pares
linhas_sem_pares <- !duplicated(nanosight_plus_sampleinfo$subjectid) & !duplicated(nanosight_plus_sampleinfo$subjectid, fromLast = TRUE)
nanosight_sem_pares <- nanosight_plus_sampleinfo[linhas_sem_pares, ] #perde 7

##Amostras com pares
nanosight_plus_sampleinfo_pares <- anti_join(nanosight_plus_sampleinfo, nanosight_sem_pares, by = "id_sample")


##Contagem de indivíduos com diferentes transtornos:
###w1
variaveis_w1 <- nanosight_w1 %>% 
  select(dcany, dcmadep, dcanyanx, dcgena, dcanyhk, dcpsych, dcptsd)
variaveis_w1_combinacao <- as.data.frame(rowSums(variaveis_w1 == "T"))
combinacoes_w1_frequencias <- variaveis_w1 %>%
  group_by(across(everything())) %>%
  summarise(Frequencia = n(), .groups = "drop")
###w2
variaveis_w2 <- nanosight_w2 %>% 
 select(dcany, dcmadep, dcanyanx, dcgena, dcanyhk, dcpsych, dcptsd)
variaveis_w2_combinacao <- as.data.frame(rowSums(variaveis_w2 == "T"))
combinacoes_w2_frequencias <- variaveis_w2 %>%
  group_by(across(everything())) %>%
  summarise(Frequencia = n(), .groups = "drop") 




