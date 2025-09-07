### Contagem de individuos com cada transtorno em cada time point
## dcany (any disorder - DSM IV)
#w1
dcany_w1 <- subset(nanosight_w1, grepl("T", dcany))
num_dcany_w1 <- nrow(dcany_w1)
num_dcany_w1
#w2
dcany_w2 <- subset(nanosight_w2, grepl("T", dcany))
num_dcany_w2 <- nrow(dcany_w2)
num_dcany_w2
## dcmadep (major depression)
#w1
dcmadep_w1 <- subset(nanosight_w1, grepl("T", dcmadep))
num_dcmadep_w1 <- nrow(dcmadep_w1)
num_dcmadep_w1
#w2
dcmadep_w2 <- subset(nanosight_w2, grepl("T", dcmadep))
num_dcmadep_w2 <- nrow(dcmadep_w2)
num_dcmadep_w2
## dcanyanx (any anxiety)
#w1
dcanyanx_w1 <- subset(nanosight_w1, grepl("T", dcanyanx))
num_dcanyanx_w1 <- nrow(dcanyanx_w1)
num_dcanyanx_w1
#w2
dcanyanx_w2 <- subset(nanosight_w2, grepl("T", dcanyanx))
num_dcanyanx_w2 <- nrow(dcanyanx_w2)
num_dcanyanx_w2
## dcgena (generalised anxiety)
#w1
dcgena_w1 <- subset(nanosight_w1, grepl("T", dcgena))
num_dcgena_w1 <- nrow(dcgena_w1)
num_dcgena_w1
#w2
dcgena_w2 <- subset(nanosight_w2, grepl("T", dcgena))
num_dcgena_w2 <- nrow(dcgena_w2)
num_dcgena_w2
## dcanyhk (ADHD)
#w1
dcanyhk_w1 <- subset(nanosight_w1, grepl("T", dcanyhk))
num_dcanyhk_w1 <- nrow(dcanyhk_w1)
num_dcanyhk_w1
#w2
dcanyhk_w2 <- subset(nanosight_w2, grepl("T", dcanyhk))
num_dcanyhk_w2 <- nrow(dcanyhk_w2)
num_dcanyhk_w2
## dcpsych (psychosis)
#w1
dcpsych_w1 <- subset(nanosight_w1, grepl("T", dcpsych))
num_dcpsych_w1 <- nrow(dcpsych_w1)
num_dcpsych_w1
#w2
dcpsych_w2 <- subset(nanosight_w2, grepl("T", dcpsych))
num_dcpsych_w2 <- nrow(dcpsych_w2)
num_dcpsych_w2
## dcmania (mania/bip)
#w1
dcmania_w1 <- subset(nanosight_w1, grepl("T", dcmania))
num_dcmania_w1 <- nrow(dcmania_w1)
num_dcmania_w1
#w2
dcmania_w2 <- subset(nanosight_w2, grepl("T", dcmania))
num_dcmania_w2 <- nrow(dcmania_w2)
num_dcmania_w2
## dcptsd 
#w1
dcptsd_w1 <- subset(nanosight_w1, grepl("T", dcptsd))
num_dcptsd_w1 <- nrow(dcptsd_w1)
num_dcptsd_w1
#w1
dcptsd_w2 <- subset(nanosight_w2, grepl("T", dcptsd))
num_dcptsd_w2 <- nrow(dcptsd_w2)
num_dcptsd_w2

# Pacotes necessários
library(ggplot2)
library(dplyr)

# Criando o data frame com os dados de contagem
df_transtornos <- data.frame(
  Disorder = rep(c("Any disorder", "Major depression", "Anxiety disorder", 
                     "Generalized anxiety", "Hyperactivity and Deficit Attention Disorder", "Psychosis", "Mania/Bipolar disorder", "Post Traumatic Stress Disorder"), each = 2),
  Timepoint = rep(c("t1", "t2"), times = 8),
  N = c(num_dcany_w1, num_dcany_w2,
        num_dcmadep_w1, num_dcmadep_w2,
        num_dcanyanx_w1, num_dcanyanx_w2,
        num_dcgena_w1, num_dcgena_w2,
        num_dcanyhk_w1, num_dcanyhk_w2,
        num_dcpsych_w1, num_dcpsych_w2,
        num_dcmania_w1, num_dcmania_w2,
        num_dcptsd_w1, num_dcptsd_w2)
)

# Organiza fatores para ordem correta
df_transtornos$Disorder <- factor(df_transtornos$Disorder,
                                    levels = c("Any disorder", "Major depression", "Anxiety disorder", 
                     "Generalized anxiety", "Hyperactivity and Deficit Attention Disorder", "Psychosis", "Mania/Bipolar disorder", "Post Traumatic Stress Disorder"))
# Gráfico de barras agrupadas
ggplot(df_transtornos, aes(x = Disorder, y = N, fill = Timepoint)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6) +
  geom_text(aes(label = N), position = position_dodge(width = 0.8), vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("t1" = "#FFB6C1", "t2" = "#81D8D0"), labels = c("t 1", "t 2")) +
  labs(title = "Number of individuals with each disorder by time point",
       x = "Mental condition",
       y = "Number of individuals",
       fill = "Timepoint") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )
# Salvar o gráfico
ggsave("contagem_transtornos_por_timepoint.png", width = 10, height = 8, dpi = 300)

#### Incident #######################
### Contagem de indivíduos com cada transtorno
# Lista de colunas com os transtornos
colunas_transtornos_incident <- c(
  "dcany", "dcmadep", "dcanyanx", "dcgena", "dcanyhk", 
  "dcpsych", "dcmania", "dcptsd"
)
# Nomes para o gráfico (em ordem)
nomes_transtornos_incident <- c(
  "Any disorder", "Major depression", "Anxiety disorder", 
  "Generalized anxiety", "Hyperactivity and Deficit Attention Disorder", 
  "Psychosis", "Mania/Bipolar disorder", "Post Traumatic Stress Disorder"
)
# Conta o número de ocorrências de "T" em cada coluna usando sapply
contagens_incident <- sapply(colunas_transtornos_incident, function(colunas_transtornos_incident) {
  sum(nanosight_intersect_Incident[[colunas_transtornos_incident]] == "T")
})
# Criando o data frame de forma mais eficiente
df_transtornos_incident <- data.frame(
  Disorder = nomes_transtornos_incident,
  N = as.numeric(contagens_incident)
)
# Organiza os fatores para a ordem correta no gráfico
df_transtornos_incident$Disorder <- factor(df_transtornos_incident$Disorder,
                                           levels = nomes_transtornos_incident)
# Gráfico de barras
ggplot(df_transtornos_incident, aes(x = Disorder, y = N)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.6, fill = "#FFB6C1") +
  geom_text(aes(label = N), position = position_dodge(width = 0.8), vjust = -0.5, size = 4.5) +
  labs(title = "Number of individuals with each disorder in the Incident Trajectory",
       x = "Mental condition",
       y = "Number of individuals") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 14),
    legend.position = "none" # Remove a legenda de cor se houver
  )
# Salvar o gráfico
ggsave("contagem_transtornos_incident.png", width = 10, height = 8, dpi = 300)
