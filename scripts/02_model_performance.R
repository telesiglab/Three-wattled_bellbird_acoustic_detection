
here::i_am("scripts/02_model_performance.R")
library(rio)
library(dplyr)
library(tuneR)
library(patchwork)
library(tidyverse)
library(pROC)
library(PRROC)
library(here)

set.seed(2026)

#-------------------------------------------
# 1. Leer datos
#-------------------------------------------
pathtxt <- "data/raw/anotaciones"

datos <- data.frame()

files <- list.files(path = file.path(pathtxt), pattern = ".txt", full.names = TRUE)

sitios <- substr(basename(files), 1, 2)
as.data.frame(table(sitios))
length(sitios)

for (i in files) {
  pdf <- rio::import(i)
  datos <- rbind(datos, pdf)
}

names(datos)
table(datos$Especie)
table(datos$`Begin File`)
table(datos$Species[datos$Calidad %in% c('A','B')])

data <- datos %>%
  dplyr::rename(filename = 'Begin File', label = Species) 
#%>%
#  dplyr::select(filename, label, "File Offset (s)", Calidad, duration)

target <- "PRTR1"

data <- data %>%
  filter(label == target)


write.csv(data, file = "data/processed/anotaciones.csv", row.names = FALSE)





annotations <- rio::import(here::here("data", "raw", "raven_annotations.csv"))

model_predictions <- rio::import(here::here("data", "processed", "model_predictions.csv"))

#-------------------------------------------
# 2. Filtrar especie
#-------------------------------------------

target <- "PRTR1"

manual_labels <- annotations %>%
  filter(label == target)

#-------------------------------------------
# 3. Crear variable TRUE/FALSE por ventana
#-------------------------------------------

model_predictions$truth <- 0

for(i in 1:nrow(manual_labels)){
  
  f <- manual_labels$filename[i]
  offset <- manual_labels$`File Offset (s)`[i]
  
  idx <- which(
    model_predictions$filename == f &
      offset >= model_predictions$window_start &
      offset <= model_predictions$window_end
  )
  
  model_predictions$truth[idx] <- 1
}

#-------------------------------------------
# 4. Funcion para calcular metricas
#-------------------------------------------

metricas <- function(threshold){
  
  pred_class <- ifelse(model_predictions$logits >= threshold,1,0)
  
  TP <- sum(pred_class==1 & model_predictions$truth==1)
  FP <- sum(pred_class==1 & model_predictions$truth==0)
  FN <- sum(pred_class==0 & model_predictions$truth==1)
  TN <- sum(pred_class==0 & model_predictions$truth==0)
  
  precision <- TP/(TP+FP)
  recall <- TP/(TP+FN)
  specificity <- TN/(TN+FP)
  accuracy <- (TP+TN)/(TP+TN+FP+FN)
  
  F1 <- 2*(precision*recall)/(precision+recall)
  
  data.frame(threshold, TP, FP, FN, TN, precision, recall, specificity, accuracy, F1)
}

#-------------------------------------------
# 5. Evaluar muchos thresholds
#-------------------------------------------

thresholds <- sort(unique(model_predictions$logits))

results <- map_df(thresholds, metricas)

#-------------------------------------------
# 6. Threshold óptimo (max F1)
#-------------------------------------------

best <- results %>%
  filter(F1 == max(F1, na.rm=TRUE))

print(best)

#-------------------------------------------
# 7. Evaluar específicamente threshold = 0.01
#-------------------------------------------

res_001 <- metricas(0.01)

print(res_001)

#-------------------------------------------
# 8. Curvas ROC y Precision-Recall
#-------------------------------------------

library(patchwork)

roc_obj <- roc(model_predictions$truth, model_predictions$logits, quiet = TRUE)
pr_obj <- pr.curve(scores.class0 = model_predictions$logits[model_predictions$truth == 1],
                   scores.class1 = model_predictions$logits[model_predictions$truth == 0], curve = TRUE)

roc_df <- tibble(fpr = 1 - roc_obj$specificities, sensitivity = roc_obj$sensitivities)
pr_df <- as_tibble(pr_obj$curve)
names(pr_df) <- c("recall", "precision", "threshold")

roc_auc <- as.numeric(auc(roc_obj))
pr_auc <- pr_obj$auc.integral

p_roc <- ggplot(roc_df, aes(fpr, sensitivity)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey60") +
  geom_path(linewidth = 1.1, color = "#0072B2") +
  geom_point(data = res_001, aes(x = 1 - specificity, y = recall), inherit.aes = FALSE,
             shape = 21, size = 3, fill = "black", color = "white", stroke = 0.5) +
  annotate("text", x = 0.97, y = 0.08, label = sprintf("ROC AUC = %.3f", roc_auc),
           hjust = 1, fontface = "bold", size = 4.2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = "False-positive rate", y = "True-positive rate") +
  coord_equal() + theme_classic(base_size = 13)

p_pr <- ggplot(pr_df, aes(recall, precision)) +
  geom_path(linewidth = 1.1, color = "#D55E00") +
  geom_point(data = res_001, aes(x = recall, y = precision), inherit.aes = FALSE,
             shape = 21, size = 3, fill = "black", color = "white", stroke = 0.5) +
  annotate("text", x = 0.97, y = 0.08, label = sprintf("PR AUC = %.3f", pr_auc),
           hjust = 1, fontface = "bold", size = 4.2) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  labs(x = "Recall", y = "Precision") +
  coord_equal() + theme_classic(base_size = 13)

figure_2 <- (p_roc | p_pr) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 15))


figure_2

dir.create(here::here("figures"), showWarnings = FALSE)
ggsave(here::here("figures", "Figure_2_classifier_performance.png"), figure_2,
       width = 10, height = 5.2, units = "in", dpi = 600, bg = "white")

#-------------------------------------------
# 10. Guardar resultados
#-------------------------------------------

dir.create(here::here("outputs"), showWarnings = FALSE)
write.csv(results, here::here("outputs", "threshold_performance.csv"), row.names = FALSE)


# Métricas para el threshold utilizado
write.csv(res_001, here::here("outputs", "classifier_metrics_threshold_001.csv"), row.names = FALSE)

# Threshold con mayor F1 en el conjunto evaluado
write.csv(best, here::here("outputs", "classifier_best_threshold.csv"), row.names = FALSE)

# AUC y composición del conjunto de prueba
evaluation_summary <- data.frame(
  annotated_vocalizations = nrow(manual_labels),
  test_files = n_distinct(model_predictions$filename),
  total_windows = nrow(model_predictions),
  positive_windows = sum(model_predictions$truth == 1),
  negative_windows = sum(model_predictions$truth == 0),
  positive_window_prevalence = mean(model_predictions$truth == 1)
)

write.csv(evaluation_summary, here::here("outputs", "classifier_evaluation_summary.csv"), row.names = FALSE)

