# Preprocesado datos de ebird

# ============================================================
# ANÁLISIS eBird - Procnias tricarunculatus
# Filtrado, zerofill, elevación y modelo GAM altitudinal
# ============================================================

# -----------------------------
# 1. Paquetes
# -----------------------------

library(auk)
library(cageo)
library(dplyr)
library(sf)
library(terra)
library(lubridate)
library(ggplot2)
library(mgcv)
library(gratia)


# -----------------------------
# 2. Parámetros generales
# -----------------------------


crs_projected <- 5367
buffer_range_m <- 20000
max_elevation <- 2500
elev_band_width <- 100



f_ebd <- "data/raw/ebd_thwbel_smp_relFeb-2026.txt"
f_smp <- "data/raw/ebd_sampling_relFeb-2026.txt"


species_name <- "Procnias tricarunculatus"
taxonomy_year <- 2025

date_range <- c("2000-01-01", "2025-12-31")

ebd_out <- file.path("data/processed" , "ebird_campana_CR-PA.txt")
sampling_out <- file.path("data/processed", "ebird_campana_CR-PA_sampling.txt")

crs_projected <- 5367
buffer_range_m <- 20000
max_elevation <- 2500
elev_band_width <- 100


# -----------------------------
# 3. Área de estudio: Costa Rica y Panamá
# -----------------------------

bbox <- cageo::ca_outline_c |>
  filter(COUNTRY %in% c("Costa Rica", "Panama"))

borde <- cageo::ca_outline |>
  filter(COUNTRY %in% c("Costa Rica", "Panama")) |>
  st_union()


# -----------------------------
# 4. Filtrar datos eBird con auk
# -----------------------------

filters <- auk_ebd(f_ebd, file_sampling = f_smp) |>
  auk_species(species_name, taxonomy_version = taxonomy_year) |>
  auk_bbox(bbox) |>
  auk_date(date = date_range) |>
  auk_complete()

ebd_filtered <- auk_filter(
  filters,
  file = ebd_out,
  file_sampling = sampling_out,
  overwrite = TRUE
)


# -----------------------------
# 5. Zerofill
# -----------------------------

ebd_zf <- auk_zerofill(ebd_filtered)

obs <- collapse_zerofill(ebd_zf)



prtr_ebd <- auk::read_ebd("data/processed/ebird_campana_CR-PA.txt")
prtr_sampling <- auk::read_sampling("data/processed/ebird_campana_CR-PA_sampling.txt")

ebd_zf <- auk_zerofill(
  x = prtr_ebd,
  sampling_events = prtr_sampling)

obs <- collapse_zerofill(ebd_zf)

# -----------------------------
# 6. Convertir a sf y depurar espacialmente
# -----------------------------

obs <- obs |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) |>
  filter(lengths(st_intersects(geometry, borde)) > 0) |>
  select(
    country,
    locality_id,
    latitude,
    longitude,
    observation_date,
    species_observed,
    geometry
  ) |>
  mutate(
    observation_date = as.Date(observation_date),
    position = paste(longitude, latitude, sep = "_"),
    month = month(observation_date),
    year = year(observation_date),
    doy = yday(observation_date),
    woy = week(observation_date)
  ) |>
  st_drop_geometry()

write.csv(obs, file = "data/processed/obs_prtr_cr-pa.csv", row.names = FALSE)

# SE PUEDE INICIAR AQUI

# -----------------------------
# 7. Definir área accesible/relevante
#    Buffer de 20 km alrededor de presencias
# -----------------------------

# observaciones de Costa Rica y Panama

obs <- rio::import("data/processed/obs_prtr_cr-pa.csv")

obs <- obs |>
  st_as_sf(
    coords = c("longitude", "latitude"),
    crs = 4326,
    remove = FALSE
  ) 


# Mantener las observaciones dentro del rango de 20km alrededor de presencias reportadas

# crear el buffer
range_sf <- obs |>
  dplyr::filter(species_observed) |>
  st_transform(crs_projected) |>
  dplyr::summarise(geometry = st_union(geometry), .groups = "drop") |>
  st_buffer(dist = buffer_range_m) |>
  st_transform(st_crs(obs))

# Filtrar por el buffer
obs <- obs |>
  mutate(
    inside_range = lengths(st_within(geometry, range_sf)) > 0
  ) |>
  dplyr::filter(species_observed | inside_range)


# -----------------------------
# 8. Extraer elevación
# -----------------------------

library(stars)

elev <- methods::as(cageo::ca_elevation, "SpatRaster")

elev_values <- terra::extract(elev, terra::vect(obs)) |>
  dplyr::select(-ID) |>
  dplyr::rename(elevacion = 1)

obs <- bind_cols(obs, elev_values) |>
  mutate(
    elevacion = if_else(is.na(elevacion), 0, elevacion),
    elevacion = if_else(elevacion < 0, 0, elevacion)
  )


# -----------------------------
# 9. Agregar datos por año, semana y banda altitudinal
# -----------------------------

df_week <- obs |>
  st_drop_geometry() |>
  dplyr::filter(elevacion < max_elevation) |>
  mutate(
    elev_band = floor(elevacion / elev_band_width) * elev_band_width
  ) |>
  dplyr::group_by(year, elev_band, woy) |>
  dplyr::summarise(
    n_listas = n(),
    detecciones = sum(species_observed),
    .groups = "drop"
  ) |>
  dplyr::filter(n_listas > 0)


# modelo quasibinomial
# Ajustado con k = 12, 12, 7
mod_shift_quasi <-mgcv::gam(
  cbind(detecciones, n_listas - detecciones) ~
    te(
      woy,
      elev_band,
      year,
      bs = c("cc", "tp", "tp"),
      k = c(12, 12, 7)
    ),
  family = quasibinomial(link = "logit"),
  data = df_week,
  method = "REML",
  knots = list(
    woy = c(0.5, 53.5)
  )
)


# Evaluacion

capture.output(
  summary(mod_shift_quasi),
  file = "output/ebird_elevation_GAM_quasi_summary.txt"
)

set.seed(123)

k_quasi <- mgcv::k.check(
  mod_shift_quasi,
  subsample = 10000,
  n.rep = 1000
)

write.csv(
  k_quasi,
  "output/ebird_GAM_quasi_k_diagnostics.csv",
  row.names = TRUE
)

quasi_dispersion <- data.frame(
  diagnostic = c(
    "Estimated quasi-binomial scale",
    "Pearson dispersion"
  ),
  value = c(
    summary(mod_shift_quasi)$scale,
    sum(residuals(mod_shift_quasi, type = "pearson")^2) /
      df.residual(mod_shift_quasi)
  )
)

write.csv(
  quasi_dispersion,
  "output/ebird_GAM_quasi_dispersion.csv",
  row.names = FALSE
)

saveRDS(
  mod_shift_quasi,
  "output/ebird_elevation_GAM_quasi.rds"
)


# -----------------------------
# 11. Extraer efectos parciales
# -----------------------------

prediction_model <- mod_shift_quasi



sm <- smooth_estimates(prediction_model) |>
  mutate(year_round = floor(year)) |>
  filter(year_round %in% c(2000, 2005, 2010, 2015, 2020, 2025))



library(ggplot2)

# -----------------------------
# 12. Preparar etiquetas de meses
# -----------------------------

lineas_v <- seq(4.3, max(sm$woy), length.out = 12)
mid_points <- lineas_v - 2.15

df_meses <- data.frame(
  x = mid_points[1:12],
  y = 0,
  label = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
            "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)



# -----------------------------
# 13. Gráfico de efectos parciales
# -----------------------------

(p_effects <- ggplot(sm, aes(x = woy, y = elev_band, z = .estimate)) + 
    geom_raster(aes(fill = .estimate)) +
    geom_contour(color = "black", alpha = 0.3) +
    scale_fill_distiller(palette = "RdBu", direction = -1, limit = max(abs(sm$.estimate)) * c(-1, 1)) +
    facet_wrap(~year) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.position = "right") +
    labs(fill = "Partial effect",
         x = "Week of the year",
         y = "Elevation (masl)",
         title = "") +
    geom_text(data = df_meses, aes(x = x, y = y, label = label), 
              inherit.aes = FALSE, vjust = 0, size = 2, color = "gray40") +
    geom_vline(xintercept = lineas_v, linetype = "dashed", color = "gray80")
)


# Guardar y exportar
png(filename = "figures/efectos parciales.png", width = 25, height = 15, units = "cm", res = 300)
p_effects
dev.off()


write.csv(sm, file.path("output/efectos parciales.csv"), row.names = FALSE)



# Revision

# Datos agregados utilizados por el modelo
write.csv(
  df_week,
  "output/ebird_week_elevation_model_data.csv",
  row.names = FALSE
)

# Guardar el modelo
saveRDS(
  prediction_model,
  "output/ebird_elevation_GAM.rds"
)

# Resumen del modelo
capture.output(
  summary(prediction_model),
  file = "output/ebird_elevation_GAM_summary.txt"
)

# Concurvidad
concurvity_results <- mgcv::concurvity(prediction_model, full = TRUE)

write.csv(
  as.data.frame(concurvity_results),
  "output/ebird_elevation_GAM_concurvity.csv",
  row.names = TRUE
)

# Diagnóstico gráfico
png(
  "figures/ebird_elevation_GAM_check.png",
  width = 2400,
  height = 2000,
  res = 300
)

mgcv::gam.check(prediction_model)

dev.off()

# Información de la sesión
capture.output(
  sessionInfo(),
  file = "output/ebird_GAM_sessionInfo.txt"
)


# Resumen rango de variables
model_data_summary <- df_week |>
  dplyr::summarise(
    n_rows = dplyr::n(),
    first_year = min(year),
    last_year = max(year),
    min_week = min(woy),
    max_week = max(woy),
    min_elevation_band = min(elev_band),
    max_elevation_band = max(elev_band),
    total_checklists = sum(n_listas),
    total_detections = sum(detecciones),
    weeks_with_53 = sum(woy == 53)
  )

write.csv(
  model_data_summary,
  "output/ebird_GAM_data_summary.csv",
  row.names = FALSE
)




## 
# Grafico de elevacion

library(dplyr)
library(tidyr)
library(lubridate)
library(MASS)
library(purrr)

elev_ranges <- bind_rows(
  data.frame(
    elev_zone = "Lowlands (0–500 m)",
    elev_band = seq(0, 500, by = 100)
  ),
  data.frame(
    elev_zone = "Mid elevations (600–1100 m)",
    elev_band = seq(600, 1100, by = 100)
  ),
  data.frame(
    elev_zone = "Highlands (1200–2499 m)",
    elev_band = seq(1200, 2400, by = 100)
  )
)

newdata <- tidyr::crossing(
  year = sort(unique(df_week$year)),
  woy = 1:53,
  elev_ranges
) |>
  mutate(
    days_in_week = case_when(
      woy <= 52 ~ 7,
      lubridate::leap_year(year) ~ 2,
      TRUE ~ 1
    )
  ) |>
  group_by(year, elev_zone) |>
  mutate(
    averaging_weight = days_in_week / sum(days_in_week)
  ) |>
  ungroup()

# Matriz del predictor lineal
Xp <- predict(
  prediction_model,
  newdata = newdata,
  type = "lpmatrix"
)

# Coeficientes y matriz de covarianza
beta_hat <- coef(prediction_model)

V_beta <- vcov(
  prediction_model,
  unconditional = TRUE
)

# Simulaciones de los coeficientes
set.seed(123)

beta_sim <- MASS::mvrnorm(
  n = 2000,
  mu = beta_hat,
  Sigma = V_beta
)

# Índices de cada combinación año–zona
prediction_groups <- split(
  seq_len(nrow(newdata)),
  interaction(
    newdata$year,
    newdata$elev_zone,
    drop = TRUE
  )
)

elev_pred <- purrr::map_dfr(
  prediction_groups,
  function(idx) {
    
    X_group <- Xp[idx, , drop = FALSE]
    weights_group <- newdata$averaging_weight[idx]
    
    # Estimación puntual
    probabilities <- plogis(
      drop(X_group %*% beta_hat)
    )
    
    mean_probability <- sum(
      probabilities * weights_group
    )
    
    # Distribución simulada del promedio
    simulated_probabilities <- plogis(
      X_group %*% t(beta_sim)
    )
    
    simulated_means <- drop(
      crossprod(
        weights_group,
        simulated_probabilities
      )
    )
    
    data.frame(
      year = newdata$year[idx[1]],
      elev_zone = newdata$elev_zone[idx[1]],
      mean_pred = mean_probability,
      low = quantile(simulated_means, 0.025),
      high = quantile(simulated_means, 0.975)
    )
  }
)


elev_pred <- elev_pred |>
  mutate(
    elev_zone = factor(
      elev_zone,
      levels = c(
        "Lowlands (0–500 m)",
        "Mid elevations (600–1100 m)",
        "Highlands (1200–2499 m)"
      )
    )
  )



library(ggsci)
# Gráfico
(elev_plot <- ggplot(elev_pred, aes(x = year, y = mean_pred, color = elev_zone, fill = elev_zone)) +
    geom_line(linewidth = 1) +
    geom_ribbon(
      aes(ymin = low, ymax = high),
      alpha = 0.2,
      color = NA
    ) +
    theme_bw() +
    labs(
      x = "Year",
      y = "Mean predicted checklist reporting probability",
      color = "Elevation zone",
      fill = "Elevation zone"
    )
)

png(filename = "figures/prob_elev_range.png", width = 20, height = 10, units = "cm", res = 300)
elev_plot
dev.off()

write.csv(elev_pred, "output/prob_elev_range.csv", row.names = FALSE)

