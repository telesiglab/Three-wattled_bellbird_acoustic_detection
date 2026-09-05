# =========================================================
# ANALISIS INTEGRADO DE DETECCION DE PAJARO CAMPANA
# Variables de paisaje + acústica + ciencia ciudadana
# =========================================================

# -----------------------------
# 1. PAQUETES
# -----------------------------
library(tidyverse)
library(sf)
library(crgeo)
library(rio)
library(auk)
library(lubridate)
library(dbscan)
library(terra)
library(raster)
library(exactextractr)
library(landscapemetrics)
library(sjPlot)
library(ggsci)
library(glmmTMB)
library(purrr)
library(tidyr)
library(stars)
library(emmeans)
library(DHARMa)

# -----------------------------
# 2. PARAMETROS GENERALES
# -----------------------------

fecha_ini <- "2000-01-01"
fecha_fin <- "2025-12-31"

crs_ll <- 4326
crs_m <- 8908
buffer_dist <- 1000
cluster_eps <- 1000

path_bbox <- "data/raw/rectangulo_area.gpkg"
path_inat <- "data/raw/NCR_02_24.xlsx"
path_raster_cob <- "data/raw/Forestal 2023 8908.tif"
path_hfootprint <- "data/raw/HumanFootprint_CostaRica.tif"

#dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 3. FUNCION AUXILIAR: METRICAS DE PAISAJE
# -----------------------------
calc_metrics_buffer <- function(i, buffers_sf, raster_cat) {
  
  b <- buffers_sf[i, ]
  b_vect <- terra::vect(b)
  
  r_crop <- terra::crop(raster_cat, b_vect)
  r_mask <- terra::mask(r_crop, b_vect)
  
  if (is.null(r_mask) || all(is.na(terra::values(r_mask)))) {
    return(
      tibble::tibble(
        point_id = b$point_id,
        metric = NA_character_,
        class = NA,
        value = NA_real_
      )
    )
  }
  
  mets <- landscapemetrics::calculate_lsm(
    landscape = r_mask,
    what = c(
      "lsm_l_pland",
      "lsm_l_np",
      "lsm_l_pd",
      "lsm_l_ed",
      "lsm_l_shdi",
      "lsm_c_ca",
      "lsm_c_area_mn",
      "lsm_c_np",
      "lsm_c_pland",
      "lsm_c_ed",
      "lsm_l_ent"
    )
  )
  
  dplyr::mutate(mets, point_id = b$point_id)
}

# -----------------------------
# 4. AREA DE ESTUDIO
# -----------------------------
bbox <- sf::read_sf(path_bbox)

amistosa <- sf::read_sf("data/raw/borde_amistosa.gpkg")

cr_clip <- sf::st_intersection(crgeo::cr_outline, bbox)

# -----------------------------
# 5. DATOS ACUSTICOS
# -----------------------------

obsmod_raw <- rio::import("data/processed/detected_days.csv")  


obsmod_sum <- obsmod_raw %>%
  dplyr::mutate(
    season = ifelse(
      mes %in% c(9, 10, 11, 12, 1, 2),
      "Set-Feb (Non-Reproductive)",
      "Mar-Aug (Reproductive)"
    ),
    source = "Acoustic"
  ) %>%
  dplyr::rename(
    latitude = ycoord,
    longitude = xcoord
  ) %>%
  dplyr::group_by(site, latitude, longitude, season, source) %>%
  dplyr::summarise(
    ndet = sum(detected_days, na.rm = TRUE),
    nlistas = sum(sampled_days, na.rm = TRUE),
    bndet = ifelse(ndet >= 1, 1, 0),
    .groups = "drop"
  )

acoustic_sites <- sf::st_as_sf(obsmod_sum, coords = c("longitude", "latitude"), crs = crs_ll, remove = FALSE) |>
  dplyr::mutate(
    site_id = paste0("AC_", site),
    source = "Acoustic"
  ) |>
  dplyr::select(
    season, source, site_id, ndet, bndet, nlistas, geometry
  )



# ebird observations

obsebird <- rio::import("data/processed/obs_prtr_cr-pa.csv")

obsebird <- obsebird |>
  dplyr::filter(country == "Costa Rica") |>
  dplyr::select(
    observation_date,
    latitude,
    longitude,
    species_observed
  )


# -----------------------------
# 7. DATOS INATURALIST
# -----------------------------

# This data is not included because there is an agreement of not sharing this data.
# if you need this information ask naturalistacr@gmail.com 
# and uncomment this section

#inat <- rio::import(path_inat)

#inaturalist_observations <- inat %>%
#  dplyr::filter(
#    quality_grade == "research",
#    !is.na(observed_on),
#    !is.na(private_latitude),
#    !is.na(private_longitude),
#    observed_on > fecha_ini,
#    observed_on <= fecha_fin,
#    positional_accuracy <= 5000
#  )

#inaturalist_observations <- inaturalist_observations[!duplicated(inaturalist_observations[which(colnames(inaturalist_observations) %in% c("observed_on", "private_latitude", "private_longitude"))]), ]

#inaturalist_observations <- inaturalist_observations[, which(
#  colnames(inaturalist_observations) %in% c(
#    "scientific_name",
#    "observed_on",
#    "private_latitude",
#    "private_longitude",
#    "positional_accuracy"
#  )
#)]

#inaturalist_observations <- inaturalist_observations %>%
#  dplyr::rename(latitude = private_latitude, longitude = private_longitude, observation_date = observed_on) %>%
#  dplyr::select(latitude, longitude, observation_date)

#inaturalist_observations$species_observed <- TRUE

#ptsinat <- sf::st_as_sf(
#  inaturalist_observations,
#  coords = c("longitude", "latitude"),
#  crs = crs_ll,
#  remove = FALSE
#)

#inside_inat <- sf::st_intersects(ptsinat, bbox, sparse = FALSE)[, 1]
#ptsinat <- ptsinat[inside_inat, ]

#inatdb <- sf::st_drop_geometry(ptsinat)

# -----------------------------
# 8. UNIR EBIRD + INATURALIST
# -----------------------------

#obscc <- base::rbind(obsebird, inatdb)

obscc <- obsebird


obscc <- obscc %>%
  dplyr::mutate(
    observation_date = lubridate::as_date(observation_date),
    month = lubridate::month(observation_date),
    season = ifelse(
      month %in% c(9, 10, 11, 12, 1, 2),
      "Set-Feb (Non-Reproductive)",
      "Mar-Aug (Reproductive)"
    ),
    species_observed = as.numeric(species_observed),
    source = "Citizen science"
  )

ptsebi <- sf::st_as_sf(
  obscc,
  coords = c("longitude", "latitude"),
  crs = crs_ll,
  remove = FALSE
)

inside_obs <- sf::st_intersects(ptsebi, bbox, sparse = FALSE)[, 1]
ptsebi <- ptsebi[inside_obs, ]

# -----------------------------
# 9. CLUSTERING ESPACIAL DE CIENCIA CIUDADANA
# -----------------------------
obs_m <- sf::st_transform(ptsebi, crs_m)

coords <- sf::st_coordinates(obs_m)

cl <- dbscan::dbscan(
  coords,
  eps = cluster_eps,
  minPts = 1
)

obs_m$site_id <- paste0("CS_", cl$cluster)

obs_n <- obs_m %>%
  dplyr::distinct(
    observation_date,
    latitude,
    longitude,
    species_observed,
    .keep_all = TRUE
  )

citizen_sites <- obs_n %>%
  dplyr::group_by(season, source, site_id) %>%
  dplyr::summarise(
    ndet = sum(species_observed, na.rm = TRUE),
    bndet = ifelse(ndet >= 1, 1, 0),
    nlistas = dplyr::n(),
    .groups = "drop"
  ) %>%
  sf::st_centroid() %>%
  sf::st_transform(crs_ll)

# Mantener sitios positivos o negativos con esfuerzo suficiente
citizen_sites <- citizen_sites %>%
  dplyr::filter(ndet >= 1 | nlistas >= 20)

# -----------------------------
# 10. UNIR ACUSTICA + CIENCIA CIUDADANA
# -----------------------------
all_sites <- dplyr::bind_rows(
  citizen_sites %>%
    dplyr::select(season, source, site_id, ndet, bndet, nlistas, geometry),
  acoustic_sites %>%
    dplyr::select(season, source, site_id, ndet, bndet, nlistas, geometry)
)

inside_all <- sf::st_intersects(all_sites, cr_clip, sparse = FALSE)[, 1]
all_sites <- all_sites[inside_all, ]

all_sites$season <- factor(all_sites$season)
all_sites$source <- factor(all_sites$source)

table(all_sites$source)
table(all_sites$bndet)
summary(all_sites$nlistas)

# -----------------------------
# 11. METRICAS DE PAISAJE EN BUFFER DE 1 KM
# -----------------------------
raster_cob <- terra::rast(path_raster_cob)

all_sites_m <- sf::st_transform(all_sites, crs_m)

raster_cob_m <- terra::project(
  raster_cob,
  paste0("EPSG:", crs_m),
  method = "near"
)

all_sites_m <- all_sites_m %>%
  dplyr::mutate(point_id = dplyr::row_number())

buf_1km <- sf::st_buffer(all_sites_m, dist = buffer_dist)

raster_cob_m <- terra::as.factor(raster_cob_m)

levels(raster_cob_m) <- data.frame(
  value = c(1, 2, 3, 4),
  cover = c(
    "bosque_primario",
    "bosque_secundario",
    "plantacion",
    "otros"
  )
)

metricas_buffers <- purrr::map_dfr(
  seq_len(nrow(buf_1km)),
  calc_metrics_buffer,
  buffers_sf = buf_1km,
  raster_cat = raster_cob_m
)

metricas_wide <- metricas_buffers %>%
  dplyr::mutate(
    class_name = dplyr::case_when(
      class == 1 ~ "primario",
      class == 2 ~ "secundario",
      class == 3 ~ "plantacion",
      class == 4 ~ "otros",
      TRUE ~ NA_character_
    ),
    metric_name = ifelse(
      is.na(class_name),
      metric,
      paste0(metric, "_", class_name)
    )
  ) %>%
  dplyr::select(point_id, metric_name, value) %>%
  tidyr::pivot_wider(
    names_from = metric_name,
    values_from = value
  )

sites_metricas <- all_sites_m %>%
  dplyr::left_join(metricas_wide, by = "point_id")

# -----------------------------
# 12. VARIABLES AMBIENTALES
# -----------------------------
sites_metricas$season <- as.factor(sites_metricas$season)
sites_metricas$source <- as.factor(sites_metricas$source)

sites_metricas$ca_primario[is.na(sites_metricas$ca_primario)] <- 0
sites_metricas$ca_secundario[is.na(sites_metricas$ca_secundario)] <- 0

sites_metricas$bosque <- sites_metricas$ca_primario + sites_metricas$ca_secundario

sites_metricas$prec <- exactextractr::exact_extract(
  raster::raster(methods::as(crgeo::cr_prec, "SpatRaster")),
  buf_1km,
  "mean"
)

sites_metricas$elevacion <- exactextractr::exact_extract(
  raster::raster(methods::as(crgeo::cr_elevation, "SpatRaster")),
  buf_1km,
  "mean"
)

h_footprint <- raster::raster(path_hfootprint)

sites_metricas$h_footprint <- exactextractr::exact_extract(
  h_footprint,
  buf_1km,
  "mean"
)

hist(sites_metricas$nlistas)

# Limitar valores extremos de esfuerzo
#sites_metricas$nlistas[sites_metricas$nlistas > 1000] <- 1000

sites_metricas <- sites_metricas |>
  dplyr::mutate(
    nlistas_original = nlistas,
    prop_det = ndet / nlistas_original,
    nlistas_modelo = pmin(nlistas_original, 1000)
  )


# -----------------------------
# 13. ESCALAR PREDICTORES
# -----------------------------

# Predictores continuos que se incluirán en el modelo
predictores <- c(
  "ca_primario",
  "ca_secundario",
  "ent",
  "np",
  "elevacion",
  "h_footprint",
  "prec"
)

df_scaled <- sites_metricas %>%
  dplyr::mutate(
    season = as.factor(season),
    source = as.factor(source)
  ) %>%
  
  # Validar respuesta, esfuerzo y predictores
  dplyr::filter(
    dplyr::if_all(dplyr::all_of(predictores), ~ !is.na(.x)),
    !is.na(ndet),
    !is.na(nlistas_original),
    nlistas_original > 0,
    ndet >= 0,
    ndet <= nlistas_original
  ) %>%
  
  # Escalar solamente los predictores
  dplyr::mutate(
    dplyr::across(
      dplyr::all_of(predictores),
      ~ as.numeric(scale(.x))
    )
  )


# -----------------------------
# 14. MODELOS
# -----------------------------

mod1 <- glmmTMB::glmmTMB(
  prop_det ~ (
    ca_primario +
      ca_secundario +
      ent +
      np +
      elevacion +
      h_footprint +
      prec
  ) * season + source,
  weights = nlistas_modelo,
  family = binomial(link = "logit"),
  data = df_scaled
)


mod_sel <- mod1

summary(mod_sel)

# -----------------------------
# 15. REVISION DEL MODELO
# -----------------------------
sjPlot::plot_model(mod_sel)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("prec [all]", "season")
)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("elevacion [all]", "season")
)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("np [all]", "season")
)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("ent [all]", "season")
)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("h_footprint [all]", "season")
)

sjPlot::plot_model(
  mod_sel,
  type = "pred",
  terms = c("elevacion [all]", "season")
)

# Diagnóstico con DHARMa
DHARMa::simulateResiduals(mod_sel, plot = TRUE)

res_sim <- DHARMa::simulateResiduals(mod_sel, n = 1000)

plot(res_sim)
DHARMa::testDispersion(res_sim)
DHARMa::testZeroInflation(res_sim)



# -----------------------------
# 16. DATOS DE PREDICCION PARA GRAFICOS
# -----------------------------
variables <- c(
  "ca_primario",
  "ca_secundario",
  "ent",
  "np",
  "h_footprint",
  "prec"
)

mdata_pred <- data.frame()

for (i in variables) {
  
  v <- sjPlot::get_model_data(
    mod_sel,
    type = "pred",
    terms = c(paste0(i, " [all]"), "season")
  )
  
  v$variable <- i
  mdata_pred <- base::rbind(mdata_pred, v)
}

mdata_pred$variable[mdata_pred$variable == "ca_primario"] <- "mature_forest_area"
mdata_pred$variable[mdata_pred$variable == "ca_secundario"] <- "secondary_forest_area"
mdata_pred$variable[mdata_pred$variable == "np"] <- "ls_patch_number"
mdata_pred$variable[mdata_pred$variable == "ent"] <- "ls_entropy"
mdata_pred$variable[mdata_pred$variable == "h_footprint"] <- "human_footprint"
mdata_pred$variable[mdata_pred$variable == "prec"] <- "precipitation"

variablesnuevas <- unique(mdata_pred$variable)

# -----------------------------
# 17. PENDIENTES POR TEMPORADA
# -----------------------------
ann_list <- lapply(variables, function(v) {
  
  tr <- emmeans::emtrends(
    object = mod_sel,
    specs = ~ season,
    var = v
  )
  
  vcpos <- which(variables == tr@roles[["trend"]])
  tr_df <- as.data.frame(summary(tr, infer = c(TRUE, TRUE)))
  
  trend_col <- grep("\\.trend$", names(tr_df), value = TRUE)
  lcl_col <- grep("LCL|lower", names(tr_df), value = TRUE)
  ucl_col <- grep("UCL|upper", names(tr_df), value = TRUE)
  
  tr_df$sig <- ifelse(
    tr_df[[lcl_col]] > 0 | tr_df[[ucl_col]] < 0,
    "*",
    ""
  )
  
  tr_df$label <- paste0(
    "Slope = ", round(tr_df[[trend_col]], 3),
    " 95% CI: [",
    round(tr_df[[lcl_col]], 3), ", ",
    round(tr_df[[ucl_col]], 3), "]",
    tr_df$sig
  )
  
  tr_df$variable <- variablesnuevas[vcpos]
  tr_df
})

ann_df <- base::do.call(base::rbind, ann_list) %>%
  dplyr::rename(group = season)

# -----------------------------
# 18. POSICIONES DE TEXTO EN FACETAS
# -----------------------------
facet_pos <- mdata_pred %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    xmin = min(x, na.rm = TRUE),
    xmax = max(x, na.rm = TRUE),
    ymax = max(conf.high, na.rm = TRUE),
    ymin = min(conf.low, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    xrng = xmax - xmin,
    yrng = ymax - ymin,
    y_pos = ymax - 0.04 * yrng,
    x_left = xmin + 0.02 * xrng,
    x_right = xmax - 0.02 * xrng
  )

ann_plot <- ann_df %>%
  dplyr::left_join(facet_pos, by = "variable") %>%
  dplyr::mutate(
    x_pos = x_left,
    hjust_pos = 0,
    y_pos2 = dplyr::case_when(
      group == "Mar-Aug (Reproductive)" ~ y_pos,
      group == "Set-Feb (Non-Reproductive)" ~ y_pos - 0.10 * yrng,
      TRUE ~ y_pos
    )
  )

# -----------------------------
# 19. GRAFICO FINAL
# -----------------------------
grafico_final <- ggplot2::ggplot(
  data = mdata_pred,
  ggplot2::aes(x = x, y = predicted, color = group, fill = group)
) +
  ggplot2::geom_line(linewidth = 1) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = conf.low, ymax = conf.high),
    alpha = 0.2,
    color = NA
  ) +
  ggplot2::facet_wrap(~ variable, scales = "free") +
  ggplot2::geom_text(
    data = ann_plot,
    ggplot2::aes(
      x = x_pos,
      y = y_pos2,
      label = label,
      color = group,
      hjust = hjust_pos
    ),
    inherit.aes = FALSE,
    vjust = 1,
    size = 2.4,
    lineheight = 0.95,
    fontface = "bold",
    show.legend = FALSE
  ) +
  ggplot2::labs(
    x = "Predictor value (scaled)",
    y = "Predicted detection probability",
    color = "Season",
    fill = "Season"
  ) +
  ggplot2::theme_bw() +
  ggplot2::scale_fill_grey(start = 0.5, end = 0.5) +
  ggsci::scale_colour_npg() +
  ggplot2::theme(
    strip.background = ggplot2::element_rect(fill = "grey95", color = "grey70"),
    panel.spacing = grid::unit(1, "lines"),
    plot.margin = margin(8, 8, 8, 8)
  )

grafico_final

png(
  filename = "figures/Figure 4 variable influence.png",
  width = 30,
  height = 15,
  units = "cm",
  res = 300
)
grafico_final
dev.off()

write.csv(mdata_pred, file = "output/Variable influence.csv", row.names = FALSE)



write.csv(
  ann_df,
  "output/Variable influence slopes.csv",
  row.names = FALSE
)

capture.output(
  summary(mod_sel),
  file = "output/Variable influence model summary.txt"
)

capture.output(
  anova(mod1, mod2, mod3),
  file = "output/Variable influence model comparison.txt"
)

saveRDS(
  mod_sel,
  "output/Variable influence model.rds"
)

capture.output(
  sessionInfo(),
  file = "output/sessionInfo.txt"
)

# Resumen
effort_summary <- sites_metricas |>
  dplyr::group_by(source, season) |>
  dplyr::summarise(
    n_observations = dplyr::n(),
    n_over_1000 = sum(nlistas_original > 1000, na.rm = TRUE),
    percent_over_1000 = 100 * mean(nlistas_original > 1000, na.rm = TRUE),
    median_effort = median(nlistas_original, na.rm = TRUE),
    q95_effort = quantile(nlistas_original, 0.95, na.rm = TRUE),
    maximum_effort = max(nlistas_original, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  effort_summary,
  "output/Effort distribution summary.csv",
  row.names = FALSE
)
