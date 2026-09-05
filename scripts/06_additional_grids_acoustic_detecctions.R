# ============================================================
# NUEVAS DETECCIONES ACÚSTICAS NO REPORTADAS EN CIENCIA CIUDADANA
# Comparación entre eBird y modelo acústico por grilla H3 y mes
# ============================================================

# -----------------------------
# 1. Paquetes
# -----------------------------

library(auk)
library(dplyr)
library(sf)
library(ggplot2)
library(rio)
library(lubridate)
library(crhexgrids)
library(crgeo)
library(cageo)


# -----------------------------
# 2. Parámetros generales
# -----------------------------

thresh_pos <- 1
thresh_neg_days <- 20

amistosa <- sf::read_sf("data/raw/borde_amistosa.gpkg")

# -----------------------------
# 3. Área de estudio
# -----------------------------

bbox <- cageo::ca_outline_c |>
  filter(COUNTRY %in% c("Costa Rica", "Panama"))

borde <- cageo::ca_outline |>
  filter(COUNTRY %in% c("Costa Rica", "Panama")) |>
  st_union()


# -----------------------------
# 4. Grilla H3 de Costa Rica
# -----------------------------

grid <- crhexgrids::cr_hex_grid_res_6

grid <- grid[
  lengths(st_intersects(grid, crgeo::cr_outline_c)) > 0,
]


# -----------------------------
# leer edb
# -----------------------------


prtr_ebd <- auk::read_ebd("data/processed/ebird_campana_CR-PA.txt")
prtr_sampling <- auk::read_sampling("data/processed/ebird_campana_CR-PA_sampling.txt")

ebd_zf <- auk_zerofill(
  x = prtr_ebd,
  sampling_events = prtr_sampling)

ebd_zf_df <- collapse_zerofill(ebd_zf)


# -----------------------------
# 7. Convertir eBird a sf y depurar espacialmente
# -----------------------------

obs <- st_as_sf(
  ebd_zf_df,
  coords = c("longitude", "latitude"),
  crs = 4326,
  remove = FALSE
)

obs <- obs[
  lengths(st_intersects(obs, borde)) > 0,
]

obs <- obs |>
  mutate(
    observation_date = as.Date(observation_date),
    month = lubridate::month(observation_date),
    position = paste(longitude, latitude, sep = "_"),
    species_observed = as.numeric(species_observed)
  )



# ***********************

# Asignar directamente cada lista eBird a una celda H3
obs_grid <- sf::st_join(
  obs,
  grid |> dplyr::select(h3_address),
  join = sf::st_within,
  left = FALSE
)

# Resumen histórico por celda y mes calendario, acumulando 2000–2025
ebird_cell <- obs_grid |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(h3_address)) |>
  dplyr::group_by(h3_address, month) |>
  dplyr::summarise(
    nlistas_ebird = dplyr::n_distinct(checklist_id),
    ndias_ebird = dplyr::n_distinct(observation_date),
    
    positive_ebird = sum(species_observed == 1, na.rm = TRUE),
    pres_ebird = as.integer(any(species_observed == 1, na.rm = TRUE)),
    
    first_positive_date = if(any(species_observed == 1, na.rm = TRUE)) {
      min(observation_date[species_observed == 1], na.rm = TRUE)
      } else {as.Date(NA)},
    last_positive_date = if(any(species_observed == 1, na.rm = TRUE)) {
      max(observation_date[species_observed == 1],na.rm = TRUE)
      } else {
      as.Date(NA)
    }, 
    .groups = "drop"
  )


# -----------------------------
# 8. Crear resumen de ciencia ciudadana
# -----------------------------
# Unidad: punto de muestreo + mes

muestreos_summary <- obs |>
  group_by(geometry, month) |>
  summarise(
    nlistas = n_distinct(checklist_id),
    ndias = n_distinct(observation_date),
    positive = sum(species_observed == 1, na.rm = TRUE),
    negative = sum(species_observed == 0, na.rm = TRUE),
    presencia = if_else(positive >= thresh_pos, "positivo", "negativo"),
    presencia_bn = if_else(positive >= thresh_pos, 1, 0),
    .groups = "drop"
  ) |>
  dplyr::filter(
    positive >= thresh_pos | ndias > thresh_neg_days
  )


# -----------------------------
# 9. Asignar ciencia ciudadana a grillas H3
# -----------------------------

ebird_grids <- st_join(
  muestreos_summary,
  grid |> dplyr::select(h3_address),
  join = st_within
)


# -----------------------------
# 10. Cargar y preparar detecciones acústicas
# -----------------------------

model_data <- rio::import("data/processed/detected_days.csv") |>
  mutate(
    prop_det = detected_days / sampled_days,
    year = as.integer(substr(month, 1, 4)),
    month_num = as.integer(substr(month, 5, 6)),
    month_date = as.Date(paste0(month, "01"), format = "%Y%m%d"),
    month_lab = format(month_date, "%b %Y"),
    season = if_else(
      month_num %in% c(9, 10, 11, 12, 1, 2),
      "Set-Feb (Non-Reproductive)",
      "Mar-Aug (Reproductive)"
    ),
    presencia_bn = if_else(detected_days >= 1, 1, 0)
  )


# -----------------------------
# 11. Convertir modelo acústico a sf y asignar grilla H3
# -----------------------------

model_sf <- st_as_sf(
  model_data,
  coords = c("xcoord", "ycoord"),
  crs = 4326,
  remove = FALSE
)

model_grids <- st_join(
  model_sf,
  grid |> dplyr::select(h3_address),
  join = st_within
)


# -----------------------------
# 12. Resumir presencia eBird por celda y mes
# -----------------------------



## Celdas de acustica nuevo

model_cell <- model_grids |>
  sf::st_drop_geometry() |>
  dplyr::filter(!is.na(h3_address)) |>
  dplyr::group_by(h3_address, month = month_num) |>
  dplyr::summarise(
    pres_acoustic = as.integer(
      any(detected_days >= 1, na.rm = TRUE)
    ),
    detected_days = sum(detected_days,na.rm = TRUE),
    sampled_days = sum(sampled_days, na.rm = TRUE),
    acoustic_years = paste(sort(unique(year[detected_days >= 1])), collapse = ";"),
    .groups = "drop"
  )





cells_combined <- model_cell |>
  dplyr::filter(pres_acoustic == 1) |>
  dplyr::left_join(ebird_cell, by = c("h3_address", "month")
  ) |>
  dplyr::mutate(
    ebird_sampling_status = dplyr::case_when(
      is.na(nlistas_ebird) ~ "No complete eBird checklists",
      pres_ebird == 1 ~ "Represented in eBird archive",
      TRUE ~ "eBird checklists without a positive record"),
    pres_ebird = dplyr::coalesce(pres_ebird, 0L),
    nlistas_ebird = dplyr::coalesce(nlistas_ebird, 0L),
    ndias_ebird = dplyr::coalesce(ndias_ebird, 0L),
    positive_ebird = dplyr::coalesce(positive_ebird,0L),
    acoustic_only_record = pres_ebird == 0
  )



## Nuevo resume

monthly_summary <- cells_combined |>
  dplyr::group_by(month) |>
  dplyr::summarise(
    acoustic_positive_cells = dplyr::n(),
    represented_in_ebird = sum(pres_ebird == 1),
    acoustic_only_cells = sum(acoustic_only_record),
    proportion_acoustic_only = acoustic_only_cells / acoustic_positive_cells,
    acoustic_only_without_ebird_sampling = sum(acoustic_only_record & ebird_sampling_status == "No complete eBird checklists"),
    acoustic_only_with_ebird_sampling = sum(acoustic_only_record & ebird_sampling_status == "eBird checklists without a positive record"),
    .groups = "drop"
  )

overall_summary <- cells_combined |>
  dplyr::summarise(
    acoustic_positive_cell_months = dplyr::n(),
    acoustic_only_cell_months = sum(acoustic_only_record),
    proportion_acoustic_only = acoustic_only_cell_months / acoustic_positive_cell_months,
    unique_acoustic_cells = dplyr::n_distinct(h3_address),
    unique_acoustic_only_cells = dplyr::n_distinct(h3_address[acoustic_only_record]),
    acoustic_only_without_ebird_sampling = sum(acoustic_only_record & ebird_sampling_status == "No complete eBird checklists"),
    acoustic_only_with_ebird_sampling = sum(acoustic_only_record & ebird_sampling_status == "eBird checklists without a positive record")
  )

acoustic_only_table <- cells_combined |>
  dplyr::filter(acoustic_only_record) |>
  dplyr::arrange(month, h3_address)






## Exports 

write.csv(cells_combined, "output/acoustic_ebird_historical_comparison.csv", row.names = FALSE)

write.csv(monthly_summary, "output/acoustic_contribution_monthly_summary.csv", row.names = FALSE)

write.csv(overall_summary, "output/acoustic_contribution_overall_summary.csv", row.names = FALSE)

write.csv(acoustic_only_table, "output/acoustic_only_cell_month_records.csv", row.names = FALSE)



##
# ============================================================
# FIGURA 7: contribución espacial y estacional de la acústica
# ============================================================



# Registros históricos de eBird por celda y mes, 2000–2025
ebird_historical <- ebird_cell |>
  filter(pres_ebird == 1) |>
  distinct(h3_address, month) |>
  transmute(
    h3_address, month,
    record_status = "Historical eBird record (2000–2025)"
  )

# Detecciones acústicas sin registro histórico equivalente en eBird
acoustic_contribution <- cells_combined |>
  filter(as.logical(acoustic_only_record)) |>
  distinct(h3_address, month) |>
  transmute(
    h3_address, month,
    record_status = "Additional acoustic record (2023–2024)"
  )

# Etiquetas para las 12 facetas
month_labels <- tibble(month = 1:12) |>
  left_join(monthly_summary, by = "month") |>
  mutate(
    across(c(acoustic_positive_cells, acoustic_only_cells), ~ replace_na(.x, 0)),
    panel_label = sprintf(
      "%s (%d additional)",
      month.abb[month], acoustic_only_cells)
  )

panel_levels <- month_labels$panel_label

# Celdas donde se realizó muestreo acústico, para definir la extensión
acoustic_extent <- grid |>
  semi_join(model_cell |> distinct(h3_address), by = "h3_address")

amistosa <- st_transform(amistosa, st_crs(grid))

# rectangulo

rectangulo = st_read("data/raw/rectangulo_area.gpkg")

library(basemaps)
library(terra)
library(tidyterra)
mapa_base <- basemap_terra(ext = rectangulo, map_service = "osm", map_type = "streets")
mapa_base <- terra::project(mapa_base, "EPSG:4326")

  
  

# Definir extensión del mapa


map_bbox_sf <- rectangulo

grid_context <- grid[lengths(st_intersects(grid, map_bbox_sf)) > 0, ]

# Unir la grilla con los registros históricos y acústicos
ebird_historical_sf <- grid |>
  inner_join(ebird_historical, by = "h3_address") |>
  st_filter(map_bbox_sf, .predicate = st_intersects) |>
  left_join(month_labels |> select(month, panel_label), by = "month") |>
  mutate(panel_label = factor(panel_label, levels = panel_levels))

acoustic_contribution_sf <- grid |>
  inner_join(acoustic_contribution, by = "h3_address") |>
  left_join(month_labels |> select(month, panel_label), by = "month") |>
  mutate(panel_label = factor(panel_label, levels = panel_levels))


figure_7 <- ggplot() +
  geom_spatraster_rgb(data = mapa_base) +
  geom_sf(
    data = grid_context,
    fill = "grey97", color = "grey88", linewidth = 0.10, alpha = 0.5) +
  geom_sf(
    data = amistosa, aes(color = "AmistOsa Biological Corridor"), fill = NA, linewidth = 0.55
    ) +
  geom_sf(
    data = ebird_historical_sf,
    aes(fill = record_status),
    color = "white", linewidth = 0.25, alpha = 0.75
  ) +
  geom_sf(
    data = acoustic_contribution_sf,
    aes(fill = record_status),
    color = "#8C510A", linewidth = 0.60, alpha = 0.95
  ) +
  facet_wrap(~ panel_label, ncol = 4, drop = FALSE) +
  scale_fill_manual(
    values = c(
      "Historical eBird record (2000–2025)" = "#56B4E9",
      "Additional acoustic record (2023–2024)" = "#E69F00"
    ),
    breaks = c(
      "Historical eBird record (2000–2025)",
      "Additional acoustic record (2023–2024)"
    )
  ) +
  scale_color_manual(
    values = c("AmistOsa Biological Corridor" = "black")
  ) +
  coord_sf(xlim = x_limits, ylim = y_limits, expand = FALSE) +
  labs(fill = "", colour = "") +
  theme_void(base_size = 10) +
  theme(
    strip.background = element_rect(
      fill = "grey92", color = "grey45", linewidth = 0.3
    ),
    strip.text = element_text(face = "bold", size = 8.5),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.1, "cm"),
    panel.spacing = grid::unit(0.35, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

figure_7

ggsave(
  "figures/Figure_7_acoustic_contribution.png",
  figure_7, width = 24, height = 18, units = "cm",
  dpi = 600, bg = "white"
)






