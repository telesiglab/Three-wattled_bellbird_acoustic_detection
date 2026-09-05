# Preprocesado y manejo de datos

# ============================================================
# 1. Paquetes
# ============================================================

library(auk)
library(cageo)
library(crgeo)
library(dplyr)
library(here)
library(lubridate)
library(rio)
library(sf)
library(stars)
library(tidyr)


# ============================================================
# 2. Parámetros generales
# ============================================================

rootpath <- "C:/Paisaje acustico UNA/Proyecto AMISTOSA/Detecciones campana modelo/Text Files Modelo"
actualizar_detecciones_desde_txt <- FALSE

f_ebd <- here("data", "raw", "ebd_thwbel_smp_relFeb-2026.txt")
f_smp <- here("data", "raw", "ebd_sampling_relFeb-2026.txt")

species_name <- "Procnias tricarunculatus"
taxonomy_year <- 2025
date_range <- c("2000-01-01", "2025-12-31")

ebd_out <- here("data", "processed", "ebird_campana_CR-PA.txt")
sampling_out <- here("data", "processed", "ebird_campana_CR-PA_sampling.txt")

dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)


# ============================================================
# 3. Funciones
# ============================================================

# Calcular semana desde el inicio del muestreo
week_number <- function(date, inicio = "20230503") {
  dia0 <- ymd(inicio)
  diai <- ymd(date)
  gap <- as.numeric(diai - dia0)
  nweek <- floor(gap/7) + 1
  return(nweek)
}


# Resumir días muestreados por sitio y mes
resumir_dias_muestreados <- function(datos) {
  sitexmes <- datos %>%
    dplyr::group_by(Site, Month) %>%
    dplyr::summarise(Dias = n_distinct(Date), .groups = "drop")

  mescol <- sitexmes %>%
    pivot_wider(
      id_cols = "Site",
      values_from = "Dias",
      names_from = "Month"
    )

  return(list(sitexmes = sitexmes, mescol = mescol))
}


# ============================================================
# 4. Preprocesado de detecciones acústicas
# ============================================================

if (actualizar_detecciones_desde_txt) {
  files <- list.files(path = rootpath, pattern = "\\.txt$", full.names = TRUE)

  if (length(files) == 0) {
    stop("No se encontraron archivos .txt en rootpath")
  }

  datos <- data.frame()
  for (i in files) {
    pdf <- rio::import(i)
    datos <- dplyr::bind_rows(datos, pdf)
  }

  datos$Site <- NA
  datos$Date <- NA
  datos$Hour <- NA
  datos$Month <- NA
  datos$Week <- NA

  for (i in seq_len(nrow(datos))) {
    pedazos <- unlist(strsplit(datos$`Begin File`[i], split = "_"))

    if (is.na(as.numeric(pedazos[1]))) {
      datos$Site[i] <- pedazos[1]
      datos$Date[i] <- pedazos[2]
      datos$Hour[i] <- sub(x = pedazos[3], pattern = "\\.mp3$", replacement = "")
      datos$Month[i] <- substr(pedazos[2], 1, 6)
    } else {
      datos$Site[i] <- pedazos[2]
      datos$Date[i] <- pedazos[3]
      datos$Hour[i] <- sub(x = pedazos[4], pattern = "\\.mp3$", replacement = "")
      datos$Month[i] <- substr(pedazos[3], 1, 6)
    }

    datos$Week[i] <- week_number(datos$Date[i])
  }

  write.csv(
    datos,
    file = here("data", "raw", "detecciones_rev_raven.csv"),
    row.names = FALSE
  )
} else {
  datos <- rio::import(here("data", "raw", "detecciones_rev_raven.csv"))
}

resumen_muestreo <- resumir_dias_muestreados(datos)
sitexmes <- resumen_muestreo$sitexmes
mescol <- resumen_muestreo$mescol

write.csv(
  sitexmes,
  file = here("data", "processed", "sitexmes.csv"),
  row.names = FALSE
)

write.csv(
  mescol,
  file = here("data", "processed", "mescol.csv"),
  row.names = FALSE
)


# ============================================================
# 5. Capas espaciales del área de estudio
# ============================================================

rectangulo <- read_sf(here("data", "raw", "rectangulo_area.gpkg"))
poligono_area <- st_as_sfc(rectangulo)

mapa_background <- st_intersection(crgeo::cr_cantons_c, poligono_area)

st_write(
  mapa_background,
  dsn = here("data", "processed", "background_cantones_amistosa.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

frontera <- cageo::ca_outline_c %>%
  filter(COUNTRY == "Panama") %>%
  st_intersection(poligono_area)

frontera$COUNTRY <- "Panamá"

st_write(
  frontera,
  dsn = here("data", "processed", "background_frontera_amistosa.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

elev_crop <- st_crop(
  cageo::ca_elevation,
  st_bbox(poligono_area)
)

elev <- elev_crop[poligono_area]

write_stars(
  elev,
  dsn = here("data", "processed", "elevacion_amistosa.tif")
)


# ============================================================
# 6. Preprocesado de datos de eBird
# ============================================================

# Área de estudio: Costa Rica y Panamá
bbox <- cageo::ca_outline_c |>
  filter(COUNTRY %in% c("Costa Rica", "Panama"))

borde <- cageo::ca_outline |>
  filter(COUNTRY %in% c("Costa Rica", "Panama")) |>
  st_union()


# Filtrar datos eBird con auk
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


# Zerofill
ebd_zf <- auk_zerofill(ebd_filtered)
obs <- collapse_zerofill(ebd_zf)


# Convertir a sf y depurar espacialmente
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

write.csv(
  obs,
  file = here("data", "processed", "obs_prtr_cr-pa.csv"),
  row.names = FALSE
)
