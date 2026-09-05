
here::i_am("scripts/03_predicted_seasonal_patterns.R")
library(dplyr)
library(sf)
library(ggplot2)
library(stars)
library(ggsci)
library(ggspatial)
library(here)



# Frecuencia de deteccion por hora

datos <- rio::import("data/raw/detecciones_rev_raven.csv") # Renombrar a validated_predicted_detections. 

datos$Hour12 <- as.numeric(floor(datos$Hour/10000))

horas <- as.data.frame(table(datos$Hour12))
horas$Var1 <- as.factor(horas$Var1)
(timeplot <- ggplot(data = horas, mapping = aes(x = Var1, y = Freq)) + geom_col(fill="#69b3a2", color = "darkgreen") + theme_minimal() + xlab("Hour") + ylab("Detection frequency"))

# Probabilidad de deteccion


df <- rio::import("data/processed/detected_days.csv")

df <- df %>%
  mutate(
    prop_det = detected_days / sampled_days,
    year = substr(month, 1, 4),
    month_num = as.integer(substr(month, 5, 6)),
    month_date = as.Date(paste0(month, "01"), format = "%Y%m%d"),
    month_lab = format(month_date, "%b %Y"),
    season = ifelse(mes %in% c(9,10,11,12,1,2), "Set-Feb (Non-Reproductive)", "Mar-Aug (Reproductive)")
  )

head(df)


landscape_data <- rio::import("data/processed/variables_sitios.csv")

df <- left_join(df, landscape_data |> 
                  dplyr::select(site, elevacion))

library(mgcv)

df$season <- as.factor(df$season)
df$site <- as.factor(df$site)

# gam prediction 
m <- gam(
  cbind(detected_days, sampled_days - detected_days) ~
    s(month_num, bs = "cc", k = 6) +
    s(elevacion, by = season, k = 5) +
    season +
    s(site, bs = "re"),
  family = binomial,
  data = df,
  method = "REML"
)


# Prediction over all combinations of site - month
pred_df <- df %>%
  distinct(site, xcoord, ycoord, elevacion, season, month, month_num, month_date, month_lab)
pred_df$pred <- predict(m, newdata = pred_df, type = "response")

write.csv(pred_df, "data/processed/pred_prob_det.csv", row.names = FALSE)

# Guardar objeto espacial 
#pred_pts <- st_as_sf(pred_df, coords = c("xcoord", "ycoord"), crs = 4326)
#names(pred_pts)
#st_write(obj =  pred_pts, dsn = "data/processed/pred_prob_det.gpkg", append = FALSE) 


# Figure 

# Load necessary layers
borde_amistosa <- read_sf("data/raw/borde_amistosa.gpkg")
mapa_background <- read_sf("data/processed/background_cantones_amistosa.gpkg")
frontera <- read_sf("data/processed/background_frontera_amistosa.gpkg")
elev <- read_stars("data/processed/elevacion_amistosa.tif")


pred_pts <- pred_pts %>%
  dplyr::rename(Det_prob = pred, 'Elevation (masl)' = elevacion)


(detection_plot <- ggplot() +
    geom_stars(data = elev, na.action = na.omit) +
    geom_sf(data = mapa_background, linewidth = 0.4, alpha = 0, color = "gray50") + geom_sf(data = frontera, color = "gray30", linewidth = 0.6, alpha = 0) +
    geom_sf(data = borde_amistosa, aes(linetype = "Amistosa Biological\nCorridor"), fill = "gray95", color = "darkred", linewidth = 0.6, alpha = 0.1)  +
    geom_sf_text(data = mapa_background, aes(label = NAME_2), size = 2, color = "gray40", position = "identity") +
    geom_sf_text(data = frontera, aes(label = COUNTRY), size = 4, color = "gray30", position = "identity")  + 
    geom_sf(data = pred_pts, aes(size = Det_prob), color = "steelblue") +
    scale_linetype_manual(name = NULL, values = c("Amistosa Biological\nCorridor" = "solid")) +
    coord_sf() +
    scale_fill_material("light-green", name = "Elevation (m a.s.l)") +
    facet_wrap(~season) + 
    theme_bw()  + xlab('') + ylab("") + 
    annotation_scale(location = "bl", width_hint = 0.3, tick_height = 0.2, line_width = 0.2,height = unit(0.15, "cm"))
  )


  png(filename = "figures/Mapa_detecciones_por_temporada.png", width = 25, height = 12, units = "cm", res = 300)
  detection_plot
  dev.off()
