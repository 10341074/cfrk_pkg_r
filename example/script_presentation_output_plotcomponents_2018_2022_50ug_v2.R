## Case study concentration: GGplot of predictions of classes
# - in kriging - v2 variogrammi sferici

rm(list = ls())
## -----------------------------------------------------------------------------
Current_Dir <- "/Users/giacomo/Library/CloudStorage/OneDrive-PolitecnicodiMilano/MOX/ricerca progetti/frk grins/CFRK/example"

## Import
if (length(Current_Dir)==0) {
  library(rstudioapi)
  Current_Dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}
source(paste0(Current_Dir, "/config.R"))

library(sf)
library(ggplot2)
library(tidyr)
library(dplyr)

Global_Use_Cov_Population = TRUE
Global_Use_Cov_Elevation = TRUE
obsFs = TRUE
BausCellsize = 1600

library(rstudioapi)
Current_Dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
source(paste0(Current_Dir, "/config.R"))
OUTPUT_FOLDER_PRESENTATION = paste0(OUTPUT_MAIN_DIR, '/figures-presentation-output_obsFs_', obsFs,'_BausCellsize_',BausCellsize ,'_2018-2022-50ug')
if(!dir.exists(OUTPUT_FOLDER_PRESENTATION)) {dir.create(OUTPUT_FOLDER_PRESENTATION, recursive = T)}
DIR_FIGURES = OUTPUT_FOLDER_PRESENTATION

data_air <- readRDS(paste0(OUTPUT_MAIN_DIR, "data_air.rds"))

DIR_FIGURES_FRK_INPUT = paste0(OUTPUT_MAIN_DIR, '/figures-quantiles-2018-2022-50ug-frk')
comuni = st_read(shapefile_comuni_filepath)

path <- paste0(DIR_FIGURES_FRK_INPUT, '/FRK_Pollutant_PM10_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '_prediction_levels3_Gaussian_PCA_FALSE_obsFs_', obsFs,'_BausCellsize_',BausCellsize ,'_IndexDrop_0/prediction_cofrk_comunale_component_all.csv')
result.frk <- read.csv(paste0(DIR_FIGURES_FRK_INPUT, "/output/prediction_cofrk_comunale_component_all.csv"))

aitch.pred <- paste0("R", 1:sum(grepl("_perc.pred", colnames(result.frk))), "_perc.pred")


result.frk <- dplyr::left_join(result.frk, comuni, by="PRO_COM") %>% st_as_sf

# frk
# result.frk.r1 = st_read(paste0(DIR_FIGURES_FRK_INPUT, '/FRK_Pollutant_PM10_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '_prediction_levels3_Gaussian_PCA_FALSE_obsFs_', obsFs,'_BausCellsize_',BausCellsize ,'_IndexDrop_0/prediction_cofrk_comunale_component_1_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '.gpkg'))
# result.frk.r2 = st_read(paste0(DIR_FIGURES_FRK_INPUT, '/FRK_Pollutant_PM10_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '_prediction_levels3_Gaussian_PCA_FALSE_obsFs_', obsFs,'_BausCellsize_',BausCellsize ,'_IndexDrop_0/prediction_cofrk_comunale_component_2_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '.gpkg'))



library(ggplot2)
library(scales)

source(paste(FRK_WORKING_DIR, 'init_data_air.R',  sep='/'))
domain <- GetDomain_Nord_Simplified()
# --------------------------------------------------------------------
# 0. Load stations
stations <- Load_Stations_Paper3Methods.Nord.Aggregated() %>%
  dplyr::select(-c("Longitude", "Latitude")) %>%
  st_as_sf() %>%
  rename(geometry = geom) %>%                # rinomina la colonna
  st_set_geometry("geometry") %>%            # imposta la colonna come geometria attiva
  st_set_crs(4326)                           # imposta il CRS (EPSG:4326)

#st_crs(stations) <- "EPSG:4326"
stations = stations %>% st_transform(st_crs(domain))
# --------------------------------------------------------------------


DoOutputPlotComponents <- function(resultlayer, method_name, output_folder, index, colors.values = NULL, colors.limits = NULL, colors.length.out = 10, do_points = T, add_name_par = "", datatitle=NULL)
{
    if(!is.null(colors.values)) {add_name = "_ScaledColors"} else {add_name = ""}
    add_name = paste0(add_name, add_name_par)

    x11(width = 8, height = 7)

    datafill = as.numeric(resultlayer[[paste0("R", index, "_perc.pred")]])
    if (is.null(datatitle)) {
      datatitle = paste0(method_name, " prediction of relative frequencies \n in interval I", as.character(i),
                         " =  ug/m3 on municipal base")
    }


    library(viridisLite)

    if(!is.null(colors.values))
    {
        # Genera 10 colori dalla scala "C"
        breaks = seq(min(colors.values, na.rm = T), max(colors.values, na.rm = T), length.out = colors.length.out)
        limits = c(min(colors.values, na.rm = T), max(colors.values, na.rm = T))
    } else {
        limits = colors.limits
        colors.values = seq(limits[1], limits[2], length.out = colors.length.out)
    }

    colors <- viridis(length(colors.values), option = "C", direction=1)

    gcomuni = ggplot() + geom_sf(data = resultlayer$geometry, aes(fill=datafill), linewidth = 0.1) +
        ggtitle(datatitle) +
        guides(fill=guide_legend(title="Frequency")) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, -0.4), "cm"),
                plot.title = element_text(size = 15, hjust=0.5),
                #legend.position = c(0.92, 0.6),
                legend.text = element_text(size = 10),   # Dimensione delle etichette della legenda
                legend.title = element_text(size = 10),   # Dimensione del titolo della legenda
                axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
                axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
                axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
                axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
        ) +
        scale_fill_gradientn(
            colours = colors,  # definizione dei colori
            values = scales::rescale(colors.values),  # punti specifici della scala
            limits = limits,  # limiti della scala
            breaks = colors.values,
            labels = label_number(0.001)) +
        labs(fill = "", x = "Longitude", y = "Latitude")
    if(do_points) {gcomuni = gcomuni + geom_sf(data=stations$geometry, aes(), shape=4, size=2, color="black", stroke=0.65)}
    gcomuni

    ggsave(paste(output_folder, paste0("/prediction_comuni_", method_name, "_c", as.character(index), add_name,
                                        '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
                                        '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
    graphics.off()
    gcomuni
}


PrintLimits <- function(resultlayer, index)
{
    v = as.numeric(resultlayer[[paste0("R", index, "_perc.pred")]])
    c(min(v, na.rm = T), max(v, na.rm = T))
}


PrintLimits(result.frk, 1)

PrintLimits(result.frk, 2)


if ((Global_Use_Cov_Elevation) && (Global_Use_Cov_Population))
{
    title = "FRK probability of PM10 of measuring \n over a year below 50 ug/m3 on municipal base"
    DoOutputPlotComponents(result.frk, "FRK", OUTPUT_FOLDER_PRESENTATION, 1, datatitle=title, colors.values = c(0.75, 0.83, 0.85, 0.88, 0.95, 1))

    title <- expression("FRK estimate: probability of PM"[10] * "\n" * "exceeding 50 " * mu * "g/m"^3 * " on municipal base")
    DoOutputPlotComponents(result.frk, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.values = c(0, 0.01, 0.05, 0.13, 0.14, 0.15, 0.16, 0.17, 0.244))
    DoOutputPlotComponents(result.frk, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.limits = c(0, 0.244))

    title = "FRK probability of exceeding 50 ug/m3 \n more than 35 days in a year on municipal base for PM10"
    DoOutputPlotComponents(result.frk, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.values = c(0.0959, 0.14, 0.16, 0.176, 0.244), add_name_par = "_threshold")
}

# ===========================================================================
# cut lombardia

result.frk.lombardia = result.frk %>% filter(COD_REG == 3)
#domain <- GetDomain_Lombardia()
stations = stations[unlist(st_contains(domain$geometry, stations$geometry)),]

if ((Global_Use_Cov_Elevation) && (Global_Use_Cov_Population))
{
  title = "FRK probability of PM10 of measuring \n over a year below 50 ug/m3 on municipal base"
  DoOutputPlotComponents(result.frk.lombardia, "FRK", OUTPUT_FOLDER_PRESENTATION, 1, datatitle=title, colors.values = c(0.75, 0.83, 0.85, 0.88, 0.95, 1), add_name_par = "_cutlombardia")

  title <- expression("FRK estimate: probability of PM"[10] * "\n" * "exceeding 50 " * mu * "g/m"^3 * " on municipal base")
  DoOutputPlotComponents(result.frk.lombardia, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.values = c(0, 0.01, 0.05, 0.13, 0.14, 0.15, 0.16, 0.17, 0.244), add_name_par = "_cutlombardia")
  DoOutputPlotComponents(result.frk.lombardia, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.limits = c(0, 0.244), add_name_par = "_cutlombardia")

  title = "FRK probability of exceeding 50 ug/m3 \n more than 35 days in a year on municipal base for PM10"
  DoOutputPlotComponents(result.frk.lombardia, "FRK", OUTPUT_FOLDER_PRESENTATION, 2, datatitle=title, colors.values = c(0.0959, 0.14, 0.16, 0.176, 0.244), add_name_par = "_cutlombardia_threshold")
}

# ===========================================================================
# sovrasoglia

# CFRK
result.frk$over = (result.frk[[paste0("R", 2, "_perc.pred")]] * 360) > 35
result.frk$giorni_sforamento <- as.logical(result.frk$over)


daily_over_35 <- ggplot(result.frk) +
    geom_sf(aes(fill = giorni_sforamento), linewidth = 0.1, show.legend = FALSE) +
    scale_fill_manual(
        values = c("TRUE" = 2, "FALSE" = "grey")
    ) +
    theme_light() +
    theme(
        plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 18, face = "bold", colour = "black")
    )

filename_sovrasoglia = paste0(OUTPUT_FOLDER_PRESENTATION, "/prediction_comuni_FRK_soprasoglia",
    '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
    '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png")
ggsave(
    filename = file.path(filename_sovrasoglia),
    plot = daily_over_35,
    width = 8,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "transparent",
    limitsize = FALSE
)

# ===========================================================================
# errore
#
# colors.values = seq(max(result.frk$sd, na.rm=T), min(result.frk$sd, na.rm=T), length.out = 10)
# library(viridisLite)
# colors <- viridis(length(colors.values), option = "C", direction=1)
#
# gsd = ggplot() + geom_sf(data = result.frk$geom, aes(fill = result.frk$sd), linewidth = 0.1) +
#     ggtitle("Error Method") +
#     guides(fill=guide_legend(title="Frequency")) +
#     theme(plot.margin = unit(c(0.1, 0.1, 0.1, -0.4), "cm"),
#             plot.title = element_text(size = 15, hjust=0.5),
#             #legend.position = c(0.92, 0.6),
#             legend.text = element_text(size = 10),   # Dimensione delle etichette della legenda
#             legend.title = element_text(size = 10),   # Dimensione del titolo della legenda
#             axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
#             axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
#             axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
#             axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
#     ) +
#     scale_fill_gradientn(
#         colours = colors,  # definizione dei colori
#         values = scales::rescale(colors.values),  # punti specifici della scala
#
#         breaks = colors.values,
#         labels = label_number(0.001)) +
#     labs(fill = "", x = "Longitude", y = "Latitude") + geom_sf(data=stations$geometry, aes(), shape=4, size=2, color="black", stroke=0.65)
#
# filename_output_sd = paste0(OUTPUT_FOLDER_PRESENTATION, "/prediction_comuni_FRK_sd",
#     '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
#     '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png")
#
# ggsave(
#     filename = file.path(filename_output_sd),
#     plot = gsd,
#     width = 8,
#     height = 8,
#     units = "in",
#     dpi = 300,
#     bg = "transparent",
#     limitsize = FALSE
# )

# Dati
df_citta <- data.frame(
  query = c("Imperia", "Argenta", "Bassano del Grappa", "Montichiari"),
  lat = c(43.95838, 44.61434, 45.76691, 45.41518),
  lon = c(7.866743, 11.833842, 11.734347, 10.390675)
)

# Creare sf POINT
points_sf <- st_as_sf(df_citta, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(st_crs(domain))

g <- ggplot() + geom_sf(data = domain$geometry, aes())

# Plot con soli punti e etichette
gmappa <- g +
  geom_sf(data = points_sf, shape=17, color = "black", size = 3) +
  geom_sf_text(data = points_sf, aes(label = query), nudge_y = 16000, size=4) +  # nudge_y in metri, regola a piacere
  coord_sf(crs = st_crs(domain))
gmappa


comuni_sample  <- comuni %>% filter(COMUNE %in% c("Imperia", "Argenta", "Bassano del Grappa", "Montichiari"))

result.frk.sample <- result.frk %>% filter(PRO_COM %in% comuni_sample$PRO_COM)
#df <- as.data.frame(result.frk.sample)[,aitch.pred]

library(dplyr)  # per select(), mutate(), row_number()
library(tidyr)  # per pivot_longer()

# aggiungiamo name a df_long (riciclato da result.frk.sample)

df_long <- result.frk.sample %>%
  pivot_longer(
    cols = c(R1_perc.pred, R2_perc.pred),
    names_to = "colonna",
    values_to = "valore"
  ) %>%
  mutate(
    idx = ifelse(colonna == "R1_perc.pred", 1, 2),
    ymin = ifelse(idx == 1, Y.pred.lower.1, Y.pred.lower.2),
    ymax = ifelse(idx == 1, Y.pred.upper.1, Y.pred.upper.2)
  )


df_long <- result.frk.sample %>%
  pivot_longer(
    cols = c(R1_perc.pred, R2_perc.pred),
    names_to = "colonna",
    values_to = "valore"
  ) %>%
  mutate(
    idx = ifelse(colonna == "R1_perc.pred", 1, 2),
    ymin = case_when(
      idx == 1 ~ Y.pred.lower.1,
      idx == 2 ~ Y.pred.lower.2
    ),
    ymax = case_when(
      idx == 1 ~ Y.pred.upper.1,
      idx == 2 ~ Y.pred.upper.2
    ),
    ymin_95 = case_when(
      idx == 1 ~ Y.pred.lower.95.1,
      idx == 2 ~ Y.pred.lower.95.2
    ),
    ymax_95 = case_when(
      idx == 1 ~ Y.pred.upper.95.1,
      idx == 2 ~ Y.pred.upper.95.2
    ),
    ymin_99 = case_when(
      idx == 1 ~ Y.pred.lower.99.1,
      idx == 2 ~ Y.pred.lower.99.2
    ),
    ymax_99 = case_when(
      idx == 1 ~ Y.pred.upper.99.1,
      idx == 2 ~ Y.pred.upper.99.2
    )
  )


ghist <- ggplot(df_long, aes(x = colonna, y = valore, color = colonna)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ymin_99, ymax = ymax_99), width = 0.15, linewidth = 0.8, color="grey0") +
  geom_errorbar(aes(ymin = ymin_95, ymax = ymax_95), width = 0.15, linewidth = 0.8, color="grey70") +
  geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.15, linewidth = 0.8) +
  # etichette sugli estremi
  #geom_text(aes(y = ymin, label = round(ymin, 2)), hjust = -1, size = 5) +
  #geom_text(aes(y = ymax, label = round(ymax, 2)), hjust = -1, size = 5) +
  # etichette intervallo dipendenti dalla colonna

  # Etichette per R1 sotto ymin
  geom_text(
    data = subset(df_long, colonna == "R1_perc.pred"),
    aes(y = ymax, x = 1.09,
        label = paste0(
          "\u03BC \u00B1 \u03C3 = [", round(ymax, 2), ", ", round(ymin, 2), "]\n",
          "\u03BC \u00B1 z(0.995) \u03C3 = [", round(ymax_99, 2), ", ", round(ymin_99, 2), "]")),
    hjust = 0, vjust = 1, size = 4, color = "black"
  ) +
  # Etichette per R2 sopra ymax
  geom_text(
    data = subset(df_long, colonna == "R2_perc.pred"),
    aes(y = ymax,x = 0.5,
        label = paste0("\u03BC \u00B1 \u03C3 = [", round(ymin, 2), ", ", round(ymax, 2), "]\n",
                       "\u03BC \u00B1 z(0.995) \u03C3 = [", round(ymin_99, 2), ", ", round(ymax_99, 2), "]")),
    hjust = 0, vjust = -0.4, size = 4, color = "black"
  )  +
  # # Etichette per R1 sotto ymin
  # geom_text(
  #   data = subset(df_long, colonna == "R1_perc.pred"),
  #   aes(y = ymax_99,
  #       label = paste0("[", round(ymax_99, 2), ", ", round(ymin_99, 2), "]")),
  #   hjust = -0.2, vjust = 1.4, size = 4, color = "black"
  # ) +
  # # Etichette per R2 sopra ymax
  # geom_text(
  #   data = subset(df_long, colonna == "R2_perc.pred"),
  #   aes(y = ymax_99,
  #       label = paste0("[", round(ymin_99, 2), ", ", round(ymax_99, 2), "]")),
  #   hjust = 1.2, vjust = -1.1, size = 4, color = "black"
  # )  +
  facet_wrap(~ COMUNE, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("R1_perc.pred" = "steelblue",
                                "R2_perc.pred" = "firebrick")) +
  scale_x_discrete(labels = c("R1_perc.pred" = data_air$levels[1],
                              "R2_perc.pred" = data_air$levels[2])) +
  theme_minimal(base_size = 12) +
  theme_bw() +
  theme(
    legend.position = "none",
    # titolo asse
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    # tick asse
    axis.text.x  = element_text(size = 14, angle = 0, vjust = 1),
    axis.text.y  = element_text(size = 14),
    # titolo facet
    strip.text   = element_text(size = 15, face = "bold")
  ) +
  labs(x = expression("Concentration classes " * mu * "g/" * m^3), y = "")

ghist

filename_output = paste0(OUTPUT_FOLDER_PRESENTATION, "/prediction_disthistograms_",
    '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
    '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png")

ggsave(
    filename = file.path(filename_output),
    plot = ghist,
    width = 8,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "transparent",
    limitsize = FALSE
)
cat("Saved ", filename_output, "\n")


## verifica che ilr su comuni è quella che
m <- as.data.frame(result.frk)[,c("R1_perc.pred", "R2_perc.pred")] %>% ilr()

g <- ggplot() + geom_sf(data=result.frk$geometry, aes(fill=as.numeric(m))) +
  geom_col() +
  scale_fill_viridis_c(option = "plasma", name = "ilr-values") +
  theme_minimal()


filename_output = paste0(OUTPUT_FOLDER_PRESENTATION, "/prediction_ilrval_",
                         '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
                         '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png")

ggsave(
  filename = file.path(filename_output),
  plot = g,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "transparent",
  limitsize = FALSE
)
cat("Saved ", filename_output, "\n")


g <- ggplot(data.frame(m = m), aes(x = m)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  theme_minimal()

filename_output = paste0(OUTPUT_FOLDER_PRESENTATION, "/prediction_distribilr_",
                         '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
                         '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png")

ggsave(
  filename = file.path(filename_output),
  plot = g,
  width = 8,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "transparent",
  limitsize = FALSE
)
cat("Saved ", filename_output, "\n")

