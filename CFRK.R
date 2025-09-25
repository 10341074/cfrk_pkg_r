source(paste(CFRK_PACKAGE_DIR, 'CFRK_utils.R',  sep='/'))

# CFRK method
no_iterations = 20

PlotSymbols = function(dataset, colname, Do.Log = T, Use.Pop = T, title="", labels=c(), layer_popolazione=NULL) {
  plotsym_vals_shape = c(1, 1, 2, 5, 5)
  plotsym_vals_size = c(3, 1, 1, 1, 3)

  dataset$interval = as.integer(findInterval(dataset[[colname]], c(-Inf,unlist(quantile(dataset[[colname]], c(0.2,0.4,0.6, 0.8))), Inf)))
  dataset$size = plotsym_vals_size[dataset$interval]
  dataset$shape = as.factor(plotsym_vals_shape[dataset$interval])
  if (Use.Pop) {
    g <- ggplot() + geom_sf(data=layer_popolazione$geom, aes(fill=1e6 * layer_popolazione$Popolazione.P1 / as.numeric(st_area(layer_popolazione$geom))), color = "white", linewidth = 0.001, alpha=0.7) +
      scale_fill_distiller(palette = "Blues", direction=-1)
  } else {
    g <- ggplot() + geom_sf(data=regioni$geometry[3], aes())
  }
  g <- g + geom_point(data=dataset, aes(x, y), shape=dataset$shape, size=dataset$size, color="black", stroke=0.5) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust=0.5)) +
    guides(fill=guide_legend(title=""))
  if (length(labels) > 0) {
    g <- g + geom_text(aes(x=dataset$x, y=dataset$y, label = labels), vjust = -1)
  }
  g
}


CFRK <- function(data.comp, flag_use_pca, flag_bau_comunali, flag_obs_fs, baus.cellsize, frk.args,
  layer_popolazione, raster.elevazione,
  domain, comuni, cov.P1 = NULL, GridBAUs = NULL,
  DestName.IndexDrop = 0, DestName.Pollutant = "", fitted.locations = NULL, obs.test = NULL,
  comp.levels = c(), colnames.aitch = c(), colnames.euclid = c(), figures_folder = "./figures/",
  Global_Use_Cov_Population = TRUE, Global_Use_Cov_Elevation = FALSE, Use_Global_Colors = TRUE, Colors_Breaks = NULL, DoPLOT = TRUE, FLAG_NO_PLOT = FALSE)
{
  colnames_pred = paste0(colnames.euclid, '.pred')
  colnames.aitch.pred = paste0(colnames.aitch, '.pred')

  folder_output = GetFRKFolderOutput (figures_folder,
   DestName.Pollutant, Global_Use_Cov_Population, Global_Use_Cov_Elevation,
   frk.args, flag_use_pca, flag_obs_fs, baus.cellsize, 0)

  cat("\nDoFixedRankCokriging: Using PCA = ", flag_use_pca, "\n")
  cat("Folder output = ", folder_output, "\n")

  # ----------------------------------------------------------------------------

  library(magrittr)
  library(FRK) # for carrying out FRK

  # 1. PCA # define 'data_air' before
  data_for_kriging <- PrepareData_For_Kriging(
    data.comp, colnames.euclid, flag_use_pca)

  data.krig <- data_for_kriging[["data.krig"]]
  pca_output <- data_for_kriging[["pca_output"]]

  # ------------------------------------- save plot ----------------------------
  if (!FLAG_NO_PLOT) {
    x11(width = 8, height = 7)
    for (i in (1:length(colnames.euclid)))
    {
      PlotSymbols(as.data.frame(data.krig), colnames.euclid[i], Use.Pop = T,
        title=paste0("Score ", as.character(i), " in Euclid"),
        layer_popolazione = layer_popolazione) # score 1 pca
      ggsave(paste0(folder_output, "/obs_sym_Score", as.character(i), "_in_Euclid.png"))
    }
    graphics.off()
  }

  # ============================================================================
  cor(as.data.frame(data.krig)[colnames.euclid])
  # ----------------------------------------------------------------------------

  resultFRK <- DoFRKPrediction(domain, GridBAUs, frk.args, data.comp,
    data_for_kriging[['data.krig']], colnames.euclid, comuni,
    Global_Use_Cov_Population, Global_Use_Cov_Elevation, no_iterations)



  predictionbaus.values = resultFRK$predictionbaus.values
  fitted.values = resultFRK$fitted.values
  comuni.values = resultFRK$comuni.values

  # ----------------------------------------------------------------------------
  # Invert PCA
  if(flag_use_pca)
  {
    originalbasis.data <- ScoresToVars(
      as.matrix(as.data.frame(predictionbaus.values)[,colnames_pred]),
      pca_result$pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)

    originalbasis.fitted.values <- ScoresToVars(
      as.matrix(as.data.frame(fitted.values)[,colnames_pred]),
      pca_result$pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)

    originalbasis.municipal.values <- ScoresToVars(
      as.matrix(as.data.frame(comuni.values)[,colnames_pred]),
      pca_result$pca_output, colnames=colnames.euclid) %>%
      set_colnames(colnames_pred)


    predictionbaus.values <-  SpatialPixelsDataFrame(coordinates(predictionbaus.values), data = originalbasis.data)
    fitted.values <- SpatialPointsDataFrame(coordinates(fitted.locations), data=originalbasis.fitted.values)
  } else {
    originalbasis.municipal.values = comuni.values
  }

  # ------------------------------------------------------------------------------
  # INVERT ILR
  # 1/3 invert data grid
  data.to.revert.prediction = as.data.frame(as.data.frame(predictionbaus.values)[,colnames_pred])
  colnames(data.to.revert.prediction) <- colnames_pred
  prediction.conc = RevertILR(data.to.revert.prediction, predictionbaus.values$x, predictionbaus.values$y)
  # seleziono solo i pixel dentro il dominio (es regione lombardia)
  prediction.conc <- prediction.conc[unlist(st_contains(data_air$domain,
    st_as_sf(prediction.conc, coords=c("x","y"), crs=st_crs(data.comp)$epsg)$geometry)),]

  # 2/3 invert data fitted
  data.to.revert.fitted = as.data.frame(as.data.frame(fitted.values)[colnames_pred])
  colnames(data.to.revert.fitted) <- colnames_pred
  comp.fitted.df = RevertILR(data.to.revert.fitted, fitted.values$x, fitted.values$y)
  comp.fitted.df = cbind(comp.fitted.df, fitted.values[colnames_pred])

  # 3/3 invert data comuni
  data.to.revert.comuni = as.data.frame(as.data.frame(originalbasis.municipal.values)[,colnames_pred])
  colnames(data.to.revert.comuni) <- colnames_pred
  comp.comuni.conc = comuni
  comp.comuni.conc[,colnames.aitch.pred] = ilrInv(data.to.revert.comuni)

  list(
    comp.prediction.values = prediction.conc,
    comp.fitted.values = comp.fitted.df,
    comp.comuni.values = comp.comuni.conc,
    Slist=resultFRK$Slist)
}


CFRK_GetSD <- function(result, data_air_dataset, comuni, index = 1) {
  prediction = SpatialPixelsDataFrame(
    as.matrix(as.data.frame(result$Slist[[index]]$BAUs_prediction_df)[,c('x','y')]),
    data = data.frame(Y.pred.sd = result$Slist[[index]]$BAUs_prediction_df$sd),
    proj4string = CRS(paste("EPSG:", st_crs(data_air_dataset)$epsg)))

  GetAvg(prediction, comuni, "Y.pred.sd", data_air_dataset)
}


CFRK_GetErrorInterval <- function (result, data_air_dataset, comuni, index = 1) {
  prediction = SpatialPixelsDataFrame(
    as.matrix(as.data.frame(result$Slist[[index]]$BAUs_prediction_df)[,c('x','y')]),
    data = data.frame(
      Y.pred.upper = result$Slist[[index]]$BAUs_prediction_df$mu + result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower = result$Slist[[index]]$BAUs_prediction_df$mu - result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.upper.95 = result$Slist[[index]]$BAUs_prediction_df$mu + qnorm(0.975) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower.95 = result$Slist[[index]]$BAUs_prediction_df$mu - qnorm(0.975) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.upper.99 = result$Slist[[index]]$BAUs_prediction_df$mu + qnorm(0.995) * result$Slist[[index]]$BAUs_prediction_df$sd,
      Y.pred.lower.99 = result$Slist[[index]]$BAUs_prediction_df$mu - qnorm(0.995) * result$Slist[[index]]$BAUs_prediction_df$sd
      ),
    proj4string = CRS(paste("EPSG:", st_crs(data_air_dataset)$epsg)))


  comp.comuni.conc = comuni
  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper.95", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper.95"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower.95", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower.95"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.upper.99", data_air_dataset)
  comp.comuni.conc[,"Y.pred.upper.99"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  originalbasis.municipal.values <- GetAvg(prediction, comuni, "Y.pred.lower.99", data_air_dataset)
  comp.comuni.conc[,"Y.pred.lower.99"] = ilrInv(matrix(originalbasis.municipal.values, ncol = 1))

  list(prediction.baus = prediction, prediction.comuni = comp.comuni.conc)
}

