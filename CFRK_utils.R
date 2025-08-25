library(magrittr)
library(compositions)

## algebric methods

InterpThisVector <- function(nodes.coords, nodes.v, eval.coords, extrap = FALSE){
  library(fields)
  library(akima)
  # Interpolazione lineare su griglia regolare
  interp_res <- interp(x = nodes.coords[,1], y = nodes.coords[,2], z = nodes.v,
                      linear = TRUE, extrap = extrap)
  # Interpolazione nei punti arbitrari
  valori_interp <- interp.surface(interp_res, eval.coords)
  valori_interp
}

## PCA Transform

VarsToScores <- function(dataframe, colnames) {
  pca_output <- prcomp(as.data.frame(dataframe)[,colnames])
  list(scores=pca_output$x, pca_output=pca_output)
}


ScoresToVars <- function(dataframe_scores, pca_result, colnames) {
  original_data <- as.data.frame(dataframe_scores %*% t(pca_result$rotation) +
                                   matrix(rep(pca_result$center,
                                              nrow(dataframe_scores)), nrow = nrow(dataframe_scores), byrow = TRUE)) %>% set_colnames(colnames)
  original_data
}


## ILR Transform

TransformILR <- function(df_cx, cnames, df_stations, domain) {
  colnames <- sapply(1:(length(cnames) - 1), FUN=function(i){paste0("Y", as.character(i))})
  df_comp = ilr(df_cx[,cnames]) %>% set_colnames(colnames)

  df1 <- df_comp %>% as.data.frame %>% tibble::rownames_to_column("id")
  df2 <- df_stations %>% tibble::rownames_to_column("id")
  df_comp <- left_join(df1, df2, by = "id")

  df_coord <- st_as_sf(df_comp)
  st_crs(df_coord) <- st_crs(df_stations)
  df_coord <- st_transform(df_coord, st_crs(domain))

  # aggiungo colonne di inversione (normalizzate per somma)
  df_comp <- as.data.frame(setNames(ilrInv(ilr(df_cx[,cnames])), paste0(cnames, "_perc"))) %>% cbind(df_comp)
  #coordinates(df_comp) <- c('Longitude','Latitude') # mimic SpatialPointsDataFrame
  df_comp_out = SpatialPointsDataFrame(data = as.data.frame(df_comp),
    coords = st_coordinates(df_coord), proj4string= CRS(paste("EPSG:", st_crs(domain)$epsg)))
  coordnames(df_comp_out) <- c('x', 'y')
  df_comp_out
}


RevertILR <- function(dataframe, x, y) {
 as.data.frame(matrix(as.numeric(ilrInv(dataframe)), ncol=dim(dataframe)[[2]]+1)) %>%
    set_colnames(paste("c",1:(dim(dataframe)[[2]]+1), sep='')) %>%
    cbind(x=x, y=y)
}

## Logit Transform

logit <- function(p) {
  log(p / (1 - p))
}

inverse_logit <- function(x) {
  1 / (1 + exp(-x))
}

ilr_two_comp <- function(x1, x2) {
  sqrt(1/2) * log(sqrt(x1*x2)/x2)
}


TransformLogit <- function(df_cx, cnames, df_stations, domain, data_p) {
  #colnames <- sapply(1:(length(cnames) - 1), FUN=function(i){paste0("Y", as.character(i))})
  #df_comp = ilr(df_cx[,cnames]) %>% set_colnames(colnames)
  df_comp <- cbind(data_p, logit=log(data_p$p / (1 - data_p$p)))

  df1 <- df_comp %>% as.data.frame %>% tibble::rownames_to_column("id")
  df2 <- df_stations %>% tibble::rownames_to_column("id")
  df_comp <- left_join(df1, df2, by = "id")

  df_coord <- st_as_sf(df_comp)
  st_crs(df_coord) <- st_crs(df_stations)
  df_coord <- st_transform(df_coord, st_crs(domain))

  # aggiungo colonne di inversione (normalizzate per somma)
  #df_comp <- data.frame(logit=log(p / (1 - p))) %>% cbind(df_comp)
  #coordinates(df_comp) <- c('Longitude','Latitude') # mimic SpatialPointsDataFrame
  df_comp_out = SpatialPointsDataFrame(data = as.data.frame(df_comp),
    coords = st_coordinates(df_coord), proj4string= CRS(paste("EPSG:", st_crs(domain)$epsg)))
  coordnames(df_comp_out) <- c('x', 'y')
  df_comp_out
}


RevertLogit <- function(logit, x, y) {
 as.data.frame(cbind(1-inverse_logit(logit), inverse_logit(logit))) %>%
    set_colnames(c("c1", "c2")) %>%
    cbind(x=x, y=y)
}

## FRK objects

GetGrid <- function(domain.geometry, grid.xnum=80, grid.ynum = 80) {
  library(sp)
  sfbbox = GetBBox_FromGeometry(domain.geometry)
  grid = GetSpatialGrid_FromBBox(sfbbox, xnum=grid.xnum, ynum=grid.ynum)
  grid.df = data.frame(grid)
  grid = SpatialPoints(
    grid.df[unlist(st_contains(domain.geometry, st_as_sf(x=grid.df, coords=c('x','y'), crs = st_crs(domain.geometry)))),],
    proj4string = CRS(paste("EPSG:", st_crs(sfbbox)$epsg)))
  grid
}

GetGridBAUs <- function(dataframe, cellsize) {
  set.seed(1)
  # 1. set of BAUs
  GridBAUs1 <- auto_BAUs(manifold = plane(), # 2D plane
                          cellsize = c(cellsize, cellsize), # BAU cellsize
                          type = "grid", # grid (not hex)
                          data = dataframe, # data around which to create BAUs
                          convex = -0.05, # border buffer factor
                          nonconvex_hull = FALSE) # convex hull


  # 2. BAU-specific variation
  # variation at BAU level for hetereoscedascity:
  # for the ith BAU, we also need to attribute the element vi that describes
  # the hetereoscedascity of the fine-scale variation for that BAU
  # (can be tipically altitude)
  GridBAUs1$fs <- 1 # fine-scale variation at BAU level
  GridBAUs1
}

GetFRKFolderOutput <- function(
  figures_folder, DestName.Pollutant,
  Global_Use_Cov_Population, Global_Use_Cov_Elevation,
  frk.args, flag_use_pca, flag_obs_fs, baus.cellsize, DestName.IndexDrop) {
    folder_output = paste0(figures_folder, "/FRK_Pollutant_", DestName.Pollutant,
    "_Global_Use_Cov_Population_", Global_Use_Cov_Population,  "_Global_Use_Cov_Elevation_", Global_Use_Cov_Elevation,
    "_prediction_levels", as.character(frk.args["basis_nres"]) , "_", frk.args["basis_type"],
    "_PCA_", as.character(flag_use_pca), "_obsFs_", as.character(flag_obs_fs), "_BausCellsize_", baus.cellsize,
    "_IndexDrop_", as.character(DestName.IndexDrop))
  folder_output
}

GetBasis <- function(dataframe, frk.args) {
  if ("basis_type" %in% names(frk.args)) {
    type = frk.args[["basis_type"]]
  } else {type = "Gaussian"}
  if ("basis_nres" %in% names(frk.args)) {
    nres = frk.args[["basis_nres"]]
  } else {nres = 3}
  cat(">> AUTO BASIS: using Type = ", type,'\n')
  cat(">> AUTO BASIS: using Nres = ", nres,'\n')

  # 3. basis functions for covariance modeling
  G <- auto_basis(manifold = plane(), # 2D plane
                  data = dataframe, # meuse data
                  nres = nres, # number of resolutions
                  type = type, # type of basis function
                  regular = 0) # place regularly in domain
  G
}

## core FRK functions

FitAndPredict_withFormula <- function(formula, dataframe, GridBAUs, BasisG, pars.fit.no_iterations, tol=0.01) {
    library(FRK) # for carrying out FRK
    S <- SRE(f = as.formula(formula), # formula
                data = list(dataframe), # list of datasets
                #method="EM",
                BAUs = GridBAUs, # BAUs
                basis = BasisG, # basis functions
                est_error = TRUE, # estimation measurement error
                average_in_BAU = FALSE) # do not average data over BAUs

    # EM method
    S <- SRE.fit(S, # SRE model
                    n_EM = pars.fit.no_iterations, # max. no. of EM iterations
                    tol = tol, # tolerance at which EM is assumed to have converged
                    print_lik=TRUE) # print log-likelihood at each iteration

    BAUs_prediction <- predict(S, obs_fs = flag_obs_fs)
    BAUs_prediction_df <- as(BAUs_prediction,"data.frame")
    list(BAUs_prediction_df=BAUs_prediction_df, S=S)
}

FitAndPredict_withFormula_K <- function(formula, dataframe, GridBAUs, BasisG, pars.fit.no_iterations, pars.fit.Ktype) {
    cat("K_type ", pars.fit.Ktype, "\n")
    library(FRK) # for carrying out FRK
    S <- SRE(f = as.formula(formula), # formula
                data = list(dataframe), # list of datasets
                #method="EM",
                BAUs = GridBAUs, # BAUs
                basis = BasisG, # basis functions
                est_error = TRUE, # estimation measurement error
                average_in_BAU = FALSE) # do not average data over BAUs

    # EM method
    S <- SRE.fit(S, # SRE model
                    n_EM = pars.fit.no_iterations, # max. no. of EM iterations
                    tol = 0.01, # tolerance at which EM is assumed to have converged
                    print_lik=TRUE) # print log-likelihood at each iteration

    BAUs_prediction <- predict(S, obs_fs = flag_obs_fs)
    BAUs_prediction_df <- as(BAUs_prediction,"data.frame")
    list(BAUs_prediction_df=BAUs_prediction_df, S=S)
}

## FRK steps functions

PrepareData_For_Kriging <- function(data.comp, colnames.euclid, flag_use_pca) {
  if (flag_use_pca)
  {
    cat('YES pca\n')
    pca_result <- VarsToScores(as.data.frame(data.comp)[,colnames.euclid], colnames.euclid)
    pca_output <- pca_result$pca_output

    data.krig = SpatialPointsDataFrame(
      coordinates(data.comp),
      as.data.frame(pca_result$scores) %>%
        set_colnames(colnames.euclid),
      proj4string = CRS(paste("EPSG:",st_crs(data.comp)$epsg)))
  } else {
    cat('NO pca\n')
    data.for.spdf = as.data.frame(as.data.frame(data.comp)[, colnames.euclid])
    colnames(data.for.spdf) <- colnames.euclid
    data.krig = SpatialPointsDataFrame(
      coordinates(data.comp),
      data.for.spdf,
      proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)), match.ID = FALSE)
    pca_output = {}
  }
  list(data.krig = data.krig, pca_output = pca_output)
}

DoFRKPrediction <- function(domain, GridBAUs, frk.args, data.comp, data.krig,
  colnames.euclid, comuni,
  Global_Use_Cov_Population, Global_Use_Cov_Elevation, no_iterations,   fitted.locations = NULL, tol=0.01) {
  colnames_pred = paste0(colnames.euclid, '.pred')

  if (is.null(GridBAUs)) {
    prediction.grid = data_air$grid
    GridBAUs1 <- GetGridBAUs(prediction.grid, cellsize=baus.cellsize)

    if (is.null(cov.P1)) {
      # aggiungo popolazione a GridBaus
      layer_popolazione.Spat = as(st_transform(layer_popolazione, st_crs(GridBAUs1)), "Spatial")
      cov.P1 = over(SpatialPoints(GridBAUs1, proj4string = CRS("+proj=utm +zone=32 +datum=WGS84 +units=m +no_defs")), layer_popolazione.Spat)$Popolazione.P1
    }

    grid_r <- st_transform(st_as_sf(as.data.frame(coordinates(GridBAUs1)),
      coords = c('x', 'y'), crs = st_crs(prediction.grid)),
                crs = st_crs(raster.elevazione))
    # conversione covariata per metodo
    cov.Elev = raster::extract(raster.elevazione,
      grid_r)

    cat("cov.P1", length(cov.P1), "\n")
    cat("cov.Elev", length(cov.Elev), "\n")
    discard_tiles = is.na(cov.P1) | is.na(cov.Elev)
    cat("dim(GridBAUs1)", dim(GridBAUs1[!discard_tiles,]), "\n")

    GridBAUs = GridBAUs1[!discard_tiles,]

    GridBAUs$Elevazione = cov.Elev[!discard_tiles]
    GridBAUs$Popolazione.P1 = cov.P1[!discard_tiles]
  }

  if (!FLAG_NO_PLOT) {
    coogrid = as.data.frame(coordinates(GridBAUs1))[ ,c('x', 'y')]
    ggplot() + geom_sf(data=domain, aes()) + geom_point(data=coogrid, aes(x, y))
    ggsave(paste(folder_output, "baus.png", sep="/"))
  }

  #spatial_domain <- as(domain, "Spatial")
  domain.prediction.grid <- GetGrid(domain) # TODO construct GetGridBAUs without domain.prediction.gridv
  G <- GetBasis(domain.prediction.grid, frk.args)
  cat("done basis G\n")

  # ---------------------------------------------------------------------------
  # ONLY SHOW
  if (!FLAG_NO_PLOT) {
    g <- ggplot() + geom_sf(data=domain, aes())
    show_basis(G, g = g) + # coord_fixed() + # fix aspect ratio
      xlab("Easting (m)") + # x-label
      ylab("Northing (m)") # y-label
    ggsave(paste(folder_output, "basis_functions.png", sep="/"))
  }

  # ONLY OUTPUT
  cat("INFO: BASE\n")
  cat("Numero Funzioni Base = ", attributes(G)$n, "\n\n")
  cat("Numero osservazioni = ", dim(data.comp)[1], "\n\n")
  cat('\n')


  # ============================================================================
  # 1   FIT AND PREDICTION
  # 1.1 Prediction su BAUS

  # componente 1 per inizializzare
  i = 1
  formula = paste0(colnames.euclid[i], " ~ 1")
  result <- FitAndPredict_withFormula(formula, data.krig, GridBAUs, G, pars.fit.no_iterations=no_iterations)


  prediction = SpatialPixelsDataFrame(
    as.matrix(as.data.frame(result$BAUs_prediction_df)[,c('x','y')]),
    data = data.frame(Y1.pred.dummy=result$BAUs_prediction_df$mu),
    proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)))

  Slist = list()
  # tutte le predizioni
  for (i in (1:length(colnames.euclid)))
  {
    formula = paste0(colnames.euclid[i], " ~ 1")
    if(Global_Use_Cov_Population)
    {
      formula = paste0(colnames.euclid[i], " ~ Popolazione.P1")
    }
    if((Global_Use_Cov_Population) & (Global_Use_Cov_Elevation))
    {
      formula = paste0(colnames.euclid[i], " ~ Popolazione.P1 + Elevazione")
    }
    cat("using formula in frk: ", formula, "\n")
    cat("using tolerance in frk: ", tol, "\n")

    if(!FLAG_NO_PLOT) {
      x11(width = 8, height = 8)
      png(paste0(folder_output, '/convergence_loglik_index_', i, '_Global_Use_Cov_Population_', Global_Use_Cov_Population, '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '.png'))
    }

    if (exists("Flag_K_type_exponential") && Flag_K_type_exponential) {
      result <- FitAndPredict_withFormula_K(formula, data.krig, GridBAUs, G, pars.fit.no_iterations=no_iterations, pars.fit.Ktype = "block-exponential", tol=tol)
    } else {
      result <- FitAndPredict_withFormula(formula, data.krig, GridBAUs, G, pars.fit.no_iterations=no_iterations, tol=tol)
    }
    Slist = append(Slist, list(result))
    result.fit.pred = result
    dev.off()
    prediction[[paste0(colnames.euclid[i],".pred")]] = result$BAUs_prediction_df$mu
    prediction[[paste0(colnames.euclid[i],".sd")]] = result$BAUs_prediction_df$sd
  }

  GetAvg = function(dataframe.prediction, comuni, field) {
    library(sp)
    library(sf)
    library(raster)
    spdf = SpatialPixelsDataFrame(coordinates(dataframe.prediction)[, c('x', 'y')], data = as.data.frame(dataframe.prediction[[field]]))
    raster_obj <- raster(spdf)
    crs(raster_obj) <- CRS(paste("EPSG:", st_crs(data.comp)$epsg))
    library(exactextractr)
    exact_extract(raster_obj, comuni, fun="mean")
  }

  geometries = st_transform(comuni$geometry, st_crs(GridBAUs1))

  comuniresultall.meaneucl = comuni

  for (i in (1:length(colnames_pred)))
  {
    datafill = GetAvg(prediction, comuni, colnames_pred[i])
    comuniresultall.meaneucl[colnames_pred[i]] = datafill
  }


  # 1.3 Prediction on observed locations (fitted)
  Extract <- function(prediction, fitted.locations, field) {
    library(raster)
    raster::extract(raster(SpatialPixelsDataFrame(
      coordinates(prediction),
      data = as.data.frame(prediction[[field]]))), fitted.locations )
  }


  if (is.null(fitted.locations)) {
    fitted.locations <- SpatialPoints(as.data.frame(data.comp)[c("x", "y")], proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)))
  }


  fitted.values = sapply(colnames_pred, FUN=function(field){Extract(prediction, fitted.locations, field)})
  fitted.values = fitted.values %>% matrix(nrow = length(fitted.locations)) %>% data.frame %>% set_colnames(colnames_pred)
  fitted.values = SpatialPointsDataFrame(coordinates(fitted.locations), data = fitted.values, proj4string = CRS(paste("EPSG:", st_crs(data.comp)$epsg)), match.ID = FALSE)

  list(
    predictionbaus.values = prediction,
    fitted.values = fitted.values,
    comuni.values = comuniresultall.meaneucl,
    Slist=Slist)
}

## Plot

DoFRKPlotPrediction <- function(prediction.conc, comp.fitted.df, comp.comuni.conc,
colnames.aitch, colnames.aitch.pred, DoPLOT) {
  # numbers
  prediction.conc.n = prediction.conc[, paste0("c", as.character(1:length(colnames.aitch)))] %>% set_colnames(colnames.aitch.pred)
  prediction.conc.n$x = prediction.conc$x
  prediction.conc.n$y = prediction.conc$y

  comp.fitted.df.n = comp.fitted.df[, paste0("c", as.character(1:length(colnames.aitch)))] %>% set_colnames(colnames.aitch.pred)
  comp.fitted.df.n$x = comp.fitted.df$x
  comp.fitted.df.n$y = comp.fitted.df$y

  if (DoPLOT) # PLOT su  BAUS
  {
    cat("Printing maps on BAUS\n")
    library(scales)  # Per la funzione label_number()


    x11(width = 8, height = 7)
    for (i in (1:length(colnames.aitch.pred)))
    {
      datafill = prediction.conc.n[[colnames.aitch.pred[i]]]
      datatitle = paste0("Prediction of interval I", as.character(i))
      ggpred =  ggplot() + geom_sf(data=domain$geom, aes()) +
        geom_tile(data = prediction.conc.n, aes(x,y,fill=get(colnames.aitch.pred[i]))) + ggtitle(datatitle) +
        guides(fill=guide_legend(title="")) +
        theme(plot.margin = unit(c(0.8, 0.5, 0.5, 1), "cm"),
              plot.title = element_text(size = 20, hjust=0.5),
              legend.text = element_text(size = 12),   # Dimensione delle etichette della legenda
              axis.title.x = element_blank(),  # Nascondere l'etichetta dell'asse X
              axis.title.y = element_blank(),   # Nascondere l'etichetta dell'asse Y
              #legend.title = element_text(size = 14),   # Dimensione del titolo della legenda
              #axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
              #axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
              axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
              axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
        ) +
        scale_fill_viridis_c(option = "C", direction = 1,
                             breaks = seq(min(datafill), max(datafill),
                            length.out = 10), labels = label_number(0.0001)) # legend name
      ggpred
      ggsave(paste(folder_output, paste0("prediction_cofrk_bausnopoints_component_", as.character(i),
          '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
          '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
      obs_name = paste0("R", as.character(i), "_perc")
      ggpred = ggpred + geom_point(data = data.frame(data.comp), aes(x,y, fill=get(obs_name)), colour="black", pch=21, size=4)
      ggpred + geom_text(aes(x = data.comp$x, y = data.comp$y, label = data.comp$stationid), vjust = -1, size=2)
      ggsave(paste(folder_output, paste0("prediction_cofrk_baus_component_", as.character(i),
          '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
          '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
    }
    graphics.off()
  }

  if (TRUE) # PLOT su base COMUNALE
  {
    for (i in (1:length(colnames.aitch.pred)))
    {
      datafill = as.numeric(as.data.frame(comp.comuni.conc)[,colnames.aitch.pred[i]])
      comuniresult = comuni
      comuniresult$datafill = datafill
      if (DoPLOT) {
        x11(width = 8, height = 7)
        breaks = seq(min(datafill, na.rm = T), max(datafill, na.rm = T), length.out = 10)
        limits = c(min(datafill, na.rm = T), max(datafill, na.rm = T))
        if (Use_Global_Colors)
        {
          colors.pred.breaks = Colors_Breaks
          breaks = seq(colors.pred.breaks[[i]][1], colors.pred.breaks[[i]][2], length.out = 10)
          limits = c(colors.pred.breaks[[i]][1], colors.pred.breaks[[i]][2])
        }
        datatitle = paste0("CFRK prediction of relative frequencies \n in interval I", as.character(i),
                          " = ",  data_air$levels[i], " μg/m3 on municipal base")
        gcomuni = ggplot() + geom_sf(data = comp.comuni.conc$geometry, aes(fill=datafill), linewidth = 0.1) + ggtitle(datatitle) +
          guides(fill=guide_legend(title="Frequency")) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 1), "cm"),
                plot.title = element_text(size = 20, hjust=0.5),
                legend.text = element_text(size = 14),   # Dimensione delle etichette della legenda
                legend.title = element_text(size = 14),   # Dimensione del titolo della legenda
                axis.title.x = element_text(size = 15),  # Dimensione del titolo dell'asse x
                axis.title.y = element_text(size = 15),  # Dimensione del titolo dell'asse y
                axis.text.x = element_text(size = 12),   # Dimensione dei valori sull'asse x
                axis.text.y = element_text(size = 12)    # Dimensione dei valori sull'asse y
          ) + scale_fill_viridis_c(option = "C", direction = 1,
                                  breaks = breaks,
                                  limits = limits,
                                  labels = label_number(0.0001)) +# legend name
          labs(fill = "", x = "Longitude", y = "Latitude")

        gcomuni

        ggsave(paste(folder_output, paste0("prediction_cofrk_comunale_component_", as.character(i),
                                            '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
                                            '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, ".png"), sep="/"))
      }
      st_write(comuniresult, paste0(folder_output,'/prediction_cofrk_comunale_component_', i,
                                    '_Global_Use_Cov_Population_', Global_Use_Cov_Population,
                                    '_Global_Use_Cov_Elevation_', Global_Use_Cov_Elevation, '.gpkg'), append = FALSE)
    }
    graphics.off()
  }

}