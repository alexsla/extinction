library(tidyverse)
library(rnaturalearth)
library(sf)
library(terra)
library(biomod2)
library(blockCV)
library(parallel)
library(pbapply)

unlink("occ.bird", recursive = TRUE)

sf_use_s2(FALSE)

cores = detectCores()
cl <- 20 ## set to use all but one thread - replace if necessary

# import Australian bird shapefiles
bird_dist <- readRDS("data/bird_dist.rds")

# create map of Australia
oz_map <-
  ne_countries(country = "Australia",
               scale = 10,
               returnclass = "sf") %>%
  st_transform(st_crs(bird_dist)) %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.01)

# import point locality of Australian bird occurrences
bird_points <- readRDS("data/bird_points.RDS") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = st_crs(oz_map))

# generate vector of bird taxa
bird_tax <- unique(bird_points$Binomial)
bird_tax <- bird_tax[1:length(bird_tax)]

# import current environmental layers
currentEnv <- rast("data/SDM/Current.tif")
names(currentEnv) <-
  c(
    "bio1",
    "bio4",
    "bio10",
    "bio11",
    "bio12",
    "bio15",
    "bio16",
    "bio17",
    "elev",
    "AWC",
    "BDW",
    "CLY",
    "pHc",
    "SLT",
    "SND",
    "SOC"
  )

# import future environmental layers
futureLayers <- dir("data/SDM", full.names = T)
futureLayers <-
  futureLayers[is.na(str_match(futureLayers, "Current"))]
futureEnv <- list()
for (i in 1:length(futureLayers)) {
  futureEnv[[i]] <- rast(futureLayers[[i]])
  names(futureEnv[[i]]) <-
    c(
      "bio1",
      "bio4",
      "bio10",
      "bio11",
      "bio12",
      "bio15",
      "bio16",
      "bio17",
      "elev",
      "AWC",
      "BDW",
      "CLY",
      "pHc",
      "SLT",
      "SND",
      "SOC"
    )
}

# generate vectors with names of models and years of projection
models <-
  c(
    "SSP1.26",
    "SSP1.26",
    "SSP1.26",
    "SSP1.26",
    "SSP5.85",
    "SSP5.85",
    "SSP5.85",
    "SSP5.85"
  )
years <- c(
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100"
)

# generate background points
# set seed for reproducibility
set.seed(42)
points_bg_bird <-
  bird_points %>% sample_frac(0.01)

points_bg_bird <-
  dismo::gridSample(st_coordinates(points_bg_bird), currentEnv, n = 1)

oz_grid <- st_make_grid(oz_map, cellsize = res(currentEnv) * 10)

# generate new output files if don't already exist
if(!file.exists("output/SDM/bird_summary.csv")){
  tibble(
    AUC = numeric(),
    BOYCE = numeric(),
    TSS = numeric(),
    block_size = numeric(),
    k = numeric(),
    threshold = numeric(),
    n = numeric(),
    Binomial = character()
  ) %>%
    write_csv("output/SDM/bird_summary.csv")
}

if(!file.exists("output/SDM/bird_failed.csv")){
  tibble(
    Binomial = character(),
    reason = character()
  ) %>%
    write_csv("output/SDM/bird_failed.csv")
}

# iterate loop over all bird taxa
for (i in 1:length(bird_tax)) {
  if(bird_tax[[i]] %in% read_csv("output/SDM/bird_summary.csv")$Binomial)
    next
  
  # cleanup
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  
  # filter points for ith species and generate background points
  points_occ <-
    st_coordinates(st_as_sf(
      bird_points %>% dplyr::filter(Binomial == bird_tax[[i]]),
      crs = st_crs(oz_map)
    )) %>%
    as_tibble()
  
  if (nrow(points_occ) < 10) {
    tibble(Binomial = bird_tax[[i]],
           reason = "not enough points") %>%
      write_csv("output/SDM/bird_failed.csv", append = T)
    tibble(AUC = NA,
           BOYCE = NA,
           TSS = NA,
           block_size = NA,
           k = NA,
           threshold = NA,
           n = nrow(points_occ),
           Binomial = bird_tax[[i]]) %>%
      write_csv("output/SDM/bird_summary.csv", append = T)
    next
  }
  
  pnts_occ <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_occ),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  # drop from background points all presence points
  pnts_bg <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_bg_bird),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  pnts_logical = lengths(pnts_bg) > 0 & lengths(pnts_occ) == 0
  
  oz_bg <- st_as_sf(oz_grid, "POLYGON")[pnts_logical,]
  
  points_bg <-
    points_bg_bird[lengths(st_intersects(st_as_sf(
      as.data.frame(points_bg_bird),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ), oz_bg)) > 0, ]
  
  sp_training <-
    bind_rows(as_tibble(points_occ) %>%
                mutate(occ = 1),
              as_tibble(points_bg) %>%
                mutate(occ = 0)) %>%
    relocate(occ) %>%
    st_as_sf(coords = c("X", "Y"), crs = st_crs(oz_map))
  
  sp_training <-
    sp_training %>%
    bind_cols(extract(currentEnv, sp_training)[,-1]) %>%
    filter_all(~!is.na(.))
  
  training <-
    sp_training %>%
    as_tibble() %>%
    dplyr::select(-geometry)
  
  biomod_data <- BIOMOD_FormatingData(resp.var = training$occ,
                                      expl.var = training[,-1],
                                      resp.xy = sf::st_coordinates(sp_training),
                                      resp.name = "occ.bird",
                                      na.rm = TRUE)
  
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  wt <- ifelse(training$occ == 1, 1, prNum / bgNum)
  
  # Try cv_spatial_autocor, but continue if it fails
  cv_test_result <- tryCatch({
    cv_spatial_autocor(r = currentEnv, x = sp_training, column = "occ")
  }, error = function(e) {
    cat("cv_spatial_autocor failed, using default size\n")
    return(NULL)
  })
  
  # Set size based on whether cv_spatial_autocor succeeded
  if (is.null(cv_test_result)) {
    block_size <- bird_dist %>% filter(Binomial == bird_tax[[i]]) %>% st_area() %>% as.numeric() %>% sqrt() %>% round(0)
  } else {
    min_area <- bird_dist %>% filter(Binomial == bird_tax[[i]]) %>% st_area() %>% as.numeric() %>% sqrt() %>% round(0)
    block_size <- min(ceiling(cv_test_result$range/1000)*1000, min_area)
  }
  
  # generate CV blocks
  # Start with k = 4
  k_value <- 4
  success <- FALSE
  warning_messages <- character(0)
  
  # Create a function to run cv_spatial and capture warnings
  run_cv_spatial <- function(k) {
    # This will store any warnings
    warnings_captured <- NULL
    
    # Run cv_spatial with warning handler
    result <- withCallingHandlers(
      tryCatch({
        cv_spatial(
          x = sp_training,
          column = "occ",
          r = currentEnv,
          k = k, # number of folds
          size = block_size, # size of blocks in metres
          selection = "random", # random blocks-to-fold
          iteration = 50, # find evenly dispersed folds
          progress = FALSE, # turn off progress bar
          biomod2 = TRUE, # also create folds for biomod2
          seed = 42
        )
      }, error = function(e) {
        # Return the error message if there's an error
        return(list(result = NULL, error = conditionMessage(e)))
      }),
      warning = function(w) {
        # Capture warnings
        warnings_captured <<- c(warnings_captured, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    
    # Return both the result and any warnings
    return(list(result = result, warnings = warnings_captured))
  }
  
  # Loop until successful or k reaches 2
  while (!success && k_value >= 2) {
    cat(paste0("Trying with k = ", k_value, "\n"))
    
    # Run cv_spatial and capture result and warnings
    output <- run_cv_spatial(k_value)
    scv1 <- output$result
    current_warnings <- output$warnings
    
    # Check if there was a zero records warning
    zero_records_warning <- grep("Fold.* .* ha(s|ve) class\\(es\\) with zero records", current_warnings)
    
    if (length(zero_records_warning) == 0) {
      # No zero records warning, we're successful
      success <- TRUE
      cat("Successfully created spatial CV blocks.\n")
    } else {
      # Got zero records warning, reduce k and try again
      warning_messages <- c(warning_messages, 
                            paste0("With k = ", k_value, ": ", 
                                   current_warnings[zero_records_warning[1]]))
      k_value <- k_value - 1
    }
  }
  
  # Output results
  if (success) {
    cat(paste0("Final k value used: ", k_value, "\n"))
  } else {
    cat("Failed to create spatial CV blocks even with k = 2.\n")
    tibble(Binomial = bird_tax[[i]],
           reason = "failed to create spatial CV blocks") %>%
      write_csv("output/SDM/bird_failed.csv", append = T)
    tibble(AUC = NA,
           BOYCE = NA,
           TSS = NA,
           block_size = NA,
           k = NA,
           threshold = NA,
           n = nrow(points_occ),
           Binomial = bird_tax[[i]]) %>%
      write_csv("output/SDM/bird_summary.csv", append = T)
    next
  }
  
  # The final scv1 object contains the successful spatial cross-validation result
  
  # use generated folds from cv_spatial in previous section
  spatial_cv_folds <- scv1$biomod_table
  # the new update of the package biomod2 (v4.2-3 <) requires the names to be as below
  colnames(spatial_cv_folds) <- paste0("_allData_RUN", 1:ncol(spatial_cv_folds))
  
  # generate SDMs (use bigboss settings)
  biomod_model_out <- tryCatch({
    BIOMOD_Modeling(biomod_data,
                    models = c('GBM','GLM','GAM','MAXENT','RF'),
                    OPT.strategy = 'bigboss',
                    CV.strategy = 'user.defined',
                    CV.user.table = spatial_cv_folds,
                    metric.eval = 'TSS',
                    weights = wt,
                    seed.val = 42,
                    nb.cpu = cl)
  }, error = function(e) {
    cat("Failed to fit SDMs.\n")
    tibble(Binomial = bird_tax[[i]],
           reason = "failed to fit SDMs") %>%
      write_csv("output/SDM/bird_failed.csv", append = T)
    tibble(AUC = NA,
           BOYCE = NA,
           TSS = NA,
           block_size = NA,
           k = NA,
           threshold = NA,
           n = nrow(points_occ),
           Binomial = bird_tax[[i]]) %>%
      write_csv("output/SDM/bird_summary.csv", append = T)
    return(NULL)
  })
  
  # Check if modeling failed and skip to next iteration
  if (is.null(biomod_model_out)) {
    unlink("occ.bird", recursive = TRUE)
    next
  }
  
  # generate ensemble model
  ens_model <- BIOMOD_EnsembleModeling(biomod_model_out,
                                       em.by = "all",
                                       em.algo = 'EMwmean',
                                       metric.eval = 'TSS',
                                       metric.select.thresh = 0)

  if(any(ens_model@em.models_kept == "none")) {
    cat("All TSS values below 0.\n")
    tibble(Binomial = bird_tax[[i]],
           reason = "All TSS values below 0") %>%
      write_csv("output/SDM/bird_failed.csv", append = T)
    tibble(AUC = NA,
           BOYCE = NA,
           TSS = NA,
           block_size = NA,
           k = NA,
           threshold = NA,
           n = nrow(points_occ),
           Binomial = bird_tax[[i]]) %>%
      write_csv("output/SDM/bird_summary.csv", append = T)
    unlink("occ.bird", recursive = TRUE)
    next
  }
  
  # plot partial response curves
  ggsave(paste("output/SDM/Response Curves/", bird_tax[[i]], ".pdf", sep = ""),
         bm_PlotResponseCurves(ens_model, do.plot = F)$plot,
         width = 7,
         height = 7,
         units = "in")
  
  # projections
  mod.dir <- dir("occ.bird/models", full.names = T)
  mods <- lapply(dir(mod.dir, full.names = T)[!str_detect(dir(mod.dir), "outputs|merged")],
                 function(x)
                   get(load(x)))
  names(mods) <- dir(mod.dir)[!str_detect(dir(mod.dir), "outputs|merged")]
  
  mods <- mods[names(mods) %in% ens_model@em.models_kept]
  
  weights <- sapply(1:length(mods), function(x) mods[[x]]@model_evaluation[,6])
  weights <- weights/sum(weights, na.rm = TRUE)
  
  # Create temporary directory for parallel predictions
  temp_pred_dir <- tempdir()
  temp_pred_files <- file.path(temp_pred_dir, paste0("pred_", 1:length(mods), ".tif"))
  
  # Use mclapply with disk-based storage to avoid pointer issues
  mclapply(1:length(mods),
           function(x){
             temp_workdir = NULL
             if(str_detect(names(mods)[[x]], "MAXENT")){
               temp_workdir = mods[[x]]@model_output_dir
             }
             pred <- predict(mods[[x]], currentEnv, on_0_1000 = T,
                             temp_workdir = temp_workdir, overwrite = T,
                             seedval = 42, mod.name = names(mods)[[x]])
             # Write prediction to disk
             terra::writeRaster(pred, temp_pred_files[x], overwrite = TRUE)
             return(temp_pred_files[x])
           },
           mc.cores = cl
  )
  
  # Read predictions back from disk
  ens_pred <- lapply(temp_pred_files, function(f) terra::rast(f))
  pr <<- weighted.mean(do.call(c, ens_pred), w = weights)
  names(pr) = "layer"
  
  pr_df = as.data.frame(pr/1000, xy = T)
  
  # Clean up temporary files
  file.remove(temp_pred_files)
  
  # plot habitat suitability probability
  ggsave(paste("output/SDM/Current/Maps/", bird_tax[[i]], "_prob.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_raster(data = pr_df %>% drop_na(), aes(
             x = x, y = y, fill = layer
           )) +
           theme_bw() +
           scale_fill_gradientn(colours = viridis::viridis(99),
                                name = "Predicted Habitat Suitability") +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  # extract model estimated suitability
  est = terra::extract(pr, sp_training)[,2]
  # detect threshold for presence/absence
  thr = bm_FindOptimStat(metric.eval = "TSS", obs = sp_training$occ, fit = est)
  # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
  pr_thr = pr > thr$cutoff
  
  # convert current presence-absence map to polygon
  current <-
    as.polygons(pr_thr) %>%
    st_as_sf() %>%
    filter(layer == 1) %>%
    st_cast("POLYGON")
  
  # find which areas of suitable habitat do not overlap occurrences
  pnts_sf <-
    st_intersects(current,
                  st_as_sf(
                    bird_points %>% dplyr::filter(Binomial == bird_tax[[i]]),
                    coords = c("X", "Y"),
                    crs = st_crs(oz_map)
                  ))
  pnts_logical = lengths(pnts_sf) > 0
  
  # filter out areas of suitable habitat which do not overlap landmasses with occurrences
  current <- current[pnts_logical, ] %>%
    st_union() %>%
    st_as_sf() %>%
    dplyr::mutate(Binomial = bird_tax[[i]])
  
  # create a buffer around current suitable habitat to allow for dispersal
  current_buff <- current %>%
    st_transform(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs") %>%
    st_buffer(dist = 150) %>% # dispersal buffer
    st_transform(crs = st_crs(oz_map)) %>%
    st_make_valid() %>% # fix geometry errors due to buffer
    st_intersection(oz_map) %>%
    st_cast("POLYGON")
  
  # find and filter out areas of buffer that do not occur on the same landmasses
  current_oz <- st_intersects(current_buff, current)
  current_logical = lengths(current_oz) > 0
  
  extant_buff <- current_buff[current_logical, ] %>%
    st_union() %>%
    st_as_sf() %>%
    dplyr::mutate(Binomial = bird_tax[[i]])
  
  # identify landmasses
  landmasses <- st_intersects(oz_map %>%
                                st_cast("POLYGON"),
                              st_as_sf(
                                bird_points %>% dplyr::filter(Binomial == bird_tax[[i]]),
                                coords = c("X", "Y"),
                                crs = st_crs(oz_map)
                              ))
  landmasses <- st_cast(oz_map, "POLYGON")[lengths(landmasses) > 0, ] %>%
    st_union() %>%
    st_as_sf()
  
  # export habitat suitability under current climate conditions
  current_predictions <- c(pr/1000,
                           as.numeric(pr_thr))
  names(current_predictions) <-
    c("predicted_probability", "presence_absence")
  writeRaster(current_predictions,
              paste("output/SDM/Current/Rasters/", bird_tax[[i]], ".tif", sep = ""),
              overwrite = T)
  write_sf(current,
           paste("output/SDM/Current/Shapefiles/", bird_tax[[i]], ".shp", sep = ""))
  
  # plot presence/absence based on threshold under current climate conditions
  ggsave(paste("output/SDM/Current/Maps/", bird_tax[[i]], "_pa.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_sf(
             data = current,
             colour = NA,
             fill = "purple"
           ) +
           geom_sf(
             data = st_as_sf(
               dplyr::filter(bird_points, Binomial == bird_tax[[i]]),
               crs = st_crs(oz_map)
             ),
             size = .1
           ) +
           theme_bw() +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  # make future predictions
  # over all 8 climate projections
  future_predictions <- list()
  for(k in 1:length(futureEnv)) {
    # make predictions under future conditions
    # Create temporary directory for parallel predictions
    temp_pred_dir <- tempdir()
    temp_pred_files <- file.path(temp_pred_dir, paste0("pred_", 1:length(mods), ".tif"))
    
    # Use mclapply with disk-based storage to avoid pointer issues
    mclapply(1:length(mods),
             function(x){
               temp_workdir = NULL
               if(str_detect(names(mods)[[x]], "MAXENT")){
                 temp_workdir = mods[[x]]@model_output_dir
               }
               pred <- predict(mods[[x]], futureEnv[[k]], on_0_1000 = T,
                               temp_workdir = temp_workdir, overwrite = T,
                               seedval = 42, mod.name = names(mods)[[x]])
               # Write prediction to disk
               terra::writeRaster(pred, temp_pred_files[x], overwrite = TRUE)
               return(temp_pred_files[x])
             },
             mc.cores = cl
    )
    
    # Read predictions back from disk
    ens_predf <- lapply(temp_pred_files, function(f) terra::rast(f))
    prf <<- weighted.mean(do.call(c, ens_predf), w = weights)
    names(prf) = "layer"
    
    # Clean up temporary files
    file.remove(temp_pred_files)
    
    # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
    prf_thr = prf > thr$cutoff
    
    # save habitat suitability under future climate conditions
    future_preds <- c(prf/1000,
                      as.numeric(prf_thr))
    names(future_preds) <-
      c("predicted_probability", "presence_absence")
    
    future_predictions[[k]] <- future_preds
  }
  
  names(future_predictions) <-
    str_remove_all(futureLayers, "data/SDM/|.tif")
  
  # extract all future predicted presence-absence maps
  dir.create("temp_rast")
  
  lapply(1:length(future_predictions),
         function(k) writeRaster(future_predictions[[k]][[2]],
                                 filename = paste0("temp_rast/", k, ".tif"),
                                 overwrite = T))
  

  
  future_all <- mclapply(1:3,
                        function(scenario) {
                          
                          future_n <- list()
                          
                          # iterate over all future projections
                          kmax = 4
                          for (m in 1:2) {
                            
                            if(scenario == 1){
                              ## standard dispersal model ##
                              # for each GCM and SSP set the current buffered suitable habitat to be the extant
                              current_buff <- extant_buff
                            }
                            
                            if(scenario == 2){
                              ## standard dispersal model ##
                              # for each GCM and SSP set the current suitable habitat (no buffer) to be the extant
                              current_buff <- current
                            }
                            
                            for (k in (kmax - 3):kmax) {
                              future <- 
                                rast(paste0("temp_rast/", k, ".tif"))
                              
                              # for each time-step
                              if(scenario %in% c(1, 2)){
                                ## standard dispersal model ##
                                # find and filter out areas of suitable habitat which do not intersect with buffered suitable habitat for current time-step
                                future_oz <- mask(future, current_buff)
                              }
                              
                              if(scenario == 3){
                                ## max dispersal ##
                                future_oz <- mask(future, landmasses)
                              }
                              
                              # if no suitable habitat exists in time step, save empty polygon
                              if (!any(unique(future_oz) == 1)){
                                future_next <- st_sf(x = st_sfc(st_multipolygon()),
                                                     crs = st_crs(oz_map)) %>%
                                  dplyr::mutate(Binomial = bird_tax[[i]],
                                                Model = models[[k]],
                                                Year = years[[k]])
                              } else {
                                # import presence-absence map for the next time-step and convert to polygon
                                future_next <- as.polygons(future_oz) %>%
                                  st_as_sf(crs = st_crs(oz_map)) %>%
                                  filter(presence_absence > 0) %>%
                                  st_cast("POLYGON") %>%
                                  st_union() %>%
                                  st_as_sf() %>%
                                  dplyr::mutate(Binomial = bird_tax[[i]],
                                                Model = models[[k]],
                                                Year = years[[k]])
                              }
                              
                              # save filtered suitable habitat
                              if(scenario == 1){
                                ## standard dispersal model ##
                                # update buffered current suitable habitat for next iteration
                                current_buff <- future_next %>%
                                  st_transform(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs") %>%
                                  st_buffer(dist = 150) %>% # numeric buffer of 150 km max dispersal range
                                  st_transform(crs = st_crs(oz_map)) %>%
                                  st_make_valid() %>% # fix geometry errors due to buffer
                                  st_intersection(oz_map) %>%
                                  st_cast("POLYGON")
                                
                                current_oz <- st_intersects(current_buff, future_next)
                                current_logical = lengths(current_oz) > 0
                                
                                current_buff <- current_buff[current_logical, ] %>%
                                  st_union() %>%
                                  st_as_sf() %>%
                                  dplyr::mutate(Binomial = bird_tax[[i]],
                                                Model = models[[k]],
                                                Year = years[[k]])
                                
                                future_n[[k]] <- future_next
                              }
                              
                              
                              if(scenario == 2){
                                ## no dispersal ##
                                future_n[[k]] <- current_buff <- future_next
                              }
                              
                              if(scenario == 3){
                                ## max dispersal ##
                                future_n[[k]] <- future_next
                              }
                            }
                            
                            kmax = kmax + 4
                          }
                          
                          future_n
                        },
                        mc.cores = 3)

  future <- bind_rows(future_all[[1]])
  future_maxdisp <- bind_rows(future_all[[3]])
  future_nodisp <- bind_rows(future_all[[2]])
  
  unlink("temp_rast", recursive = TRUE)
  
  future <- bind_rows(future %>% mutate(Dispersal = "Mean"),
                      future_nodisp %>% mutate(Dispersal = "Min"),
                      future_maxdisp %>% mutate(Dispersal = "Max"))
  
  # export habitat suitability under future climate conditions
  writeRaster(do.call(c, lapply(1:length(future_predictions),
                                function(x)
                                  setNames(future_predictions[[x]],
                                           paste(names(future_predictions)[x], names(future_predictions[[x]]), sep = "_")))),
              paste("output/SDM/Future/Rasters/", bird_tax[[i]], ".tif", sep = ""), overwrite = T)
  write_sf(future,
           paste("output/SDM/Future/Shapefiles/", bird_tax[[i]], ".shp", sep = ""))
  
  # combine to dataframe of predicted probability (for plotting)
  prf_df <- list()
  for (k in 1:length(future_predictions)) {
    prf_df[[k]] <-
      as.data.frame(future_predictions[[k]]$predicted_probability, xy = T) %>%
      mutate(model = models[[k]],
             year = years[[k]])
  }
  prf_df <- bind_rows(prf_df)
  
  # plot all future predicted probability maps
  ggsave(paste("output/SDM/Future/Maps/", bird_tax[[i]], "_prob.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_raster(data = prf_df %>% drop_na(), aes(
             x = x, y = y, fill = predicted_probability
           )) +
           theme_bw() +
           scale_fill_gradientn(colours = viridis::viridis(99),
                                name = "Predicted Habitat Suitability") +
           facet_grid(cols = vars(model),
                      rows = vars(year)) +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  # plot all future presence/absence maps
  ggsave(paste("output/SDM/Future/Maps/", bird_tax[[i]], "_pa_meanDispersal.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_sf(
             data = future %>%
               filter(Dispersal == "Mean"),
             colour = NA,
             fill = "purple"
           ) +
           theme_bw() +
           facet_grid(cols = vars(Model),
                      rows = vars(Year)) +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  ggsave(paste("output/SDM/Future/Maps/", bird_tax[[i]], "_pa_minDispersal.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_sf(
             data = future %>%
               filter(Dispersal == "Min"),
             colour = NA,
             fill = "purple"
           ) +
           theme_bw() +
           facet_grid(cols = vars(Model),
                      rows = vars(Year)) +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  ggsave(paste("output/SDM/Future/Maps/", bird_tax[[i]], "_pa_maxDispersal.pdf", sep = ""),
         oz_map %>%
           ggplot() +
           geom_sf(aes(geometry = geometry), colour = NA) +
           geom_sf(
             data = future %>%
               filter(Dispersal == "Max"),
             colour = NA,
             fill = "purple"
           ) +
           theme_bw() +
           facet_grid(cols = vars(Model),
                      rows = vars(Year)) +
           labs(
             y = "Latitude",
             x = "Longitude",
             title = str_replace(bird_tax[i], "_", " ")
           ),
         width = 7,
         height = 7,
         units = "in")
  
  unlink("occ.bird", recursive = TRUE)
  
  # save in results file
  stats <-
    bind_rows(bm_FindOptimStat(metric.eval = "ROC",
                               sp_training$occ,
                               est),
              bm_FindOptimStat(metric.eval = "BOYCE",
                               sp_training$occ,
                               est),
              bm_FindOptimStat(metric.eval = "TSS",
                               sp_training$occ,
                               est)) %>%
    as_tibble() %>%
    pull(best.stat)
  
  tibble(AUC = stats[[1]],
         BOYCE = stats[[2]],
         TSS = stats[[3]],
         block_size = block_size,
         k = k_value,
         threshold = thr$cutoff/1000,
         n = nrow(points_occ),
         Binomial = bird_tax[[i]]) %>%
    write_csv("output/SDM/bird_summary.csv", append = T)
}