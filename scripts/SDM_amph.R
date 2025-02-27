options(java.parameters = "-Xmx8g")

library(tidyverse)
library(rnaturalearth)
library(sf)
library(raster)
library(biomod2)
library(parallel)
library(pbapply)
library(doSNOW)

cores = detectCores()
cl <- cores[1] - 1 ## set to use all but one thread - replace if necessary

jar <-
  paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '')
require(rJava)

# import Australian amphmal shapefiles
amph_dist <- readRDS("data/amph_dist.RDS")

# create map of Australia
oz_map <-
  ne_countries(country = "Australia",
               scale = 10,
               returnclass = "sf") %>%
  st_transform(st_crs(amph_dist)) %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.01)

# import point locality of Australian amphmal occurrences
amph_points <- readRDS("data/amph_points.RDS") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = st_crs(oz_map))

# generate vector of amphmal taxa
amph_tax <- unique(amph_points$Binomial)
amph_tax <- amph_tax[1:length(amph_tax)]

# import current environmental layers
currentEnv <- stack("data/SDM/Current.tif")
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
  futureEnv[[i]] <- stack(futureLayers[[i]])
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
points_bg_amph <-
  amph_points[sample(1:nrow(amph_points), nrow(amph_points) / 10), ]

points_bg_amph <-
  dismo::gridSample(st_coordinates(points_bg_amph), currentEnv, n = 1)

oz_grid <- st_make_grid(oz_map, cellsize = res(currentEnv) * 10)

amph_sdm <-
  tibble(
    AUC_GAM = character(),
    AUC_GBM = character(),
    AUC_GLM = character(),
    AUC_MAXENT = character(),
    AUC_RF = character(),
    BOYCE_GAM = character(),
    BOYCE_GBM = character(),
    BOYCE_GLM = character(),
    BOYCE_MAXENT = character(),
    BOYCE_RF = character(),
    TSS_GAM = character(),
    TSS_GBM = character(),
    TSS_GLM = character(),
    TSS_MAXENT = character(),
    TSS_RF = character(),
    threshold = character(),
    Binomial = character()
  )
failed <- c()

# iterate loop over all amphmal taxa
for (i in 1:length(amph_tax)) {
  # set seed for reproducibility
  set.seed(42)
  
  # filter points for ith species and generate background points
  points_occ <-
    st_coordinates(st_as_sf(
      amph_points %>% dplyr::filter(Binomial == amph_tax[[i]]),
      crs = st_crs(oz_map)
    )) %>%
    as_tibble()
  
  if (nrow(points_occ) == 1) {
    failed <- c(failed, amph_tax[[i]])
    amph_sdm <- bind_rows(amph_sdm,
                         setNames(c(rep(NA, 16), amph_tax[[i]]), names(amph_sdm)))
    next
  }
  
  # perform spatial thinning on points
  points_occ <-
    dismo::gridSample(as.matrix(points_occ[, 1:2]), currentEnv, n = 1)
  
  pnts_occ <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_occ),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  if (nrow(points_occ) < 10) {
    failed <- c(failed, amph_tax[[i]])
    amph_sdm <- bind_rows(amph_sdm,
                         setNames(c(rep(NA, 16), amph_tax[[i]]), names(amph_sdm)))
    next
  }
  
  # drop from background points all presence points
  pnts_bg <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_bg_amph),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  pnts_logical = lengths(pnts_bg) > 0 & lengths(pnts_occ) == 0
  
  oz_bg <- st_as_sf(oz_grid, "POLYGON")[pnts_logical,]
  
  points_bg <-
    points_bg_amph[lengths(st_intersects(st_as_sf(
      as.data.frame(points_bg_amph),
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
    bind_cols(extract(currentEnv, as(sp_training, "Spatial"))) %>%
    filter_all(~!is.na(.))
  
  training <-
    sp_training %>%
    as_tibble() %>%
    dplyr::select(-geometry)
  
  biomod_data <- BIOMOD_FormatingData(resp.var = training$occ,
                                      expl.var = training[,-1],
                                      resp.xy = sf::st_coordinates(sp_training),
                                      resp.name = "occ.amph",
                                      na.rm = TRUE)
  
  prNum <- as.numeric(table(training$occ)["1"]) # number of presences
  bgNum <- as.numeric(table(training$occ)["0"]) # number of backgrounds
  wt <- ifelse(training$occ == 1, 1, prNum / bgNum)
  
  # generate SDMs
  biomod_options <- BIOMOD_ModelingOptions(GBM = list(distribution = "bernoulli",
                                                      n.trees = 2500,
                                                      interaction.depth = 3,
                                                      n.minobsinnode = 5,
                                                      learning.rate = 0.01,
                                                      bag.fraction = 0.75,
                                                      train.fraction = 1,
                                                      keep.data = FALSE,
                                                      verbose = FALSE,
                                                      n.cores = 1),
                                           GLM = list(type = 'quadratic',
                                                      interaction.level = 0,
                                                      myFormula = NULL,
                                                      test = 'AIC',
                                                      family = binomial(link = 'logit'),
                                                      mustart = 0.5,
                                                      control = glm.control(epsilon = 1e-08, maxit = 50, trace = FALSE)),
                                           MAXENT = list(path_to_maxent.jar = jar, 
                                                         memory_allocated = 512,
                                                         initial_heap_size = NULL,
                                                         maximum_heap_size = NULL,
                                                         background_data_dir = 'default',
                                                         maximumbackground = 'default',
                                                         maximumiterations = 200,
                                                         visible = FALSE,
                                                         linear = TRUE,
                                                         quadratic = TRUE,
                                                         product = TRUE,
                                                         threshold = TRUE,
                                                         hinge = TRUE,
                                                         lq2lqptthreshold = 80,
                                                         l2lqthreshold = 10,
                                                         hingethreshold = 15,
                                                         beta_threshold = -1,
                                                         beta_categorical = -1,
                                                         beta_lqp = -1,
                                                         beta_hinge = -1,
                                                         betamultiplier = 1,
                                                         defaultprevalence = 0.5))
  
  biomod_model_out <- BIOMOD_Modeling(biomod_data,
                                      models = c('GBM','GLM','GAM','MAXENT','RF'),
                                      bm.options = biomod_options,
                                      CV.strategy = 'block',
                                      metric.eval = c('ROC', 'TSS'),
                                      weights = wt,
                                      seed.val = 42,
                                      nb.cpu = cl)
  
  # generate ensemble model
  ens_model <- BIOMOD_EnsembleModeling(biomod_model_out,
                                       em.by = "all",
                                       em.algo = 'EMwmean',
                                       metric.eval = c('ROC', 'TSS'))
  
  # plot partial response curves
  pdf(paste("output/SDM/Response Curves/", amph_tax[[i]], ".pdf", sep = ""),
      useDingbats = F)
  bm_PlotResponseCurves(ens_model,
                        models.chosen = get_built_models(ens_model)[[2]])
  dev.off()
  
  # projections
  ens_pred <- BIOMOD_EnsembleForecasting(bm.em = ens_model,
                                         models.chosen = get_built_models(ens_model)[[2]],
                                         proj.name = "ens",
                                         new.env = currentEnv,
                                         build.clamping.mask = FALSE,
                                         output.format = ".tif",
                                         nb.cpu = cl)
  
  pr <- raster::raster(get_predictions(ens_pred))
  names(pr) = "layer"
  
  pr_df = as.data.frame(pr/1000, xy = T)
  
  # plot habitat suitability probability
  pdf(paste("output/SDM/Current/Maps/", amph_tax[[i]], "_prob.pdf", sep = ""),
      useDingbats = F)
  print(
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
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # extract model estimated suitability
  est = raster::extract(pr, st_coordinates(sp_training))
  # detect threshold for presence/absence
  thr = bm_FindOptimStat(metric.eval = "TSS", obs = sp_training$occ, fit = est)
  # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
  pr_thr = pr > thr$cutoff
  
  # save in results file
  amph_sdm <-
    amph_sdm %>%
    bind_rows(
      setNames(c(BIOMOD_PresenceOnly(biomod_model_out) %>%
                   as_tibble() %>%
                   mutate(metric.eval = case_when(metric.eval == "ROC" ~ "AUC",
                                                  T ~ metric.eval)) %>%
                   group_by(metric.eval, algo) %>%
                   summarise(validation = mean(validation)) %>%
                   filter(metric.eval %in% c("AUC", "TSS", "BOYCE")) %>%
                   pull(validation),
                 thr$cutoff,
                 amph_tax[[i]]),
               names(amph_sdm))
    )
  
  # convert current presence-absence map to polygon
  current <-
    rasterToPolygons(
      pr_thr,
      fun = function(x) {
        x == 1
      },
      dissolve = T
    ) %>%
    st_as_sf() %>%
    st_cast("POLYGON")
  
  # find which areas of suitable habitat do not overlap occurrences
  pnts_sf <-
    st_intersects(current,
                  st_as_sf(
                    amph_points %>% dplyr::filter(Binomial == amph_tax[[i]]),
                    coords = c("X", "Y"),
                    crs = st_crs(oz_map)
                  ))
  pnts_logical = lengths(pnts_sf) > 0
  
  # filter out areas of suitable habitat which do not overlap landmasses with occurrences
  current <- current[pnts_logical, ] %>%
    st_union() %>%
    st_as_sf() %>%
    dplyr::mutate(Binomial = amph_tax[[i]])
  
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
    dplyr::mutate(Binomial = amph_tax[[i]])
  
  # identify landmasses
  landmasses <- st_intersects(oz_map %>%
                                st_cast("POLYGON"),
                              st_as_sf(
                                amph_points %>% dplyr::filter(Binomial == amph_tax[[i]]),
                                coords = c("X", "Y"),
                                crs = st_crs(oz_map)
                              ))
  landmasses <- st_cast(oz_map, "POLYGON")[lengths(landmasses) > 0, ] %>%
    st_union()
  
  # export habitat suitability under current climate conditions
  current_predictions <- list(pr/1000,
                              pr_thr)
  names(current_predictions) <-
    c("predicted_probability", "presence_absence")
  writeRaster(stack(current_predictions),
              paste("output/SDM/Current/Rasters/", amph_tax[[i]], ".tif", sep = ""))
  write_sf(current,
           paste("output/SDM/Current/Shapefiles/", amph_tax[[i]], ".shp", sep = ""))
  
  # plot presence/absence based on threshold under current climate conditions
  pdf(paste("output/SDM/Current/Maps/", amph_tax[[i]], "_pa.pdf", sep = ""),
      useDingbats = F)
  print(
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
          dplyr::filter(amph_points, Binomial == amph_tax[[i]]),
          crs = st_crs(oz_map)
        ),
        size = .1
      ) +
      theme_bw() +
      labs(
        y = "Latitude",
        x = "Longitude",
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # make future predictions
  # over all 8 climate projections
  future_predictions <- list()
  for(k in 1:length(futureEnv)) {
    # make predictions under future conditions
    ens_predf <- BIOMOD_EnsembleForecasting(bm.em = ens_model,
                                            models.chosen = get_built_models(ens_model)[[2]],
                                            proj.name = "ens",
                                            new.env = futureEnv[[k]],
                                            build.clamping.mask = FALSE,
                                            output.format = ".tif",
                                            nb.cpu = cl)
    
    prf <- raster::raster(get_predictions(ens_predf))
    names(prf) = "layer"
    
    # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
    prf_thr = prf > thr$cutoff
    
    # save habitat suitability under future climate conditions
    future_preds <- list(prf/1000,
                         prf_thr)
    names(future_preds) <-
      c("predicted_probability", "presence_absence")
    
    future_predictions[[k]] <- future_preds
  }
  
  names(future_predictions) <-
    str_remove_all(futureLayers, "data/SDM/|.tif")
  
  # extract all future predicted presence-absence maps
  future <- lapply(future_predictions, "[[", "presence_absence")
  
  cl <- makeCluster(3)
  registerDoSNOW(cl)
  
  future_all <- foreach(scenario = 1:3,
                        .packages = c("tidyverse", "sf", "raster"),
                        .export = c("future",
                                    "oz_map", "extant_buff",
                                    "amph_tax", "models", "years"),
                        .inorder = TRUE,
                        .errorhandling = "pass",
                        .verbose = FALSE) %dopar% {
                          
                          # iterate over all future projections
                          kmax = 4
                          for (m in 1:2) {
                            if(scenario == 1){
                              ## standard dispersal model ##
                              # for each GCM and SSP set the current buffered suitable habitat to be the extant
                              current_buff <- extant_buff
                              for (k in (kmax - 3):kmax) {
                                # for each time-step
                                # if no suitable habitat exists in time step, save empty polygon
                                if (maxValue(future[[k]]) == 0)
                                  future_next <- st_sf(x = st_sfc(st_multipolygon()),
                                                       crs = st_crs(oz_map)) %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                else {
                                  # import presence-absence map for the next time-step and convert to polygon
                                  # capture instance where only one pixel remains
                                  if (length(which(values(future[[k]]) == 1)) == 1) {
                                    future_next <- rasterToPolygons(
                                      future[[k]],
                                      dissolve = T
                                    ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      filter(layer > 0) %>%
                                      st_cast("POLYGON")
                                  } else {
                                    future_next <-
                                      rasterToPolygons(
                                        future[[k]],
                                        fun = function(x) {
                                          x == 1
                                        },
                                        dissolve = T
                                      ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      st_cast("POLYGON")
                                  }
                                  
                                  # find and filter out areas of suitable habitat which do not intersect with buffered suitable habitat for current time-step
                                  future_oz <- st_intersects(future_next, current_buff)
                                  future_logical = lengths(future_oz) > 0
                                  
                                  future_next <- future_next[future_logical, ] %>%
                                    st_union() %>%
                                    st_as_sf() %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                }
                                
                                # save filtered suitable habitat
                                future[[k]] <- future_next
                                
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
                                  dplyr::mutate(Binomial = amph_tax[[i]],
                                                Model = models[[k]],
                                                Year = years[[k]])
                              }
                            }
                            
                            if(scenario == 2){
                              ## no dispersal ##
                              for (k in (kmax - 3):kmax) {
                                # for each time-step
                                # if no suitable habitat exists in time step, save empty polygon
                                if (maxValue(future[[k]]) == 0)
                                  future_next <- st_sf(x = st_sfc(st_multipolygon()),
                                                       crs = st_crs(oz_map)) %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                else {
                                  # import presence-absence map for the next time-step and convert to polygon
                                  # capture instance where only one pixel remains
                                  if (length(which(values(future[[k]]) == 1)) == 1) {
                                    future_next <- rasterToPolygons(
                                      future[[k]],
                                      dissolve = T
                                    ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      filter(layer > 0) %>%
                                      st_cast("POLYGON")
                                  } else {
                                    future_next <-
                                      rasterToPolygons(
                                        future[[k]],
                                        fun = function(x) {
                                          x == 1
                                        },
                                        dissolve = T
                                      ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      st_cast("POLYGON")
                                  }
                                  
                                  # find and filter out areas of suitable habitat which do not intersect with buffered suitable habitat for current time-step
                                  future_oz <- st_intersects(future_next, current)
                                  future_logical = lengths(future_oz) > 0
                                  
                                  future_next <- future_next[future_logical, ] %>%
                                    st_union() %>%
                                    st_as_sf() %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                }
                                
                                # save filtered suitable habitat
                                future[[k]] <- current <- future_next
                              }
                            }
                            
                            if(scenario == 3){
                              ## max dispersal ##
                              for (k in (kmax - 3):kmax) {
                                # for each time-step
                                # if no suitable habitat exists in time step, save empty polygon
                                if (maxValue(future[[k]]) == 0)
                                  future_next <- st_sf(x = st_sfc(st_multipolygon()),
                                                       crs = st_crs(oz_map)) %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                else {
                                  # import presence-absence map for the next time-step and convert to polygon
                                  # capture instance where only one pixel remains
                                  if (length(which(values(future[[k]]) == 1)) == 1) {
                                    future_next <- rasterToPolygons(
                                      future[[k]],
                                      dissolve = T
                                    ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      filter(layer > 0) %>%
                                      st_cast("POLYGON")
                                  } else {
                                    future_next <-
                                      rasterToPolygons(
                                        future[[k]],
                                        fun = function(x) {
                                          x == 1
                                        },
                                        dissolve = T
                                      ) %>%
                                      st_as_sf(crs = st_crs(oz_map)) %>%
                                      st_cast("POLYGON")
                                  }
                                  
                                  # filter out areas of suitable habitat which do not overlap landmasses with occurrences
                                  future_oz <- st_intersects(future_next, landmasses)
                                  future_logical = lengths(future_oz) > 0
                                  
                                  future_next <- future_next[future_logical, ] %>%
                                    st_union() %>%
                                    st_as_sf() %>%
                                    dplyr::mutate(Binomial = amph_tax[[i]],
                                                  Model = models[[k]],
                                                  Year = years[[k]])
                                }
                                
                                # save filtered suitable habitat
                                future[[k]] <- future_next
                              }
                            }
                            
                            kmax = kmax + 4
                          }
                          future
                        }
  stopCluster(cl)
  
  future <- bind_rows(future_all[[1]])
  future_maxdisp <- bind_rows(future_all[[3]])
  future_nodisp <- bind_rows(future_all[[2]])
  
  future <- bind_rows(future %>% mutate(Dispersal = "Mean"),
                      future_nodisp %>% mutate(Dispersal = "Min"),
                      future_maxdisp %>% mutate(Dispersal = "Max"))
  
  # export habitat suitability under future climate conditions
  writeRaster(stack(lapply(1:length(future_predictions),
                           function(x)
                             setNames(stack(future_predictions[[x]]),
                                      paste(names(future_predictions)[x], names(future_predictions[[x]]), sep = "_")))),
              paste("output/SDM/Future/Rasters/", amph_tax[[i]], ".tif", sep = ""))
  write_sf(future,
           paste("output/SDM/Future/Shapefiles/", amph_tax[[i]], ".shp", sep = ""))
  
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
  pdf(
    paste("output/SDM/Future/Maps/", amph_tax[[i]], "_prob.pdf", sep = ""),
    useDingbats = F,
    width = 20,
    height = 14
  )
  print(
    oz_map %>%
      ggplot() +
      geom_sf(aes(geometry = geometry), colour = NA) +
      geom_raster(data = prf_df %>% drop_na(), aes(
        x = x, y = y, fill = layer
      )) +
      theme_bw() +
      scale_fill_gradientn(colours = viridis::viridis(99),
                           name = "Predicted Habitat Suitability") +
      facet_grid(cols = vars(model),
                 rows = vars(year)) +
      labs(
        y = "Latitude",
        x = "Longitude",
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # plot all future presence/absence maps
  pdf(
    paste("output/SDM/Future/Maps/", amph_tax[[i]], "_pa_meanDispersal.pdf", sep = ""),
    useDingbats = F,
    width = 20,
    height = 14
  )
  print(
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
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  pdf(
    paste("output/SDM/Future/Maps/", amph_tax[[i]], "_pa_minDispersal.pdf", sep = ""),
    useDingbats = F,
    width = 20,
    height = 14
  )
  print(
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
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  pdf(
    paste("output/SDM/Future/Maps/", amph_tax[[i]], "_pa_maxDispersal.pdf", sep = ""),
    useDingbats = F,
    width = 20,
    height = 14
  )
  print(
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
        title = str_replace(amph_tax[i], "_", " ")
      )
  )
  dev.off()
  
  unlink("occ.amph", recursive = TRUE)
}

# export results file
write_csv(drop_na(amph_sdm), "output/SDM/amph_summary.csv")
