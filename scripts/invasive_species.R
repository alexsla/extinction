options(java.parameters = "-Xmx8g")

library(tidyverse)
library(sf)
library(rnaturalearth)
library(raster)
library(galah)
library(ape)
library(taxize)
library(slga)
library(foreach)
library(doParallel)
library(ENMeval)

sf_use_s2(FALSE)

jar <-
  paste(system.file(package = "dismo"), "/java/maxent.jar", sep = '')
require(rJava)

Sys.setenv(JAVA_HOME = 'C:\\Program Files\\Java\\jre1.8.0_341')

# create map of Australia
oz_map <-
  ne_countries(country = "Australia",
               scale = 10,
               returnclass = "sf") %>%
  st_transform("wgs84") %>%
  st_simplify(preserveTopology = FALSE, dTolerance = 0.01)

## create dataset of occurrence points in Australia
# generate vector of taxa
inv_tax <-
  c(
    "Rhinella marina",
    "Felis catus",
    "Vulpes vulpes",
    "Camelus dromedarius",
    "Oryctolagus cuniculus"
  )

# download occurrences per species
cores = detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)

foreach(i = 1:length(inv_tax),
        .packages = c("tidyverse", "galah", "sf")) %dopar% {
          # test if already downloaded (if so, skip to next iteration)
          if (!inv_tax[[i]] %in% str_remove(dir("data/invasive/citations/"), ".txt")) {
            # set email registered to ALA
            galah_config(email = "your_email")
            
            result <- tryCatch(
              galah_call() |>
                galah_identify(str_replace(inv_tax[[i]], "_", " ")) |>
                galah_filter(profile = "ALA", year > 1990, coordinateUncertaintyInMeters < 4500) |> # filter after 1990, uncertainty of less than 4.5 km (resolution of climate data)
                galah_select(basisOfRecord, group = "basic") |>
                atlas_occurrences(mint_doi = T),
              error = function(e)
                e
            )
            
            if (inherits(result, "error")) {
              result <- tryCatch(
                galah_call() |>
                  galah_identify(
                    dplyr::filter(
                      bind_rows(rep_synonym, bird_synonym),
                      Binomial == str_replace(inv_tax[[i]], "_", " ")
                    )$synonym
                  ) |>
                  galah_filter(
                    profile = "ALA",
                    year > 1990,
                    coordinateUncertaintyInMeters < 4500
                  ) |>
                  galah_select(basisOfRecord, group = "basic") |>
                  atlas_occurrences(mint_doi = T),
                error = function(e)
                  e
              )
            }
            
            if (!inherits(result, "error")) {
              # if there are occurrence points left after filtering
              # convert to spatial object
              result <-
                st_as_sf(
                  result %>% drop_na(decimalLatitude, decimalLongitude, scientificName),
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = st_crs(oz_map)
                )
              
              # filter out points which do not fall within Australia
              result_points <- st_join(result, oz_map, join = st_within) %>%
                dplyr::filter(!is.na(sovereignt))
              
              if (dim(result_points)[1] > 0) {
                # if there are occurrence points left after filtering
                # convert to tibble
                result_points <-
                  bind_cols(as_tibble(result_points),
                            as_tibble(st_coordinates(result_points))) %>%
                  rename(decimalLongitude = X,
                         decimalLatitude = Y) %>%
                  select(c(names(result), decimalLatitude, decimalLongitude),
                         -geometry)
                
                # export points
                write_csv(
                  result_points,
                  paste("data/invasive/occurrences/", inv_tax[[i]], ".csv", sep = "")
                )
                
                # export citation for occurrences
                write_file(
                  atlas_citation(result),
                  paste("data/invasive/citations/", inv_tax[[i]], ".txt", sep = "")
                )
              }
            }
          }
        }
stopCluster(cl)

# convert to single dataframe
inv_points <-
  bind_rows(lapply(dir("data/invasive/occurrences/"), function(x)
    read_csv(paste("data/invasive/occurrences/", x, sep = "")) %>% mutate(Binomial = str_remove(x, ".csv"))))

# export points
write_rds(inv_points, "data/inv_points.RDS")

# export citations
inv_citations <-
  tibble(Binomial = str_remove(dir("data/invasive/citations/"), ".txt"),
         citation = sapply(dir("data/invasive/citations/", full.names = T), readLines))
write_csv(inv_citations, "data/inv_citations.csv")

# SDMs
# import point locality of Australian invmal occurrences
inv_points <- readRDS("data/inv_points.RDS") %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = st_crs(oz_map))

# drop fox localities from Tasmania (erroneous)
inv_points <-
  inv_points[which(!(
    st_coordinates(inv_points)[, 2] < -40 &
      inv_points$Binomial == "Vulpes vulpes"
  )), ]

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
    "CNRM-ESM2-1 SSP1.26",
    "CNRM-ESM2-1 SSP1.26",
    "CNRM-ESM2-1 SSP1.26",
    "CNRM-ESM2-1 SSP1.26",
    "CNRM-ESM2-1 SSP2.45",
    "CNRM-ESM2-1 SSP2.45",
    "CNRM-ESM2-1 SSP2.45",
    "CNRM-ESM2-1 SSP2.45",
    "CNRM-ESM2-1 SSP5.85",
    "CNRM-ESM2-1 SSP5.85",
    "CNRM-ESM2-1 SSP5.85",
    "CNRM-ESM2-1 SSP5.85",
    "MICRO6 SSP1.26",
    "MICRO6 SSP1.26",
    "MICRO6 SSP1.26",
    "MICRO6 SSP1.26",
    "MICRO6 SSP2.45",
    "MICRO6 SSP2.45",
    "MICRO6 SSP2.45",
    "MICRO6 SSP2.45",
    "MICRO6 SSP5.85",
    "MICRO6 SSP5.85",
    "MICRO6 SSP5.85",
    "MICRO6 SSP5.85"
  )
years <- c(
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100",
  "2021-2040",
  "2041-2060",
  "2061-2080",
  "2081-2100"
)

# generate vectors with parameters for SDM building and evaulation
rm_v = c(0.1, 0.5, 1, 2, 4)
fc_v = c('L', 'LQ', 'LQP', 'H', 'Q', 'QH')
cv_method = "block"

# generate background points
# set seed for reproducibility
set.seed(42)

points_bg_inv <-
  dismo::gridSample(st_coordinates(inv_points), currentEnv, n = 1)

oz_grid <- st_make_grid(oz_map, cellsize = res(currentEnv) * 10)

inv_sdm <-
  tibble(
    fc = character(),
    rm = numeric(),
    AUC = numeric(),
    threshold = numeric(),
    Binomial = character()
  )
failed <- c()
# iterate loop over all invmalian taxa
for (i in 1:length(inv_tax)) {
  # set seed for reproducibility
  set.seed(42)
  
  # filter points for ith species and generate background points
  points_occ <-
    st_coordinates(st_as_sf(
      inv_points %>% dplyr::filter(Binomial == inv_tax[[i]]),
      crs = st_crs(oz_map)
    )) %>%
    as_tibble()
  
  # perform spatial thinning on points
  points_occ <-
    dismo::gridSample(as.matrix(points_occ[, 1:2]), currentEnv, n = 1)
  
  pnts_occ <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_occ),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  # drop from background points all presence points
  pnts_bg <-
    st_intersects(oz_grid, st_as_sf(
      as.data.frame(points_bg_inv),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ))
  
  pnts_logical = lengths(pnts_bg) > 0 & lengths(pnts_occ) == 0
  
  oz_bg <- st_as_sf(oz_grid, "POLYGON")[pnts_logical,]
  
  points_bg <-
    points_bg_inv[lengths(st_intersects(st_as_sf(
      as.data.frame(points_bg_inv),
      coords = c("X", "Y"),
      crs = st_crs(oz_map)
    ), oz_bg)) > 0, ]
  
  # run and perform cross-evaulatuion on maxent models with different parameters
  eval <- tryCatch(
    ENMevaluate(
      occs = points_occ,
      envs = currentEnv,
      bg = points_bg,
      partitions = "block",
      algorithm = 'maxent.jar',
      tune.args = list(fc = fc_v, rm = rm_v),
      parallel = TRUE,
      numCores = 8
    ),
    error = function(e)
      e
  )
  
  if (inherits(eval, "error")) {
    failed <- c(failed, inv_tax[[i]])
    inv_sdm[i, 5] <- inv_tax[[i]]
    next
  }
  
  # select best model based on AICc scores
  bestmod = which(eval@results$AICc == min(eval@results$AICc, na.rm = T))
  
  if (length(bestmod) < 1) {
    failed <- c(failed, inv_tax[[i]])
    inv_sdm[i, 5] <- inv_tax[[i]]
    next
  }
  
  if (length(bestmod) > 1) {
    # if more than one model is tied for best, choose most parsimonious (fewest feature class parameters)
    equal.models <- eval@results[bestmod, 1]
    ln.equal.models <- sapply(equal.models, str_length)
    bestmod <-
      bestmod[which(ln.equal.models == min(ln.equal.models))]
    if (length(bestmod) > 1) {
      # if more than one model is tied for most parsimonious, select one randomly
      set.seed(42)
      bestmod <- bestmod[sample(length(bestmod), 1)]
    }
  }
  
  # plot partial response curves
  pdf(paste("output/SDM/Response Curves/", inv_tax[[i]], ".pdf", sep = ""),
      useDingbats = F)
  dismo::response(eval@models[[bestmod]])
  dev.off()
  
  # generate habitat suitability predictions under current climate conditions
  pr = predict(currentEnv, eval@models[[bestmod]], type = 'cloglog')
  pr_df = as.data.frame(pr, xy = T)
  
  # plot habitat suitability probability
  pdf(paste("output/SDM/Current/Maps/", inv_tax[[i]], "_prob.pdf", sep = ""),
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
        title = str_replace(inv_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # extract model estimated suitability for occurrence localities
  est.loc = raster::extract(pr, eval@occs[, c("X", "Y")])
  # extract model estimated suitability for background localities
  est.bg = raster::extract(pr, eval@bg[, c("X", "Y")])
  # evaluate predictive ability of model
  ev = dismo::evaluate(est.loc, est.bg)
  # detect possible thresholds for presence/absence
  thr = dismo::threshold(ev)
  # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
  pr_thr = pr > thr$spec_sens
  
  # save in results file
  inv_sdm[i, 1] <- as.character(eval@results[bestmod, ]$fc)
  inv_sdm[i, 2] <- as.numeric(eval@results[bestmod, ]$rm)
  inv_sdm[i, 3] <- eval@results[bestmod, ]$auc.train
  inv_sdm[i, 4] <- thr$spec_sens
  inv_sdm[i, 5] <- inv_tax[[i]]
  
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
                    inv_points %>% dplyr::filter(Binomial == inv_tax[[i]]),
                    coords = c("X", "Y"),
                    crs = st_crs(oz_map)
                  ))
  pnts_logical = lengths(pnts_sf) > 0
  
  # filter out areas of suitable habitat which do not overlap landmasses with occurrences
  current <- current[pnts_logical, ] %>%
    st_union() %>%
    st_as_sf() %>%
    dplyr::mutate(Binomial = inv_tax[[i]])
  
  # create a buffer around current suitable habitat to allow for dispersal
  current_buff <- current %>%
    st_transform(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs") %>%
    st_buffer(dist = 150) %>% # numeric buffer of 100 km to capture points just outside IUCN ranges
    st_transform(crs = st_crs(oz_map)) %>%
    st_make_valid() %>% # fix geometry errors due to buffer
    st_intersection(oz_map) %>%
    st_cast("POLYGON")
  
  # find and filter out areas of buffer that do not occurr on the same landmasses
  current_oz <- st_intersects(current_buff, current)
  current_logical = lengths(current_oz) > 0
  
  extant_buff <- current_buff[current_logical, ] %>%
    st_union() %>%
    st_as_sf() %>%
    dplyr::mutate(Binomial = inv_tax[[i]])
  
  # export habitat suitability under current climate conditions
  current_predictions <- list(pr,
                              pr_thr)
  names(current_predictions) <-
    c("predicted_probability", "presence_absence")
  write_rds(current_predictions,
            paste("output/SDM/Current/Rasters/", inv_tax[[i]], sep = ""))
  write_sf(current,
           paste("output/SDM/Current/Shapefiles/", inv_tax[[i]], ".shp", sep = ""))
  
  # plot presence/absence based on threshold under current climate conditions
  pdf(paste("output/SDM/Current/Maps/", inv_tax[[i]], "_pa.pdf", sep = ""),
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
          dplyr::filter(inv_points, Binomial == inv_tax[[i]]),
          crs = st_crs(oz_map)
        ),
        size = .1
      ) +
      theme_bw() +
      labs(
        y = "Latitude",
        x = "Longitude",
        title = str_replace(inv_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # make future predictions
  cl <- makeCluster(12)
  registerDoParallel(cl)
  
  # parallelize over all 24 climate projections
  future_predictions <-
    foreach(
      k = 1:length(futureEnv),
      .packages = c("ENMeval", "maxnet", "raster")
    ) %dopar% {
      # make predictions under future conditions
      prf = predict(futureEnv[[k]], eval@models[[bestmod]], type = 'cloglog')
      # convert predicted probability to presence-absence based on threshold maximising the sum of sensitivity and specificity
      prf_thr = prf > thr$spec_sens
      
      # save habitat suitability under future climate conditions
      future_preds <- list(prf,
                           prf_thr)
      names(future_preds) <-
        c("predicted_probability", "presence_absence")
      
      future_preds
    }
  stopCluster(cl)
  
  names(future_predictions) <-
    str_remove_all(futureLayers, "data/SDM/|.tif")
  
  # extract all future predicted presence-absence maps
  future <- lapply(future_predictions, "[[", "presence_absence")
  
  # iterate over all future projections
  kmax = 4
  for (m in 1:6) {
    # for each GCM and SSP set the current buffered suitable habitat to be the extant
    current_buff <- extant_buff
    for (k in (kmax - 3):kmax) {
      # for each time-step
      # if no suitable habitat exists in time step, save empty polygon
      if (maxValue(future[[k]]) == 0)
        future_next <- st_sf(x = st_sfc(st_multipolygon()),
                             crs = st_crs(oz_map)) %>%
          dplyr::mutate(Binomial = inv_tax[[i]],
                        Model = models[[k]],
                        Year = years[[k]])
      else {
        # import presence-absence map for the next time-step and convert to polygon
        future_next <-
          rasterToPolygons(
            future[[k]],
            fun = function(x) {
              x == 1
            },
            dissolve = T
          ) %>%
          st_as_sf() %>%
          st_cast("POLYGON")
        
        # find and filter out areas of suitable habitat which do not intersect with buffered suitable habitat for current time-step
        future_oz <- st_intersects(future_next, current_buff)
        future_logical = lengths(future_oz) > 0
        
        future_next <- future_next[future_logical, ] %>%
          st_union() %>%
          st_as_sf() %>%
          dplyr::mutate(Binomial = inv_tax[[i]],
                        Model = models[[k]],
                        Year = years[[k]])
      }
      
      # save filtered suitable habitat
      future[[k]] <- future_next
      
      # update buffered current suitable habitat for next iteration
      current_buff <- future_next %>%
        st_transform(crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=km +no_defs") %>%
        st_buffer(dist = 150) %>% # numeric buffer of 100 km to capture points just outside IUCN ranges
        st_transform(crs = st_crs(oz_map)) %>%
        st_make_valid() %>% # fix geometry errors due to buffer
        st_intersection(oz_map) %>%
        st_cast("POLYGON")
      
      current_oz <- st_intersects(current_buff, future_next)
      current_logical = lengths(current_oz) > 0
      
      current_buff <- current_buff[current_logical, ] %>%
        st_union() %>%
        st_as_sf() %>%
        dplyr::mutate(Binomial = inv_tax[[i]],
                      Model = models[[k]],
                      Year = years[[k]])
    }
    kmax = kmax + 4
  }
  future <- bind_rows(future)
  
  # export habitat suitability under future climate conditions
  write_rds(future_predictions,
            paste("output/SDM/Future/Rasters/", inv_tax[[i]], sep = ""))
  write_sf(future,
           paste("output/SDM/Future/Shapefiles/", inv_tax[[i]], ".shp", sep = ""))
  
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
    paste("output/SDM/Future/Maps/", inv_tax[[i]], "_prob.pdf", sep = ""),
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
        title = str_replace(inv_tax[i], "_", " ")
      )
  )
  dev.off()
  
  # plot all future presence/absence maps
  pdf(
    paste("output/SDM/Future/Maps/", inv_tax[[i]], "_pa.pdf", sep = ""),
    useDingbats = F,
    width = 20,
    height = 14
  )
  print(
    oz_map %>%
      ggplot() +
      geom_sf(aes(geometry = geometry), colour = NA) +
      geom_sf(
        data = future,
        colour = NA,
        fill = "purple"
      ) +
      theme_bw() +
      facet_grid(cols = vars(Model),
                 rows = vars(Year)) +
      labs(
        y = "Latitude",
        x = "Longitude",
        title = str_replace(inv_tax[i], "_", " ")
      )
  )
  dev.off()
}

# export results file
write_csv(drop_na(inv_sdm), "output/SDM/inv_summary.csv")
