#### MAMMALS ####
mam_sdm <- tet_sdm[which(tet_sdm %in% filter(tet_iucn, Class == "Mammal")$Binomial)]

mam_current <- lapply(mam_sdm, function(x) read_sf(paste("output/SDM/Current/Shapefiles/", x, ".shp", sep = ""))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))
mam_future <- lapply(mam_sdm, function(x) read_sf(paste("output/SDM/Future/Shapefiles/", x, ".shp", sep = "")) %>%
                       mutate(SSP = sapply(str_split(Model, " "), "[[", 2),
                              Model = sapply(str_split(Model, " "), "[[", 1),
                              Year = sapply(str_split(Year, "-"), "[[", 2))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))


mam_current_biomes <- read_csv("output/tet_biomes.csv") %>%
  filter(Binomial %in% mam_current$Binomial)

mam_future_biomes <- bind_rows(lapply(1:24, function(x) read_csv(paste("output/tet_biomes_", x, ".csv", sep = "")))) %>%
  filter(Binomial %in% mam_current$Binomial)

# import mammal data
EPBC <- read_csv("data/EPBC.csv") %>% mutate(Binomial = str_replace(Binomial, " ", "_"))

mam_data <- read_csv("output/tet_lu_hp.csv") %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Mass)) %>%
  left_join(tet_dist %>% as_tibble() %>% select(-geometry)) %>%
  left_join(EPBC %>% select(-Class)) %>%
  left_join(bind_rows(mam_current_biomes, mam_future_biomes) %>% select(-area)) %>%
  mutate(category = case_when(is.na(EPBC) ~ category,
                              TRUE ~ EPBC),
         category = case_when(Year == 2010 ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-EPBC) %>%
  filter(Binomial %in% mam_current$Binomial)

mam_data[which(is.na(mam_data[,"category"])), "category"] <- "NE"
mam_data[which(mam_data$Binomial == "Chelodina_novaeguineae"), "category"] <- "LC"

# drop species with extremely wrong predicted ranges
mam_data <- mam_data %>%
  left_join(tet_dist  %>%
              st_transform("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=61.55,-10.87,-40.19,-39.4924,-32.7221,-32.8979,-9.99400000001316 +units=m +no_defs +type=crs") %>%
              mutate(area_iucn = as.numeric(units::set_units(st_area(geometry), "km2"))) %>%
              as_tibble() %>%
              select(Binomial, area_iucn)) %>%
  mutate(pred_real = log(area/area_iucn),
         category = case_when(abs(pred_real) >= 3 & Year == 2010 ~ "NE",
                              TRUE ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-c(pred_real, area_iucn)) %>%
  drop_na(Mass)

# combine all data and omit species without ranges
mam_all <- bind_rows(mam_current, mam_future) %>%
  filter(Binomial %in% mam_data$Binomial) %>%
  mutate(id = str_replace(paste(Binomial, Model, SSP, Year, sep = "_"), "NA_NA_NA", "2020"),
         Model = case_when(is.na(Model) ~ "current",
                           TRUE ~ Model),
         Year = case_when(is.na(Year) ~ "current",
                          TRUE ~ Year),
         SSP = case_when(is.na(SSP) ~ "current",
                         TRUE ~ SSP))

mam_all_sp <- st_centroid(mam_all)

empty <- which(is.na(st_coordinates(mam_all_sp)))
mam_all_sp <- mam_all_sp[-empty, ]
mam_all <- mam_all[-empty, ]

# format for xgboost
pred <-
  cbind(mam_data %>% filter(Year == 2010) %>% select(-c(Binomial, SSP, Year, Model)),
        st_coordinates(mam_all_sp %>% filter(Year == "current")))

rownames(pred) <- filter(mam_all_sp, Year == "current")$Binomial

#### RUN AUTOMATED ASSESSMENT
biomes <- list(thrt = list(tropical = c(1, 7),
                           temperate = c(8, 4, 10),
                           dry = c(12, 13)),
               cr = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               envu = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               nthrt = list(tropical = c(1, 7),
                            temperate = c(8, 4, 10),
                            dry = c(12, 13)))

# seperate into assessed and non-assessed species
pred_ne_dd <- 
  pred[pred$category %in% c("NE", "DD"),]
pred <- 
  pred[!pred$category %in% c("NE", "DD", "EX", "EW"),]

aa_thrt <- auto_ass(pred = pred[, -c(12)],
                    biomes = biomes$thrt,
                    path = "output/AA/mammals/",
                    cat = "threatened")

# predict specific categories
pred_thrt <- pred[pred$threatened == "yes", ]
pred_nthrt <- pred[pred$threatened == "no", ]

### CR vs VU+EN
cr_ind <- pred_thrt[pred_thrt$threatened == "yes", ]$category == "CR"

pred_thrt$cr <- 0
pred_thrt$cr[cr_ind] <- 1

aa_cr <- auto_ass(pred = pred_thrt[, -c(12, 14)],
                  biomes = biomes$cr,
                  path = "output/AA/mammals/",
                  cat = "cr")


pred_thrt_envu <- pred_thrt[pred_thrt$cr == 0, ]

### VU vs EN
en_ind <- pred_thrt_envu[pred_thrt_envu$threatened == "yes", ]$category == "EN"

pred_thrt_envu$envu <- 0
pred_thrt_envu$envu[en_ind] <- 1

pred_thrt_envu <- 
  pred_thrt_envu[!names(pred_thrt_envu) %in% c("cr")]

aa_envu <- auto_ass(pred = pred_thrt_envu[, -c(12, 14)],
                    biomes = biomes$envu,
                    path = "output/AA/mammals/",
                    cat = "envu")

### NT vs LC
nt_ind <- pred_nthrt[pred_nthrt$threatened == "no", ]$category == "NT"

pred_nthrt$nt <- 0
pred_nthrt$nt[nt_ind] <- 1

aa_nthrt <- auto_ass(pred = pred_nthrt[, -c(12, 14)],
                     biomes = biomes$nthrt,
                     path = "output/AA/mammals/",
                     cat = "nt")


## predict NE & DD species
ne_dd <- pred_ass(pred = rbind(pred, pred_ne_dd)[, -12],
                   biomes = biomes,
                   models = list(aa_thrt,
                                 aa_cr,
                                 aa_envu,
                                 aa_nthrt))

mam_iucn_current <- mam_data %>% filter(Year == 2010) %>%
  mutate(category = case_when(Binomial %in% ne_dd$lc ~ "LC",
                              Binomial %in% ne_dd$nt ~ "NT",
                              Binomial %in% ne_dd$vu ~ "VU",
                              Binomial %in% ne_dd$en ~ "EN",
                              Binomial %in% ne_dd$cr ~ "CR",
                              TRUE ~ category)) %>%
  select(Binomial, Year, SSP, Model, BIOME, category) %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

## predict future projections
mam_iucn_future <- list()
for(k in 1:24){
  mam_temp <- mam_data %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]]) %>% select(-c(Binomial, SSP, Year, Model))
  
  pred_future <- cbind(mam_temp %>% filter(area > 0),
                       st_coordinates(mam_all_sp %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])))
  
  rownames(pred_future) <- filter(mam_all_sp, Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])$Binomial
  
  future <- pred_ass(pred = pred_future[, -12],
                      biomes = biomes,
                      models = list(aa_thrt,
                                    aa_cr,
                                    aa_envu,
                                    aa_nthrt))
  
  mam_iucn_future[[k]] <- mam_data %>% filter(Year == years[[k]], SSP == ssps[[k]], Model == models[[k]]) %>%
    mutate(category = case_when(Binomial %in% future$lc ~ "LC",
                                Binomial %in% future$nt ~ "NT",
                                Binomial %in% future$vu ~ "VU",
                                Binomial %in% future$en ~ "EN",
                                Binomial %in% future$cr ~ "CR")) %>%
    select(Binomial, Year, SSP, Model, BIOME, category)
}
mam_iucn_future <- bind_rows(mam_iucn_future) %>% left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

mam_iucn <- bind_rows(mam_iucn_current, mam_iucn_future)

mam_accuracy <-
  lapply(lapply(list(aa_thrt, aa_cr, aa_envu, aa_nthrt),
                "[[", "mean_errors"),
         function(x)
           tibble(accuracy = x, partition = names(x), clade = "mammal"))

for(i in 1:4){
  mam_accuracy[[i]] <- mam_accuracy[[i]] %>%
    mutate(model = case_when(i == 1 ~ "Threatened/Non-threatened",
                             i == 2 ~ "CR/VE|EN",
                             i == 3 ~ "VE/EN",
                             i == 4 ~ "NT/LC"))
}

#### BIRDS ####
bird_sdm <- tet_sdm[which(tet_sdm %in% filter(tet_iucn, Class == "Bird")$Binomial)]

bird_current <- lapply(bird_sdm, function(x) read_sf(paste("output/SDM/Current/Shapefiles/", x, ".shp", sep = ""))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))
bird_future <- lapply(bird_sdm, function(x) read_sf(paste("output/SDM/Future/Shapefiles/", x, ".shp", sep = "")) %>%
                       mutate(SSP = sapply(str_split(Model, " "), "[[", 2),
                              Model = sapply(str_split(Model, " "), "[[", 1),
                              Year = sapply(str_split(Year, "-"), "[[", 2))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))

bird_current_biomes <- read_csv("output/tet_biomes.csv") %>%
  filter(Binomial %in% bird_current$Binomial)

bird_future_biomes <- bind_rows(lapply(1:24, function(x) read_csv(paste("output/tet_biomes_", x, ".csv", sep = "")))) %>%
  filter(Binomial %in% bird_current$Binomial)

# import birdmal data
EPBC <- read_csv("data/EPBC.csv") %>% mutate(Binomial = str_replace(Binomial, " ", "_"))

bird_data <- read_csv("output/tet_lu_hp.csv") %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Mass)) %>%
  left_join(tet_dist %>% as_tibble() %>% select(-geometry)) %>%
  left_join(EPBC %>% select(-Class)) %>%
  left_join(bind_rows(bird_current_biomes, bird_future_biomes) %>% select(-area)) %>%
  mutate(category = case_when(is.na(EPBC) ~ category,
                              TRUE ~ EPBC),
         category = case_when(Year == 2010 ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-EPBC) %>%
  filter(Binomial %in% bird_current$Binomial)

bird_data[which(is.na(bird_data[,"category"])), "category"] <- "NE"
bird_data[which(bird_data$Binomial == "Chelodina_novaeguineae"), "category"] <- "LC"

# drop species with extremely wrong predicted ranges
bird_data <- bird_data %>%
  left_join(tet_dist  %>%
              st_transform("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=61.55,-10.87,-40.19,-39.4924,-32.7221,-32.8979,-9.99400000001316 +units=m +no_defs +type=crs") %>%
              mutate(area_iucn = as.numeric(units::set_units(st_area(geometry), "km2"))) %>%
              as_tibble() %>%
              select(Binomial, area_iucn)) %>%
  mutate(pred_real = log(area/area_iucn),
         category = case_when(abs(pred_real) >= 3 & Year == 2010 ~ "NE",
                              TRUE ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-c(pred_real, area_iucn)) %>%
  drop_na(Mass)

# combine all data and omit species without ranges
bird_all <- bind_rows(bird_current, bird_future) %>%
  filter(Binomial %in% bird_data$Binomial) %>%
  mutate(id = str_replace(paste(Binomial, Model, SSP, Year, sep = "_"), "NA_NA_NA", "2020"),
         Model = case_when(is.na(Model) ~ "current",
                           TRUE ~ Model),
         Year = case_when(is.na(Year) ~ "current",
                          TRUE ~ Year),
         SSP = case_when(is.na(SSP) ~ "current",
                         TRUE ~ SSP))

bird_all_sp <- st_centroid(bird_all)

empty <- which(is.na(st_coordinates(bird_all_sp)))
bird_all_sp <- bird_all_sp[-empty, ]
bird_all <- bird_all[-empty, ]

# format for xgboost
pred <-
  cbind(bird_data %>% filter(Year == 2010) %>% select(-c(Binomial, SSP, Year, Model)),
        st_coordinates(bird_all_sp %>% filter(Year == "current")))

rownames(pred) <- filter(bird_all_sp, Year == "current")$Binomial

#### RUN AUTOMATED ASSESSMENT
biomes <- list(thrt = list(tropical = c(1, 7),
                           temperate = c(8, 4, 10),
                           desert = c(13),
                           med = c(12)),
               cr = list(tropical = c(1, 7),
                         non_tropical = c(8, 4, 10, 12, 13)),
               envu = list(tropical = c(1, 7),
                           non_tropical = c(8, 4, 10, 12, 13)),
               nthrt = list(tropical = c(1, 7),
                            temperate = c(8, 4, 10),
                            dry = c(12, 13)))

# seperate into assessed and non-assessed species
pred_ne_dd <- 
  pred[pred$category %in% c("NE", "DD"),]
pred <- 
  pred[!pred$category %in% c("NE", "DD", "EX", "EW"),]

aa_thrt <- auto_ass(pred = pred[, -c(12)],
                    biomes = biomes$thrt,
                    path = "output/AA/birds/",
                    cat = "threatened")

# predict specific categories
pred_thrt <- pred[pred$threatened == "yes", ]
pred_nthrt <- pred[pred$threatened == "no", ]

### CR vs VU+EN
cr_ind <- pred_thrt[pred_thrt$threatened == "yes", ]$category == "CR"

pred_thrt$cr <- 0
pred_thrt$cr[cr_ind] <- 1

aa_cr <- auto_ass(pred = pred_thrt[, -c(12, 14)],
                  biomes = biomes$cr,
                  path = "output/AA/birds/",
                  cat = "cr")


pred_thrt_envu <- pred_thrt[pred_thrt$cr == 0, ]

### VU vs EN
en_ind <- pred_thrt_envu[pred_thrt_envu$threatened == "yes", ]$category == "EN"

pred_thrt_envu$envu <- 0
pred_thrt_envu$envu[en_ind] <- 1

pred_thrt_envu <- 
  pred_thrt_envu[!names(pred_thrt_envu) %in% c("cr")]

aa_envu <- auto_ass(pred = pred_thrt_envu[, -c(12, 14)],
                    biomes = biomes$envu,
                    path = "output/AA/birds/",
                    cat = "envu")

### NT vs LC
nt_ind <- pred_nthrt[pred_nthrt$threatened == "no", ]$category == "NT"

pred_nthrt$nt <- 0
pred_nthrt$nt[nt_ind] <- 1

aa_nthrt <- auto_ass(pred = pred_nthrt[, -c(12, 14)],
                     biomes = biomes$nthrt,
                     path = "output/AA/birds/",
                     cat = "nt")


## predict NE & DD species
ne_dd <- pred_ass(pred = rbind(pred, pred_ne_dd)[, -12],
                   biomes = biomes,
                   models = list(aa_thrt,
                                 aa_cr,
                                 aa_envu,
                                 aa_nthrt))

bird_iucn_current <- bird_data %>% filter(Year == 2010) %>%
  mutate(category = case_when(Binomial %in% ne_dd$lc ~ "LC",
                              Binomial %in% ne_dd$nt ~ "NT",
                              Binomial %in% ne_dd$vu ~ "VU",
                              Binomial %in% ne_dd$en ~ "EN",
                              Binomial %in% ne_dd$cr ~ "CR",
                              TRUE ~ category)) %>%
  select(Binomial, Year, SSP, Model, BIOME, category) %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

## predict future projections
bird_iucn_future <- list()
for(k in 1:24){
  bird_temp <- bird_data %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]]) %>% select(-c(Binomial, SSP, Year, Model))
  
  pred_future <- cbind(bird_temp %>% filter(area > 0),
                       st_coordinates(bird_all_sp %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])))
  
  rownames(pred_future) <- filter(bird_all_sp, Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])$Binomial
  
  future <- pred_ass(pred = pred_future[, -12],
                      biomes = biomes,
                      models = list(aa_thrt,
                                    aa_cr,
                                    aa_envu,
                                    aa_nthrt))
  
  bird_iucn_future[[k]] <- bird_data %>% filter(Year == years[[k]], SSP == ssps[[k]], Model == models[[k]]) %>%
    mutate(category = case_when(Binomial %in% future$lc ~ "LC",
                                Binomial %in% future$nt ~ "NT",
                                Binomial %in% future$vu ~ "VU",
                                Binomial %in% future$en ~ "EN",
                                Binomial %in% future$cr ~ "CR")) %>%
    select(Binomial, Year, SSP, Model, BIOME, category)
}
bird_iucn_future <- bind_rows(bird_iucn_future) %>% left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

bird_iucn <- bind_rows(bird_iucn_current, bird_iucn_future)

bird_accuracy <-
  lapply(lapply(list(aa_thrt, aa_cr, aa_envu, aa_nthrt),
                "[[", "mean_errors"),
         function(x)
           tibble(accuracy = x, partition = names(x), clade = "bird"))

for(i in 1:4){
  bird_accuracy[[i]] <- bird_accuracy[[i]] %>%
    mutate(model = case_when(i == 1 ~ "Threatened/Non-threatened",
                             i == 2 ~ "CR/VE|EN",
                             i == 3 ~ "VE/EN",
                             i == 4 ~ "NT/LC"))
}

#### REPTILES ####
rep_sdm <- tet_sdm[which(tet_sdm %in% filter(tet_iucn, Class == "Reptile")$Binomial)]

rep_current <- lapply(rep_sdm, function(x) read_sf(paste("output/SDM/Current/Shapefiles/", x, ".shp", sep = ""))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))
rep_future <- lapply(rep_sdm, function(x) read_sf(paste("output/SDM/Future/Shapefiles/", x, ".shp", sep = "")) %>%
                       mutate(SSP = sapply(str_split(Model, " "), "[[", 2),
                              Model = sapply(str_split(Model, " "), "[[", 1),
                              Year = sapply(str_split(Year, "-"), "[[", 2))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))


rep_current_biomes <- read_csv("output/tet_biomes.csv") %>%
  filter(Binomial %in% rep_current$Binomial)

rep_future_biomes <- bind_rows(lapply(1:24, function(x) read_csv(paste("output/tet_biomes_", x, ".csv", sep = "")))) %>%
  filter(Binomial %in% rep_current$Binomial)

# import repmal data
EPBC <- read_csv("data/EPBC.csv") %>% mutate(Binomial = str_replace(Binomial, " ", "_"))

rep_data <- read_csv("output/tet_lu_hp.csv") %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Mass)) %>%
  left_join(tet_dist %>% as_tibble() %>% select(-geometry)) %>%
  left_join(EPBC %>% select(-Class)) %>%
  left_join(bind_rows(rep_current_biomes, rep_future_biomes) %>% select(-area)) %>%
  mutate(category = case_when(is.na(EPBC) ~ category,
                              TRUE ~ EPBC),
         category = case_when(Year == 2010 ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-EPBC) %>%
  filter(Binomial %in% rep_current$Binomial)

rep_data[which(is.na(rep_data[,"category"])), "category"] <- "NE"
rep_data[which(rep_data$Binomial == "Chelodina_novaeguineae"), "category"] <- "LC"

# drop species with extremely wrong predicted ranges
rep_data <- rep_data %>%
  left_join(tet_dist  %>%
              st_transform("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=61.55,-10.87,-40.19,-39.4924,-32.7221,-32.8979,-9.99400000001316 +units=m +no_defs +type=crs") %>%
              mutate(area_iucn = as.numeric(units::set_units(st_area(geometry), "km2"))) %>%
              as_tibble() %>%
              select(Binomial, area_iucn)) %>%
  mutate(pred_real = log(area/area_iucn),
         category = case_when(abs(pred_real) >= 3 & Year == 2010 ~ "NE",
                              TRUE ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-c(pred_real, area_iucn)) %>%
  drop_na(Mass)

# combine all data and omit species without ranges
rep_all <- bind_rows(rep_current, rep_future) %>%
  filter(Binomial %in% rep_data$Binomial) %>%
  mutate(id = str_replace(paste(Binomial, Model, SSP, Year, sep = "_"), "NA_NA_NA", "2020"),
         Model = case_when(is.na(Model) ~ "current",
                           TRUE ~ Model),
         Year = case_when(is.na(Year) ~ "current",
                          TRUE ~ Year),
         SSP = case_when(is.na(SSP) ~ "current",
                         TRUE ~ SSP))

rep_all_sp <- st_centroid(rep_all)

empty <- which(is.na(st_coordinates(rep_all_sp)))
rep_all_sp <- rep_all_sp[-empty, ]
rep_all <- rep_all[-empty, ]

# format for xgboost
pred <-
  cbind(rep_data %>% filter(Year == 2010) %>% select(-c(Binomial, SSP, Year, Model)),
        st_coordinates(rep_all_sp %>% filter(Year == "current")))

rownames(pred) <- filter(rep_all_sp, Year == "current")$Binomial

#### RUN AUTOMATED ASSESSMENT
biomes <- list(thrt = list(tropical = c(1, 7),
                           temperate = c(8, 4, 10),
                           desert = c(13),
                           med = c(12)),
               cr = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               envu = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               nthrt = list(all = c(1, 4, 7, 8, 10, 12, 13)))

# seperate into assessed and non-assessed species
pred_ne_dd <- 
  pred[pred$category %in% c("NE", "DD"),]
pred <- 
  pred[!pred$category %in% c("NE", "DD", "EX", "EW"),]

aa_thrt <- auto_ass(pred = pred[, -c(12)],
                    biomes = biomes$thrt,
                    path = "output/AA/reptiles/",
                    cat = "threatened")

# predict specific categories
pred_thrt <- pred[pred$threatened == "yes", ]
pred_nthrt <- pred[pred$threatened == "no", ]

### CR vs VU+EN
cr_ind <- pred_thrt[pred_thrt$threatened == "yes", ]$category == "CR"

pred_thrt$cr <- 0
pred_thrt$cr[cr_ind] <- 1

aa_cr <- auto_ass(pred = pred_thrt[, -c(12, 14)],
                  biomes = biomes$cr,
                  path = "output/AA/reptiles/",
                  cat = "cr",
                  partition = F)


pred_thrt_envu <- pred_thrt[pred_thrt$cr == 0, ]

### VU vs EN
en_ind <- pred_thrt_envu[pred_thrt_envu$threatened == "yes", ]$category == "EN"

pred_thrt_envu$envu <- 0
pred_thrt_envu$envu[en_ind] <- 1

pred_thrt_envu <- 
  pred_thrt_envu[!names(pred_thrt_envu) %in% c("cr")]

aa_envu <- auto_ass(pred = pred_thrt_envu[, -c(12, 14)],
                    biomes = biomes$envu,
                    path = "output/AA/reptiles/",
                    cat = "envu",
                    partition = F)

### NT vs LC
nt_ind <- pred_nthrt[pred_nthrt$threatened == "no", ]$category == "NT"

pred_nthrt$nt <- 0
pred_nthrt$nt[nt_ind] <- 1

aa_nthrt <- auto_ass(pred = pred_nthrt[, -c(12, 14)],
                     biomes = biomes$nthrt,
                     path = "output/AA/reptiles/",
                     cat = "nt")


## predict NE & DD species
ne_dd <- pred_ass(pred = rbind(pred, pred_ne_dd)[, -12],
                   biomes = biomes,
                   models = list(aa_thrt,
                                 aa_cr,
                                 aa_envu,
                                 aa_nthrt))

rep_iucn_current <- rep_data %>% filter(Year == 2010) %>%
  mutate(category = case_when(Binomial %in% ne_dd$lc ~ "LC",
                              Binomial %in% ne_dd$nt ~ "NT",
                              Binomial %in% ne_dd$vu ~ "VU",
                              Binomial %in% ne_dd$en ~ "EN",
                              Binomial %in% ne_dd$cr ~ "CR",
                              TRUE ~ category)) %>%
  select(Binomial, Year, SSP, Model, BIOME, category) %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

## predict future projections
rep_iucn_future <- list()
for(k in 1:24){
  rep_temp <- rep_data %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]]) %>% select(-c(Binomial, SSP, Year, Model))
  
  pred_future <- cbind(rep_temp %>% filter(area > 0),
                       st_coordinates(rep_all_sp %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])))
  
  rownames(pred_future) <- filter(rep_all_sp, Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])$Binomial
  
  future <- pred_ass(pred = pred_future[, -12],
                      biomes = biomes,
                      models = list(aa_thrt,
                                    aa_cr,
                                    aa_envu,
                                    aa_nthrt))
  
  rep_iucn_future[[k]] <- rep_data %>% filter(Year == years[[k]], SSP == ssps[[k]], Model == models[[k]]) %>%
    mutate(category = case_when(Binomial %in% future$lc ~ "LC",
                                Binomial %in% future$nt ~ "NT",
                                Binomial %in% future$vu ~ "VU",
                                Binomial %in% future$en ~ "EN",
                                Binomial %in% future$cr ~ "CR")) %>%
    select(Binomial, Year, SSP, Model, BIOME, category)
}
rep_iucn_future <- bind_rows(rep_iucn_future) %>% left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

rep_iucn <- bind_rows(rep_iucn_current, rep_iucn_future)

rep_accuracy <-
  lapply(lapply(list(aa_thrt, aa_cr, aa_envu, aa_nthrt),
                "[[", "mean_errors"),
         function(x)
           tibble(accuracy = x, partition = names(x), clade = "reptile"))

for(i in 1:4){
  rep_accuracy[[i]] <- rep_accuracy[[i]] %>%
    mutate(model = case_when(i == 1 ~ "Threatened/Non-threatened",
                             i == 2 ~ "CR/VE|EN",
                             i == 3 ~ "VE/EN",
                             i == 4 ~ "NT/LC"))
}

#### AMPHIBIANS ####
amph_sdm <- tet_sdm[which(tet_sdm %in% filter(tet_iucn, Class == "Amphibian")$Binomial)]

amph_current <- lapply(amph_sdm, function(x) read_sf(paste("output/SDM/Current/Shapefiles/", x, ".shp", sep = ""))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))
amph_future <- lapply(amph_sdm, function(x) read_sf(paste("output/SDM/Future/Shapefiles/", x, ".shp", sep = "")) %>%
                       mutate(SSP = sapply(str_split(Model, " "), "[[", 2),
                              Model = sapply(str_split(Model, " "), "[[", 1),
                              Year = sapply(str_split(Year, "-"), "[[", 2))) %>%
  bind_rows() %>%
  st_transform(st_crs(tet_dist))

amph_current_biomes <- read_csv("output/tet_biomes.csv") %>%
  filter(Binomial %in% amph_current$Binomial)

amph_future_biomes <- bind_rows(lapply(1:24, function(x) read_csv(paste("output/tet_biomes_", x, ".csv", sep = "")))) %>%
  filter(Binomial %in% amph_current$Binomial)

# import amphmal data
EPBC <- read_csv("data/EPBC.csv") %>% mutate(Binomial = str_replace(Binomial, " ", "_"))

amph_data <- read_csv("output/tet_lu_hp.csv") %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Mass)) %>%
  left_join(tet_dist %>% as_tibble() %>% select(-geometry)) %>%
  left_join(EPBC %>% select(-Class)) %>%
  left_join(bind_rows(amph_current_biomes, amph_future_biomes) %>% select(-area)) %>%
  mutate(category = case_when(is.na(EPBC) ~ category,
                              TRUE ~ EPBC),
         category = case_when(Year == 2010 ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-EPBC) %>%
  filter(Binomial %in% amph_current$Binomial)

amph_data[which(is.na(amph_data[,"category"])), "category"] <- "NE"
amph_data[which(amph_data$Binomial == "Chelodina_novaeguineae"), "category"] <- "LC"

# drop species with extremely wrong predicted ranges
amph_data <- amph_data %>%
  left_join(tet_dist  %>%
              st_transform("+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=61.55,-10.87,-40.19,-39.4924,-32.7221,-32.8979,-9.99400000001316 +units=m +no_defs +type=crs") %>%
              mutate(area_iucn = as.numeric(units::set_units(st_area(geometry), "km2"))) %>%
              as_tibble() %>%
              select(Binomial, area_iucn)) %>%
  mutate(pred_real = log(area/area_iucn),
         category = case_when(abs(pred_real) >= 3 & Year == 2010 ~ "NE",
                              TRUE ~ category),
         threatened = case_when(category %in% c("LC", "NT") ~ "no",
                                category %in% c("VU", "EN", "CR") ~ "yes")) %>%
  select(-c(pred_real, area_iucn)) %>%
  drop_na(Mass)

# combine all data and omit species without ranges
amph_all <- bind_rows(amph_current, amph_future) %>%
  filter(Binomial %in% amph_data$Binomial) %>%
  mutate(id = str_replace(paste(Binomial, Model, SSP, Year, sep = "_"), "NA_NA_NA", "2020"),
         Model = case_when(is.na(Model) ~ "current",
                           TRUE ~ Model),
         Year = case_when(is.na(Year) ~ "current",
                          TRUE ~ Year),
         SSP = case_when(is.na(SSP) ~ "current",
                         TRUE ~ SSP))

amph_all_sp <- st_centroid(amph_all)

empty <- which(is.na(st_coordinates(amph_all_sp)))
amph_all_sp <- amph_all_sp[-empty, ]
amph_all <- amph_all[-empty, ]

# format for xgboost
pred <-
  cbind(amph_data %>% filter(Year == 2010) %>% select(-c(Binomial, SSP, Year, Model)),
        st_coordinates(amph_all_sp %>% filter(Year == "current")))

rownames(pred) <- filter(amph_all_sp, Year == "current")$Binomial

#### RUN AUTOMATED ASSESSMENT
biomes <- list(thrt = list(tropical = c(1, 7),
                           non_tropical = c(8, 4, 10, 12, 13)),
               cr = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               envu = list(all = c(1, 4, 7, 8, 10, 12, 13)),
               nthrt = list(tropical = c(1, 7),
                            non_tropical = c(8, 4, 10, 12, 13)))

# seperate into assessed and non-assessed species
pred_ne_dd <- 
  pred[pred$category %in% c("NE", "DD"),]
pred <- 
  pred[!pred$category %in% c("NE", "DD", "EX", "EW"),]

aa_thrt <- auto_ass(pred = pred[, -c(12)],
                    biomes = biomes$thrt,
                    path = "output/AA/amphibians/",
                    cat = "threatened")

# predict specific categories
pred_thrt <- pred[pred$threatened == "yes", ]
pred_nthrt <- pred[pred$threatened == "no", ]

### CR vs VU+EN
cr_ind <- pred_thrt[pred_thrt$threatened == "yes", ]$category == "CR"

pred_thrt$cr <- 0
pred_thrt$cr[cr_ind] <- 1

aa_cr <- auto_ass(pred = pred_thrt[, -c(12, 14)],
                  biomes = biomes$cr,
                  path = "output/AA/amphibians/",
                  cat = "cr",
                  partition = F)


pred_thrt_envu <- pred_thrt[pred_thrt$cr == 0, ]

### VU vs EN
en_ind <- pred_thrt_envu[pred_thrt_envu$threatened == "yes", ]$category == "EN"

pred_thrt_envu$envu <- 0
pred_thrt_envu$envu[en_ind] <- 1

pred_thrt_envu <- 
  pred_thrt_envu[!names(pred_thrt_envu) %in% c("cr")]

aa_envu <- auto_ass(pred = pred_thrt_envu[, -c(12, 14)],
                    biomes = biomes$envu,
                    path = "output/AA/amphibians/",
                    cat = "envu",
                    partition = F)

### NT vs LC
nt_ind <- pred_nthrt[pred_nthrt$threatened == "no", ]$category == "NT"

pred_nthrt$nt <- 0
pred_nthrt$nt[nt_ind] <- 1

aa_nthrt <- auto_ass(pred = pred_nthrt[, -c(12, 14)],
                     biomes = biomes$nthrt,
                     path = "output/AA/amphibians/",
                     cat = "nt")


## predict NE & DD species
ne_dd <- pred_ass(pred = rbind(pred, pred_ne_dd)[, -12],
                   biomes = biomes,
                   models = list(aa_thrt,
                                 aa_cr,
                                 aa_envu,
                                 aa_nthrt))

amph_iucn_current <- amph_data %>% filter(Year == 2010) %>%
  mutate(category = case_when(Binomial %in% ne_dd$lc ~ "LC",
                              Binomial %in% ne_dd$nt ~ "NT",
                              Binomial %in% ne_dd$vu ~ "VU",
                              Binomial %in% ne_dd$en ~ "EN",
                              Binomial %in% ne_dd$cr ~ "CR",
                              TRUE ~ category)) %>%
  select(Binomial, Year, SSP, Model, BIOME, category) %>%
  left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

## predict future projections
amph_iucn_future <- list()
for(k in 1:24){
  amph_temp <- amph_data %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]]) %>% select(-c(Binomial, SSP, Year, Model))
  
  pred_future <- cbind(amph_temp %>% filter(area > 0),
                       st_coordinates(amph_all_sp %>% filter(Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])))
  
  rownames(pred_future) <- filter(amph_all_sp, Year == years[[k]], Model == models[[k]], SSP == ssps[[k]])$Binomial
  
  future <- pred_ass(pred = pred_future[, -12],
                      biomes = biomes,
                      models = list(aa_thrt,
                                    aa_cr,
                                    aa_envu,
                                    aa_nthrt))
  
  amph_iucn_future[[k]] <- amph_data %>% filter(Year == years[[k]], SSP == ssps[[k]], Model == models[[k]]) %>%
    mutate(category = case_when(Binomial %in% future$lc ~ "LC",
                                Binomial %in% future$nt ~ "NT",
                                Binomial %in% future$vu ~ "VU",
                                Binomial %in% future$en ~ "EN",
                                Binomial %in% future$cr ~ "CR")) %>%
    select(Binomial, Year, SSP, Model, BIOME, category)
}
amph_iucn_future <- bind_rows(amph_iucn_future) %>% left_join(read_csv("data/tet_data.csv") %>% select(Binomial, Class))

amph_iucn <- bind_rows(amph_iucn_current, amph_iucn_future)

amph_accuracy <-
  lapply(lapply(list(aa_thrt, aa_cr, aa_envu, aa_nthrt),
                "[[", "mean_errors"),
         function(x)
           tibble(accuracy = x, partition = names(x), clade = "amphibian"))

for(i in 1:4){
  amph_accuracy[[i]] <- amph_accuracy[[i]] %>%
    mutate(model = case_when(i == 1 ~ "Threatened/Non-threatened",
                             i == 2 ~ "CR/VE|EN",
                             i == 3 ~ "VE/EN",
                             i == 4 ~ "NT/LC"))
}