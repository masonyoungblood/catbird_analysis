
# SYLLABLE-LEVEL LOGISTICS ------------------------------------------------

#load libraries and set working directory
library(dplyr)
library(brms)
library(parallelDist)
library(parallel)
library(geosphere)
setwd("/Users/masonyoungblood/Documents/Work/Spring_2025/Catbirds/catbird_analysis")

#load in data
unit_df <- read.csv("data/unit_features.csv")
context <- read.csv("data/tweetynet_parsing.csv")
sex <- read.csv("data/sexing.csv")
data <- read.csv("data/clusters_features.csv")
locations <- read.csv("data/lat_lon.csv")

#reformat filenames
unit_df$filename_base <- unit_df$filename
context$basename <- gsub(".wav$", "", context$filename)

#create features object
features <- unit_df %>%
  left_join(context, by = c("filename_base" = "basename")) %>%
  rename(individual = ID) %>%
  left_join(sex %>% rename(individual = ColorID, sex = Sex), by = "individual") %>%
  mutate(context_simple = case_when(
    context == "Intrapair_Courtship" ~ "Intrapair",
    context == "Territorial" ~ "Territorial",
    context == "TerrResponse" ~ "TerrResponse",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(context_simple)) %>%
  mutate(context_simple = factor(context_simple, levels = c("Territorial", "Intrapair", "TerrResponse")))

#subset, z-score, individual to factor, sex to binary, and excursion and concavity divide by duration
features <- features[, c(5, 6, 7, 8, 9, 10, 12, 16, 42, 43)]
features$excursion <- features$excursion/features$duration
features$concavity <- features$concavity/features$duration
features <- cbind(scale(features[, 1:7]), features[, 8:10])
features$individual <- factor(features$individual)
features$sex <- ifelse(features$sex == "M", 1, 0)

#save features object
save(features, file = "proc/syl_features.RData")

#check for multicollinearity
performance::check_collinearity(lm(as.numeric(context_simple) ~ duration + min_freq + max_freq + bandwidth + concavity + excursion + entropy, data = features))

#restructure data for analysis
data$duration <- data$offset - data$onset
data$song <- sub("_[^_]+$", "", data[, which(colnames(data) == "unit")])

#add year and context and location data
data$context <- context$context[match(data$song, gsub(".wav", "", context$filename))]
data$individual <- context$ID[match(data$song, gsub(".wav", "", context$filename))]
data$year <- context$year[match(data$song, gsub(".wav", "", context$filename))]
locations <- aggregate(cbind(Longitude, Latitude) ~ ID, locations, FUN = mean)
locations$individual <- do.call(rbind, strsplit(locations$ID, "_"))[, 1]
locations$year <- do.call(rbind, strsplit(locations$ID, "_"))[, 2]
data$latitude <- NA
data$longitude <- NA
for(i in 1:nrow(locations)){
  data$latitude[which(data$individual == locations$individual[i] & data$year == locations$year[i])] <- locations$Latitude[i]
  data$longitude[which(data$individual == locations$individual[i] & data$year == locations$year[i])] <- locations$Longitude[i]
}

#sort data by song (using stringr order which handles numbers in strings correctly)
data <- data[stringr::str_order(data$unit, numeric = TRUE), ]

#remove songs without geographic info
data <- data[-which(is.na(data$latitude)), ]

#subset to assess (full dataset is WAY too large)
#not just for matrix calculation, for brms to analyze the data
set.seed(123)
inds <- sample(which(!is.na(data$latitude)), 5000, replace = FALSE)
sampled_data <- data[inds, ]

#get distance matrix
syllable_distances <- parDist(as.matrix(sampled_data[, grep("feat_", colnames(sampled_data))]), method = "cosine", threads = detectCores()-1)
syllable_distances <- as.matrix(syllable_distances)

#get other matrices
geo_dist_matrix <- distm(sampled_data[, c("longitude", "latitude")], fun = distHaversine)
year_dist_matrix <- abs(outer(sampled_data$year, sampled_data$year, "-"))

#identify final to store, masking syllables from same individual
#go ahead and compute different year mask for later
upper_tri_indices <- upper.tri(syllable_distances)
same_id_mask <- outer(sampled_data$individual, sampled_data$individual, "==")

#construct and save data object
syl_dist_data <- data.frame(
  syl_dist = syllable_distances[upper_tri_indices & !same_id_mask],
  year_dist = year_dist_matrix[upper_tri_indices & !same_id_mask],
  geo_dist = geo_dist_matrix[upper_tri_indices & !same_id_mask]
)
save(syl_dist_data, file = "proc/syl_dist_data.RData")

#do the same except getting syllable distances within individuals across years
matrix_indices <- which(upper_tri_indices & same_id_mask, arr.ind = TRUE)
row_indices <- matrix_indices[, "row"]
col_indices <- matrix_indices[, "col"]
syl_dist_w_in_data <- data.frame(
  individual = sampled_data$individual[row_indices],
  year1 = sampled_data$year[row_indices],
  year2 = sampled_data$year[col_indices],
  syl_dist = syllable_distances[upper_tri_indices & same_id_mask],
  year_dist = year_dist_matrix[upper_tri_indices & same_id_mask]
)
save(syl_dist_w_in_data, file = "proc/syl_dist_w_in_data.RData")

# CONTEXT MODEL -----------------------------------------------------------

#load data and libraries
library(brms)
load("proc_data/syl_features.RData")

#set priors for more efficient sampling
syl_context_prior <- c(
  prior(normal(0, 2), class = "Intercept", dpar = "muIntrapair"),
  prior(normal(0, 2), class = "Intercept", dpar = "muTerrResponse"),
  prior(normal(0, 0.5), class = "b", dpar = "muIntrapair"),
  prior(normal(0, 0.5), class = "b", dpar = "muTerrResponse"),
  prior(normal(0, 1), lb = 0, class = "sd", dpar = "muIntrapair"),
  prior(normal(0, 1), lb = 0, class = "sd", dpar = "muTerrResponse")
)

#fit model
#individual included as random effect, because context varies plenty across it
#note that context and sex cannot appear as multivariate outcomes, because mv brms is only for gaussian and student dists
syl_context_model <- brm(
  context_simple ~ duration + min_freq + max_freq + bandwidth + concavity + excursion + entropy + (1|individual),
  data = features, prior = syl_context_prior, family = categorical(),
  #backend = "cmdstanr", 
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(syl_context_model, file = "models/syl_features_context.RData")

# SEX MODEL ---------------------------------------------------------------

#load data and libraries
library(brms)
load("proc/syl_features.RData")

#set priors for more efficient sampling
syl_sex_prior <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b")
)

#fit model
#not including individual as random effect, because only two females (extreme imbalance at group level)
syl_sex_model <- brm(
  sex ~ duration + min_freq + max_freq + bandwidth + concavity + excursion + entropy,
  data = features, prior = syl_sex_prior, family = bernoulli(),
  #backend = "cmdstanr", 
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(syl_sex_model, file = "models/syl_features_sex.RData")

# DISTANCE MODELS ---------------------------------------------------------

#load data and libraries
library(brms)
load("proc/syl_dist_data.RData")
load("proc/syl_dist_w_in_data.RData")

#set distance priors
syl_dist_prior <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b")
)

#run main distance model
syl_dist_model <- brm(
  scale(syl_dist) ~ scale(year_dist) + scale(geo_dist),
  data = syl_dist_data, prior = syl_dist_prior, family = gaussian(),
  #backend = "cmdstanr",
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(syl_dist_model, file = "models/syl_dists.RData")

#run model of distances within individuals model
syl_dist_w_in_model <- brm(
  scale(syl_dist) ~ scale(year_dist) + (1|individual),
  data = syl_dist_w_in_data, prior = syl_dist_prior, family = gaussian(),
  #backend = "cmdstanr",
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(syl_dist_w_in_model, file = "models/syl_dists_w_in.RData")

# SONG-LEVEL LOGISTICS ----------------------------------------------------

#load libraries and set working directory
library(dplyr)
library(brms)
library(parallelDist)
library(parallel)
library(geosphere)
library(dtw)
setwd("/Users/masonyoungblood/Documents/Work/Spring_2025/Catbirds/catbird_analysis")

#load in data
context <- read.csv("data/tweetynet_parsing.csv")
sex <- read.csv("data/sexing.csv")
data <- read.csv("data/clusters_features.csv")
locations <- read.csv("data/lat_lon.csv")

#restructure data for analysis
data$duration  <- data$offset - data$onset
data$song <- sub("_[^_]+$", "", data[, which(colnames(data) == "unit")])

#add year and context and location data
data$context <- context$context[match(data$song, gsub(".wav", "", context$filename))]
data$song_duration <- context$length_s[match(data$song, gsub(".wav", "", context$filename))]
data$individual <- context$ID[match(data$song, gsub(".wav", "", context$filename))]
data$year <- context$year[match(data$song, gsub(".wav", "", context$filename))]
locations <- aggregate(cbind(Longitude, Latitude) ~ ID, locations, FUN = mean)
locations$individual <- do.call(rbind, strsplit(locations$ID, "_"))[, 1]
locations$year <- do.call(rbind, strsplit(locations$ID, "_"))[, 2]
data$latitude <- NA
data$longitude <- NA
for(i in 1:nrow(locations)){
  data$latitude[which(data$individual == locations$individual[i] & data$year == locations$year[i])] <- locations$Latitude[i]
  data$longitude[which(data$individual == locations$individual[i] & data$year == locations$year[i])] <- locations$Longitude[i]
}

#sort data by song (using stringr order which handles numbers in strings correctly)
data <- data[stringr::str_order(data$unit, numeric = TRUE), ]

#remove songs without geographic info
data <- data[-which(is.na(data$latitude)), ]

#compute path lengths of songs for context and sex comparison
songs <- unique(data$song)
context_data <- do.call(rbind, mclapply(songs, function(x){
  indices <- which(data$song == x)
  if(length(indices) > 1){
    temp <- data[indices, grep("feat", colnames(data))]
    data.frame(
      path = sum(sapply(2:length(indices), function(y){dist(rbind(as.numeric(temp[y, ]), as.numeric(temp[y-1, ])), method = "cosine")}))/length(indices),
      length = length(indices),
      rate = length(indices)/data$song_duration[indices[1]],
      context = data$context[indices[1]],
      individual = data$individual[indices[1]]
    )
  } else{
    data.frame(
      path = NA,
      length = length(indices),
      rate = length(indices)/data$song_duration[indices[1]],
      context = data$context[indices[1]],
      individual = data$individual[indices[1]]
    )
  }
}, mc.cores = detectCores()-1))
context_data$sex <- ifelse(sex$Sex[match(context_data$individual, sex$ColorID)] == "M", 1, 0)
save(context_data, file = "proc_data/song_contexts.RData")

#get metadata object, and split features up into a matrix per song, as a list
first_song_indices <- !duplicated(data$song)
song_metadata <- data[first_song_indices, c("song", "individual", "year", "latitude", "longitude")]
feature_cols <- grep("feat_", colnames(data))
feature_data <- data[, feature_cols]
song_features_list <- split(feature_data, data$song)
song_features_list <- lapply(song_features_list, as.matrix)

#create all unique pairs of sampled songs
song_pairs <- combn(song_metadata$song, 2, simplify = FALSE)

#run dtw in parallel and assemble matrix
song_dists <- mclapply(song_pairs, function(pair){
  dtw(song_features_list[[pair[1]]], song_features_list[[pair[2]]], dist.method = "cosine")$distance
}, mc.cores = detectCores() - 1)
song_dists <- unlist(song_dists)
song_dists_matrix <- matrix(NA, nrow = nrow(song_metadata), ncol = nrow(song_metadata))
rownames(song_dists_matrix) <- colnames(song_dists_matrix) <- song_metadata$song
song_dists_matrix[upper.tri(song_dists_matrix)] <- song_dists
song_dists_matrix[lower.tri(song_dists_matrix)] <- t(song_dists_matrix)[lower.tri(song_dists_matrix)]
diag(song_dists_matrix) <- 0

#create other matrices
geo_dist_matrix <- distm(song_metadata[, c("longitude", "latitude")], fun = distHaversine)
year_dist_matrix <- abs(outer(song_metadata$year, song_metadata$year, "-"))

#create masks for filtering
upper_tri_indices <- upper.tri(song_dists_matrix)
same_id_mask <- outer(song_metadata$individual, song_metadata$individual, "==")

#create data for songs from different individuals
song_dist_data <- data.frame(
  song_dist = song_dists_matrix[upper_tri_indices & !same_id_mask],
  year_dist = year_dist_matrix[upper_tri_indices & !same_id_mask],
  geo_dist = geo_dist_matrix[upper_tri_indices & !same_id_mask]
)
save(song_dist_data, file = "proc/song_dist_data.RData")

#create data for songs within individuals
matrix_indices <- which(upper_tri_indices & same_id_mask, arr.ind = TRUE)
row_indices <- matrix_indices[, "row"]
col_indices <- matrix_indices[, "col"]
song_dist_w_in_data <- data.frame(
  individual = song_metadata$individual[row_indices],
  year1 = song_metadata$year[row_indices],
  year2 = song_metadata$year[col_indices],
  song_dist = song_dists_matrix[upper_tri_indices & same_id_mask],
  year_dist = year_dist_matrix[upper_tri_indices & same_id_mask]
)
save(song_dist_w_in_data, file = "proc/song_dist_w_in_data.RData")

# CONTEXT MODEL -----------------------------------------------------------

#load data and libraries
library(brms)
load("proc/song_contexts.RData")

#convert to factor
context_data$context[which(context_data$context == "Intrapair_Courtship")] <- "Intrapair"
context_data$context <- factor(context_data$context, levels = c("Territorial", "Intrapair", "TerrResponse"))

#set priors for more efficient sampling
song_context_prior <- c(
  prior(normal(0, 2), class = "Intercept", dpar = "muIntrapair"),
  prior(normal(0, 2), class = "Intercept", dpar = "muTerrResponse"),
  prior(normal(0, 0.5), class = "b", dpar = "muIntrapair"),
  prior(normal(0, 0.5), class = "b", dpar = "muTerrResponse"),
  prior(normal(0, 1), lb = 0, class = "sd", dpar = "muIntrapair"),
  prior(normal(0, 1), lb = 0, class = "sd", dpar = "muTerrResponse")
)

#fit model
#individual included as random effect, because context varies plenty across it
#note that context and sex cannot appear as multivariate outcomes, because mv brms is only for gaussian and student dists
song_context_model <- brm(
  context ~ scale(path) + scale(length) + scale(rate) + (1|individual),
  data = context_data, prior = song_context_prior, family = categorical(),
  #backend = "cmdstanr", 
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(song_context_model, file = "models/song_context.RData")

# SEX MODEL ---------------------------------------------------------------

#load data and libraries
library(brms)
load("proc/song_contexts.RData")

#set priors for more efficient sampling
song_sex_prior <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b")
)

#fit model
#not including individual as random effect, because only two females (extreme imbalance at group level)
#note that context and sex cannot appear as multivariate outcomes, because mv brms is only for gaussian and student dists
song_sex_model <- brm(
  sex ~ scale(path) + scale(length) + scale(rate),
  data = context_data, prior = song_sex_prior, family = bernoulli(),
  #backend = "cmdstanr", 
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(song_sex_model, file = "models/song_sex.RData")

# DISTANCE MODELS ---------------------------------------------------------

#load data and libraries
library(brms)
load("proc/song_dist_data.RData")
load("proc/song_dist_w_in_data.RData")

#set distance priors
song_dist_prior <- c(
  prior(normal(0, 2), class = "Intercept"),
  prior(normal(0, 0.5), class = "b")
)

#run main distance model
song_dist_model <- brm(
  scale(log(song_dist)) ~ scale(year_dist) + scale(geo_dist),
  data = song_dist_data, prior = song_dist_prior, family = gaussian(),
  #backend = "cmdstanr",
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(song_dist_model, file = "models/song_dists.RData")

#run model of distances within individuals model
song_dist_w_in_model <- brm(
  scale(log(song_dist)) ~ scale(year_dist) + (1|individual),
  data = song_dist_w_in_data, prior = song_dist_prior, family = gaussian(),
  #backend = "cmdstanr",
  #algorithm = "laplace"
  iter = 2000, chains = 4, cores = 4, threads = threading(3)
)

#save model
save(song_dist_w_in_model, file = "models/song_dists_w_in.RData")

#example posterior predictive check of songs distances with x number of years of distance
hist(rowMeans(posterior_predict(song_dist_w_in_model, newdata = data.frame(year_dist = 2, individual = unique(song_dist_w_in_data$individual)))), breaks = 100)
