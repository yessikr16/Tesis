#                   Species Distribution Modelling Script 1

#                        National University of Loja                       
# Tropical research center for environment and biodiversity (CITIAB)          

# Authors:  - M. Sc. Erick Angamarca (UNL)
#           - M. Sc. Juan Maita (UNL)

## Topics

#    Download gbif data
#    Filter data
#    Calibration areas
#    Masking variables
#    Jackknife analysis 
#    Correlation analysis

## 1. Libraries   ##############################################################

pacman::p_load(raster, sf, kuenm, tidyverse, datasets, sp, geodata,
  data.table, terra, grinnell, stringr, ellipsenm)

## 2. Variables   ##############################################################

#LIMPIAR AREA DE TRABAJO
rm(list = ls(all.names = TRUE))&cat("\014")&graphics.off()&gc()

#--------------------------------------------------------------editing variables
name_sp <- "Polylepis incana" #scientific name
umbral_min <- 3000  # minimum
umbral_max <- 3700  # maximum altitude range

#---------------------------------------------------------------------no editing
split_name <- str_split(name_sp, " ", n = 2)
genus <- split_name[[1]][1]
species <- split_name[[1]][2]
nm_sp <- paste0(substring(genus, 1, 3), "_", substring(species, 1, 3))
dis_het <- "_5_1km" #herogeneity distance
mxpath <- "C:/maxent" #maxent
mask <- "countries"

## 3. Download gbif data  ######################################################

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))
plot(mask_file$geometry)

#download
occ <- sp_occurrence(genus = genus, species = species, ext = mask_file, 
                     geo = T, download = T, fixnames = F, 
                     args = c("occurrenceStatus=PRESENT"))

#folder
sapply("4_data_csv/paso_1", function(x) if (!dir.exists(x)) 
  dir.create(x, recursive = T))

#save
write.csv(occ, paste0("4_data_csv/paso_1/", nm_sp, "_gbif", ".csv"), 
          row.names = F)

## 4. Relocation data ##########################################################

#bndb
file.copy(from = paste0("0_bndb_temp/", nm_sp, "_bndb.csv"), 
          to = paste0("4_data_csv/paso_1/", nm_sp, "_bndb", ".csv"))
#file.remove(paste0(nm_sp, "_bndb", ".csv"))

## 5. Union data base ##########################################################

#load data
df1 <- read.csv(paste0("4_data_csv/paso_1/", nm_sp, "_gbif", ".csv"))
df2 <- read.csv(paste0("4_data_csv/paso_1/", nm_sp, "_bndb", ".csv"))

#rename columns
df1 <- setnames(df1, old = c("scientificName", 'decimalLongitude',
                             'decimalLatitude'), new = c('species', 'longitude',
                                                         'latitude'))
df2 <- setnames(df2, old = c("scientificName", 'decimalLongitude',
                             'decimalLatitude'), new = c('species', 'longitude',
                                                         'latitude'))

df1 <- df1[ , c('species','longitude','latitude')]
df2 <- df2[ , c('species', 'longitude','latitude')]

#union data base
df_merge <- rbind(df1, df2)

#species name
df_merge$species <- name_sp

#plot
plot(mask_file$geometry)
points(df_merge[, 2:3])

#save
sapply("4_data_csv/paso_2", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))
write.csv(df_merge, paste0("4_data_csv/paso_2/", nm_sp, "_merge", ".csv"), 
          row.names = F)

## 6. Filter data  #############################################################

#load data
occ <- read.csv(paste0("4_data_csv/paso_2/", nm_sp, "_merge", ".csv"), 
                header = T, sep = ",", dec = ".") 

#without coordinates filter
occ_1 <- occ[!is.na(occ$longitude) & !is.na(occ$latitude), ] 

#duplicates filter
occ_1$code <-  paste(occ_1$species, occ_1$longitude, 
                   occ_1$latitude, sep = "_")  
occ_2 <- occ_1[!duplicated(occ_1$code), 1:4] 

#zero filter
occ_3 <- occ_2[occ_2$longitude != 0 & occ_2$latitude != 0, 1:3]

#folder
sapply("4_data_csv/paso_3", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#save
write.csv(occ_3, paste0("4_data_csv/paso_3/", nm_sp, "_filt", ".csv"), 
          row.names = FALSE)  

## 7. Altitude range  filter ################################################

#load data
occ_filt <- read.csv(paste0("4_data_csv/paso_3/", nm_sp, "_filt", ".csv"), 
                     header = T, sep = ",", dec = ".") 
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#shape points
spatial_pts <- st_as_sf(occ_filt, coords = c("longitude","latitude"), 
                        crs = st_crs(4326))

#intersection
spatial_pts <- suppressWarnings(st_intersection(spatial_pts, mask_file))
spatial_pts$ID <- NULL
plot(mask_file$geometry)
plot(spatial_pts, add = T)

#folder
sapply("2_vector/registros_sp", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#save
st_write(spatial_pts, paste0("2_vector/registros_sp/", nm_sp, "_filt.shp"), 
         driver = "ESRI Shapefile", delete_layer = T)

#load data
spatial_pts_filt <- st_read(paste0("2_vector/registros_sp/", nm_sp, "_filt.shp"))
alt <- rast("3_raster/wc_alt_countries.tif")

#extract altitude information
data <- data.frame(spatial_pts_filt$species, st_coordinates(spatial_pts_filt),
                   terra::extract(alt, spatial_pts_filt, ID = F))

#update colum names
names(data) <- c("species", "longitude", "latitude", "alt")
names(data)

#remove NA values
data <- na.omit(data)

#plot
plot(alt)
points(data[, 2:3], col = "red", pch = 20)

#print
data %>% arrange(desc(alt)) %>% head(30)
data %>% arrange(desc(alt)) %>% tail(30) 

#remove atipic data
data_umb_min <- data[data$alt > umbral_min, ]
data_umb_max <- data_umb_min[data_umb_min$alt < umbral_max, ]

#folder
sapply("4_data_csv/boxplot/", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#save
png(paste0("4_data_csv/boxplot/", nm_sp, "_bxp_alt.png"), width = 720, 
    height = 400, units = "px")
boxplot(data$alt, horizontal=T)
stripchart(data$alt, method = "jitter", pch = 1, add = T, col = "blue")
stripchart(data_umb_max$alt, method = "jitter", pch = 20, add = T, col = "red")
dev.off()

#plot
boxplot(data$alt, horizontal=T)
stripchart(data$alt, method = "jitter", pch = 1, add = T, col = "blue")
stripchart(data_umb_max$alt, method = "jitter", pch = 20, add = T, col = "red")

#folder
sapply("4_data_csv/paso_4", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#save
write.csv(data_umb_max, paste0("4_data_csv/paso_4/", nm_sp, "_alt", ".csv"), 
          row.names = F)

## 8. Data to heterogeneity climate analysis  ##################################

# Note: check projection data

#load data
data_alt_csv <- read.csv(paste0("4_data_csv/paso_4/", nm_sp, "_alt", ".csv"), 
                          header = T, sep = ",", dec = ".")

#shape points
data_alt_sf <- st_as_sf(data_alt_csv, coords = c("longitude", "latitude"), 
                        crs = st_crs(4326), remove = F)

#save
st_write(data_alt_sf, paste0("2_vector/registros_sp/", nm_sp, "_alt", ".shp"),
         driver = "ESRI Shapefile", delete_layer = T)

#folders
sapply(paste0("5_heterogeneidad/paso_1"), function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))
sapply(paste0("5_heterogeneidad/paso_2"), function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))
sapply(paste0("5_heterogeneidad/", nm_sp, dis_het), function(x)if
       (!dir.exists(x)) dir.create(x, recursive = T))

## 9. Heterogeneity climate analysis ######################################

## ArcMap and SDMTools

## 10. Calibration area  ###################################################

#folders
sapply("4_data_csv/paso_5", function(x)if(!dir.exists(x))
  dir.create(x, recursive = T))
sapply(paste0("6_calibracion/",  nm_sp, dis_het), function(x)if(!dir.exists(x))
  dir.create(x, recursive = T))

#directory
directory <- paste0("6_calibracion/",  nm_sp, dis_het, "/m_grinnell")

#prepare csv
occ_het <- read.csv(paste0("5_heterogeneidad/", nm_sp, dis_het, "/", 
                           nm_sp, dis_het,  "_rarefied_points.csv"), 
                    header = T)
occ_het$alt <- NULL

#save
write.csv(occ_het, paste0("4_data_csv/paso_5/", nm_sp, dis_het,
                          "_rarefied_points.csv"), row.names = F)

#load data
bioclim_current <- rast(list.files(path = paste0("3_raster/var_cur_2.5"), 
                                   pattern = ".asc$", 
                                   full.names = T))
bioclim_lgm <- rast(list.files(path = paste0("3_raster/var_lgm_2.5"),
                               pattern = ".asc$", full.names = T))

occ_het <- read.csv(paste0("4_data_csv/paso_5/", nm_sp, dis_het, 
                           "_rarefied_points.csv"), header = T)

#check extent
#ext(bioclim_current)
#ext(bioclim_lgm)
#ext(bioclim_current) == ext(bioclim_lgm)
#ext(bioclim_current) <- ext(bioclim_lgm)
#ext(bioclim_current) == ext(bioclim_lgm)

#simulation
help("M_simulationR")
M_simulationR(occ_het, current_variables = bioclim_current, project = T,
              projection_variables = bioclim_lgm, dispersal_kernel = "normal", 
              kernel_spread = 2, max_dispersers = 2, replicates = 10, 
              dispersal_events =10, simulation_period = 70, stable_lgm = 25, 
              transition_to_lgm = 10, lgm_to_current = 10, stable_current = 25, 
              scenario_span = 1, output_directory = directory, scale = T, 
              center = T, overwrite = T)

## 11. M variables mask   #############################################

#load data
bioclim_30s <- rast(list.files(paste0("3_raster/var_bio_30s"), 
                               pattern = ".asc$", full.names = T))
M_grinnell <- st_read(paste0("6_calibracion/", nm_sp, dis_het, "/m_grinnell/",
                             "accessible_area_M.shp"))

#mask
bioclim_mask <- mask(crop(bioclim_30s, M_grinnell), M_grinnell)

#plot
plot(bioclim_mask[[1]])
plot(st_geometry(M_grinnell), add = T)

#folder
suppressWarnings(dir.create(paste0("6_calibracion/", nm_sp, dis_het, "/mask_var"), 
                            recursive = T))

#directories
names_30s <- paste0("6_calibracion/", nm_sp, dis_het, "/mask_var/", 
                    names(bioclim_mask), ".asc")

#spat to stack
bioclim_mask <- stack(bioclim_mask)

#save
wr <- lapply(1:nlayers(bioclim_mask), function(x) {
  writeRaster(bioclim_mask[[x]], filename = names_30s[x], format = "ascii", 
              overwrite=T)
})

## 12. Variables selection   ###################################################

#load data
occ_het <- read.csv(paste0("5_heterogeneidad/", nm_sp, dis_het, "/", 
  nm_sp, dis_het, "_rarefied_points.csv"), header = T)

bioclim_mask <- stack(list.files(paste0("6_calibracion/", nm_sp, dis_het,
  "/mask_var"), pattern = ".asc$", full.names = T))

#jackknife test
bioclim_cont <- explore_var_contrib(occ = occ_het, M_variables = bioclim_mask,
  maxent.path = mxpath, plot = F, max.memory = 1200)
write.csv(bioclim_cont$Jackknife_results$Training_gain_with_without, 
  paste0("6_calibracion/", nm_sp, dis_het,"/jackknife_contrib", ".csv"),
  row.names = F)

#save plot
png(paste0("6_calibracion/", nm_sp, dis_het,"/jackknife_", nm_sp,
  ".png"), width = 560, height = 560, units = "px")
plot <- plot_contribution(bioclim_cont, col.cont = "gray25", col.imp = "gray25",
  col.with = "blue3", col.without = "cyan3", col.all = "black")
dev.off()

#correlation test
png(paste0("6_calibracion/", nm_sp, dis_het, "/correlation", "_", nm_sp,  
  ".png"), width = 510, height = 510, units = "px")
cor <- variable_correlation(bioclim_mask, correlation_limit = 0.8, corrplot = T,
  magnify_to = 4, save = F)
dev.off()

#bioclimatic variables selected
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
numb_vars <- c("01", "02", "03", "13","15") # change values  <<<<<<<<<<<
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#data frame
df <- data.frame(bioclim = paste0("bio_", numb_vars))
vars_select <- df %>%
  mutate(order = as.numeric(sub("bio_", "", bioclim))) %>%
  mutate(order = ifelse(order >= 10, order - 2, order))
print(vars_select)

#save
write.csv(vars_select, paste0("6_calibracion/", nm_sp, dis_het, 
  "/var_select", ".csv"), row.names = F)


