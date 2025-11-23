#              Download data to species distribution modelling

#                   Species Distribution Modelling Script 2

#                        National University of Loja                       
# Tropical research center for environment and biodiversity (CITIAB)          

# Authors:  - M. Sc. Erick Angamarca (UNL)
#           - M. Sc. Juan Maita (UNL)

## Topics

#    Install packages
#    Download Ecuador Boundaries
#    Download altitude raster
#    Download bio climatic variables 30''
#    Download bio climatic variables 2.5'

## 1. Install packages #########################################################

install.packages("raster", dependencies = T)
install.packages("rmapshaper", dependencies = T)
install.packages("sf", dependencies = T)
install.packages("sp", dependencies = T)
install.packages("terra", dependencies = T)
install.packages("tidyverse", dependencies = T)
install.packages("pacman", dependencies = T)
install.packages("remotes", dependencies = T)
install.packages("devtools", dependencies = T)
install.packages("nngeo", dependencies = T)
install.packages("rpaleoclim", dependencies = T)
install.packages("ggspatial", dependencies = T)
install.packages("scales", dependencies = T)
devtools::install_github("rspatial/geodata", force = T, dependencies = T)
devtools::install_github("marlonecobos/kuenm", force = T, dependencies = T)
devtools::install_github("marlonecobos/ellipsenm", force = T, dependencies = F)
devtools::install_github("fmachados/grinnell", force = T, dependencies = T)
remotes::install_github("rsh249/vegdistmod", force = T)
remotes::install_github("r-spatial/qgisprocess")
remotes::install_version("rgeos", version = "0.6-2")
remotes::install_version("rgdal", version = "1.5-27", dependencies = T)
remotes::install_github("marlonecobos/evniche", dependencies = T)

## 2. Libraries ###############################################################

pacman::p_load(raster, sf, kuenm, tidyverse, datasets, geodata, nngeo, 
               terra, vegdistmod, rpaleoclim, ellipsenm)

## 3. Variables   ##############################################################

#clean workspace
rm(list = ls(all.names = T))&cat("\014")&graphics.off()&gc()

#variables
countries_cod <- c("EC", "CO", "PE")
mask <- "countries" # countries or ecuador

## 4. Ecuador Boundaries #######################################################

#link
link_dpa_prv <- paste0("https://www.ecuadorencifras.gob.ec//documentos/web-inec",
                       "/Cartografia/Clasificador_Geografico/2012/SHP.zip")

#folder
sapply("2_vector/INEC",function(x)if(!dir.exists(x))dir.create(x, recursive=T))

#download
{
  options(timeout = max(1000, getOption("timeout")))
  if (!file.exists(file.path("2_vector/INEC", "SHP.zip"))) { 
    download.file(link_dpa_prv, destfile = file.path("2_vector/INEC", "SHP.zip")) 
    unzip(file.path("2_vector/INEC", "SHP.zip"), exdir = "2_vector/INEC")}
}

#load data
ecu <- st_read("2_vector/INEC/SHP/nxprovincias.shp")

#remove Galapagos
ecu_filt <- ecu[!(ecu$DPA_DESPRO == "GALAPAGOS"),]

#dissolve
ecu_diss <- ecu_filt %>% group_by() %>% summarize()%>% mutate(ID=0) 

#buffer
ecu_buffer <- st_buffer(ecu_diss, 5000, endCapStyle="ROUND") %>% 
  st_transform(crs = 4326)
plot(ecu_buffer$geometry)

#fix shape file
st_is_valid(ecu_diss, reason = T)
ecu_fix <- st_make_valid(ecu_diss) %>% st_transform(crs = 4326)
st_is_valid(ecu_fix, reason = T)

#fix ecu buffer
st_is_valid(ecu_buffer, reason = T)
ecu_fix_buffer <- st_make_valid(ecu_buffer)
st_is_valid(ecu_fix_buffer, reason = T)

#plot
plot(st_geometry(ecu_fix))
plot(st_geometry(ecu_fix_buffer), add = T)

#save
st_write(ecu_fix, dsn = "2_vector/", layer = "ecuador_diss_4326.shp",
         driver = "ESRI Shapefile", delete_layer = T)
st_write(ecu_fix_buffer, dsn = "2_vector/", layer = "ecuador_buff_4326.shp",
         driver = "ESRI Shapefile", delete_layer = T)

## 5. GADM Boundaries ##########################################################

#folder
sapply("2_vector/", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#download
countries <- gadm(countries_cod, level=1, path = "2_vector/", 
                  version="latest", resolution=1)
countries
plot(countries)

#remove island
countries <- countries[!(countries$NAME_1 == "Galápagos") 
                       & !(countries$NAME_1 == "San Andrés y Providencia"), ]
plot(countries)

#dissolve
countries_sf <- st_as_sf(countries)
countries_diss <- countries_sf %>% group_by() %>% summarize()%>% mutate(ID=0) 
plot(countries_diss$geometry)

#projection
countries_utm <- st_transform(countries_diss, crs = 32717)

#buffer
countries_buff <- st_buffer(countries_utm, 5000, endCapStyle="ROUND") %>%
  st_transform(crs = 4326)
plot(countries_buff$geometry)

#fix shape file
st_is_valid(countries_diss, reason = T)
diss_fix <- st_make_valid(countries_diss)
st_is_valid(diss_fix, reason = T)

#fix buffer
st_is_valid(countries_buff, reason = T)
buff_fix <- st_make_valid(countries_buff)
st_is_valid(buff_fix, reason = T)

#plot
plot(diss_fix$geometry)
plot(buff_fix$geometry, add = T)

#save
st_write(diss_fix, dsn = "2_vector/", layer = "countries_diss_4326.shp",
         driver = "ESRI Shapefile", delete_layer = T)
st_write(buff_fix, dsn = "2_vector/", layer = "countries_buff_4326.shp",
         driver = "ESRI Shapefile", delete_layer = T)

## 6. Altitude    ###############################################################

#folder
sapply("3_raster/download/", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#Download
alt_tile1 <- worldclim_tile(var = "elev", -75, -10, 
                            path = "3_raster/download/", version = "2.1")
alt_tile2 <- worldclim_tile(var = "elev", -75, 10, 
                            path = "3_raster/download/", version = "2.1")

#join
alt_merge <- merge(alt_tile1, alt_tile2)
plot(alt_merge)

#mask
alt_mask <- mask(crop(alt_merge, mask_file), mask_file)

#plot
plot(alt_mask)
plot(mask_file, color = "", border ="red", add = T)

#save
writeRaster(alt_mask, filename = "3_raster/wc_alt_countries.tif", 
            filetype= "GTiff", overwrite=T)

## 7. Download current 0.5 min variables    ####################################

#load data
mask_buff <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#download tiles
bio_30s_tile1 <- worldclim_tile(var = "bio", -75, -10, path = "3_raster/download/", version = "2.1")
bio_30s_tile2 <- worldclim_tile(var = "bio", -75, 10, path = "3_raster/download/", version = "2.1")

#load tiles
bio_30s_tile1 <- rast(paste0("3_raster/download/climate/wc2.1_tiles/", "tile_40_wc2.1_30s_bio.tif"))[[-c(10, 11, 18, 19)]]
bio_30s_tile2 <- rast(paste0("3_raster/download/climate/wc2.1_tiles/", "tile_28_wc2.1_30s_bio.tif"))[[-c(10, 11, 18, 19)]]

#names
names(bio_30s_tile1)
names(bio_30s_tile2)
plot(bio_30s_tile1[[4]])

#union
bio_30s_merge <- merge(bio_30s_tile1, bio_30s_tile2)
plot(bio_30s_merge[[1]])

#mask
bio_30s_mask <- mask(crop(bio_30s_merge, mask_buff), mask_buff)
plot(bio_30s_mask[[1]])

#project
bio_30s_mask <- project(bio_30s_mask, "EPSG:4326")

#update names
names(bio_30s_mask) <- c("bio_01", paste0("bio_", seq(10, 17)), paste0("bio_0", seq(2, 7)))

#folder
sapply("3_raster/var_bio_30s", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#directories
names_30s <- paste0("3_raster/var_bio_30s/", names(bio_30s_mask), ".asc")

#spat to stack
bio_30s_mask <- stack(bio_30s_mask)

#save
wr <- lapply(1:nlayers(bio_30s_mask), function(x) {
  writeRaster(bio_30s_mask[[x]], filename = names_30s[x], format = "ascii", overwrite=T)
})

## 8. Download  ssp245 0.5  variables ##########################################

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#download tiles
bio_30s_ssp245_t1 <- cmip6_tile(-75, -10, 
                            model = "HadGEM3-GC31-LL",
                            ssp = "245",
                            time = "2041-2060",
                            var = "bio",
                            path = "3_raster/download/")[[-c(8, 9, 18, 19)]]
bio_30s_ssp245_t2 <- cmip6_tile(-75, 10, 
                                model = "HadGEM3-GC31-LL",
                                ssp = "245",
                                time = "2041-2060",
                                var = "bio",
                                path = "3_raster/download/")[[-c(8, 9, 18, 19)]]

#names
names(bio_30s_ssp245_t1)
names(bio_30s_ssp245_t2)
plot(bio_30s_ssp245_t2[[4]])

#join
bio_30s_ssp245_merge <- merge(bio_30s_ssp245_t1, bio_30s_ssp245_t2)
plot(bio_30s_ssp245_merge[[1]])

#mask
bio_30s_ssp245_mask <- mask(crop(bio_30s_ssp245_merge, mask_file), 
                            mask_file)
plot(bio_30s_ssp245_mask[[1]])

#update names
names(bio_30s_ssp245_mask) <- c(paste0("bio_0", seq(1, 7)), paste0("bio_", seq(10, 17)))

#folder
sapply("3_raster/var_bio_30s_ssp245", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#directories
names_30s <- paste0("3_raster/var_bio_30s_ssp245/", names(bio_30s_ssp245_mask), 
                    ".asc")
#spat to stack
bio_30s_ssp245_mask <- stack(bio_30s_ssp245_mask)

#save
wr <- lapply(1:nlayers(bio_30s_ssp245_mask), function(x) {
  writeRaster(bio_30s_ssp245_mask[[x]], filename = names_30s[x], format = "ascii", 
              overwrite=T)
})

## 9. Download ssp585 0.5 variables ############################################

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#Download tiles
bio_30s_ssp585_t1 <- cmip6_tile(-75, -10, 
                                model = "HadGEM3-GC31-LL",
                                ssp = "585",
                                time = "2041-2060",
                                var = "bio",
                                path = "3_raster/download/")[[-c(8, 9, 18, 19)]]
bio_30s_ssp585_t2 <- cmip6_tile(-75, 10, 
                                model = "HadGEM3-GC31-LL",
                                ssp = "585",
                                time = "2041-2060",
                                var = "bio",
                                path = "3_raster/download/")[[-c(8, 9, 18, 19)]]

#names
names(bio_30s_ssp585_t1)
names(bio_30s_ssp585_t2)
plot(bio_30s_ssp585_t2[[4]])

#join
bio_30s_ssp585_merge <- merge(bio_30s_ssp585_t1, bio_30s_ssp585_t2)
plot(bio_30s_ssp585_merge[[1]])

#mask
bio_30s_ssp585_mask <- mask(crop(bio_30s_ssp585_merge, mask_file), 
                            mask_file)
plot(bio_30s_ssp585_mask[[1]])

#update names
names(bio_30s_ssp585_mask) <- c(paste0("bio_0", seq(1, 7)), paste0("bio_", seq(10, 17)))

#folder
sapply("3_raster/var_bio_30s_ssp585", function(x)if(!dir.exists(x)) 
  dir.create(x, recursive = T))

#directories
names_30s <- paste0("3_raster/var_bio_30s_ssp585/", names(bio_30s_ssp585_mask), 
                    ".asc")
#spat to stack
bio_30s_ssp585_mask <- stack(bio_30s_ssp585_mask)

#GUARDAR BIOCLIMAS
wr <- lapply(1:nlayers(bio_30s_ssp585_mask), function(x) {
  writeRaster(bio_30s_ssp585_mask[[x]], filename = names_30s[x], format = "ascii", 
              overwrite=T)
})


## 8. Download current 2.5 min  variables   ####################################

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#link
link_2_5_bio_cur <- "https://geodata.ucdavis.edu/climate/worldclim/1_4/grid/cur/bio_2-5m_bil.zip"

#Download
{
  options(timeout = max(1000, getOption("timeout"))) 
  if(!file.exists(file.path("3_raster/download/bio_2-5m_bil", "bio1.bil"))){
    sapply("3_raster/download/bio_2-5m_bil", function(x) if (!dir.exists(x)) dir.create(x, recursive = T))
    download.file(link_2_5_bio_cur, destfile = "3_raster/download/bio_2-5m_bil/bio_2-5m_bil.zip",
      method = "libcurl", mode = "wb", quiet = F)
    unzip("3_raster/download/bio_2-5m_bil/bio_2-5m_bil.zip", exdir = "3_raster/download/bio_2-5m_bil")
  }
  bio_cur_2.5 <- rast(list.files(path = "3_raster/download/bio_2-5m_bil/", pattern = ".bil$", 
    full.names = TRUE))[[-c(10, 11, 18, 19)]]
}

#update names
names(bio_cur_2.5)
names(bio_cur_2.5) <- c("bio_01", paste0("bio_", seq(10, 17)), paste0("bio_0", seq(2, 7)))
plot(bio_cur_2.5[[1]])

#masked
bio_cur_2.5_mask <- raster::mask(crop(bio_cur_2.5, mask_buff), mask_buff)
plot(bio_cur_2.5_mask[[1]])

#project
bio_cur_2.5_mask <- project(bio_cur_2.5_mask, "EPSG:4326")

#directories
suppressWarnings(dir.create("3_raster/var_cur_2.5", recursive = T))
names_2_5_cur <- paste0("3_raster/var_cur_2.5/", names(bio_cur_2.5_mask), ".asc")

#spat to stack
bio_cur_2.5_mask <- stack(bio_cur_2.5_mask)

#save
wr <- lapply(1:nlayers(bio_cur_2.5_mask), function(x) {
  writeRaster(bio_cur_2.5_mask[[x]], filename = names_2_5_cur[x], format = "ascii", 
    overwrite=T)
})

## 9. Download lgm 2.5 variables   #############################################

#load data
mask_file <- st_read(paste0("2_vector/", mask, "_buff_4326.shp"))

#link
link_2_5_bio_lgm <- "https://geodata.ucdavis.edu/climate/cmip5/lgm/cclgmbi_2-5m.zip"

#download
{
  options(timeout = max(1000, getOption("timeout"))) 
  if(!file.exists(file.path("3_raster/download/cclgmbi_2-5m", "cclgmbi1.tif"))){
    sapply("3_raster/download/cclgmbi_2-5m", function(x) if (!dir.exists(x)) dir.create(x, recursive = T))
    download.file(link_2_5_bio_lgm, destfile = "3_raster/download/cclgmbi_2-5m/cc2.1_world.zip",
      method = "libcurl", mode = "wb", quiet = F)
    unzip("3_raster/download/cclgmbi_2-5m/cc2.1_world.zip", exdir = "3_raster/download/cclgmbi_2-5m")
  }
  bio_lgm_2.5 <- rast(list.files(path = "3_raster/download/cclgmbi_2-5m/", pattern = ".tif$", 
    full.names = TRUE))[[-c(10, 11, 18, 19)]]
}

#update names
names(bio_lgm_2.5)
names(bio_lgm_2.5) <- c("bio_01", paste0("bio_", seq(10, 17)), paste0("bio_0", seq(2, 7)))
plot(bio_lgm_2.5[[1]])

#masked
bio_lgm_2.5_mask <- raster::mask(crop(bio_lgm_2.5, mask_buff), mask_buff)
plot(bio_lgm_2.5_mask[[1]])

#project
bio_lgm_2.5_mask <- project(bio_lgm_2.5_mask, "EPSG:4326")

#directories
suppressWarnings(dir.create("3_raster/var_lgm_2.5", recursive = T))
names_2_5_lgm <- paste0("3_raster/var_lgm_2.5/", names(bio_lgm_2.5_mask), ".asc")

bio_lgm_2.5_mask <- stack(bio_lgm_2.5_mask)

#save
wr <- lapply(1:nlayers(bio_lgm_2.5_mask), function(x) {
  writeRaster(bio_lgm_2.5_mask[[x]], filename = names_2_5_lgm[x], format = "ascii", overwrite=T)
})

