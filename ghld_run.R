
# Script associated with the following publication:
# Norder and van Gijn. 2026. Global hotspots of linguistic diversity: mapping richness and endemism. 
# Journal: Linguistic Typology

library(glottospace) # dev version: devtools::install_github("glottospace/glottospace")
# library(geospace)
library(magrittr)
library(rnaturalearth)
library(ggplot2)
library(tmap)
library(sfdep)
library(dplyr)
library(sf)
library(tidyr)
set.seed(42) # for phoible data in glottospace
sf_use_s2(FALSE)
tmap_mode(mode = "plot")
constrain <- c(TRUE, FALSE)[2]

csdf <- data.frame(cs = c(500000, 250000, 100000, 50000, 25000), csl = c("500k", "250k", "100k", "50k", "25k"))

x <- 1
cs <- csdf[x, "cs"] # cell size in km2
csl <- csdf[x, "csl"] # cell size label


# Data preparation: get polygons for linguistic units and generate grids -----------------------------------------------------

# Optionally, take only glottocodes present in both world atlas language polygons (4382) and phoible (2186)

if(constrain == TRUE){
phoible_raw <- glottospace::glottoget("phoible_raw")
phoible_glottocodes <- unique(phoible_raw$id)
glottopols <- st_read("data/worldatlas/cldf/traditional/languages.geojson")
colnames(glottopols)[4] <- "glottocode"
intersub <- intersect(phoible_glottocodes, st_drop_geometry(glottopols$glottocode))
intersublab <- "intersub" # to label output
} else{
  intersublab <- ""
}

# Manually downloaded geojson files (traditional for families) from zenodo: https://zenodo.org/doi/10.5281/zenodo.15287258
# Files are stored in data/worldatlas

# FAMILIES: Obtain polygons for each family
if (!file.exists("data/glottopols/glottopols_fam.gpkg")) {
glottopols <- st_read("data/worldatlas/cldf/traditional/families.geojson")
colnames(glottopols)[4] <- "glottocode"
glottopols <- glottopols %>% group_by(glottocode) %>% summarise() %>% st_cast("MULTIPOLYGON")
glottopols <- sf::st_transform(glottopols, crs = 'ESRI:54012') # World Eckert IV
st_write(obj = glottopols, dsn = paste0("data/glottopols/glottopols_fam.gpkg"), append = F)
# tm_shape(glottopols) + 
#   tm_fill(col = "glottocode", legend.show = F, alpha = 0.6, palette = "Paired") 
}

# LANGUAGES: Obtain polygons for each language 
if (!file.exists(paste0("data/glottopols/glottopols_lang", intersublab, ".gpkg"))) {
glottopols <- st_read("data/worldatlas/cldf/traditional/languages.geojson")
colnames(glottopols)[4] <- "glottocode"
glottopols <- glottopols %>% group_by(glottocode) %>% summarise() %>% st_cast("MULTIPOLYGON")
glottopols <- sf::st_transform(glottopols, crs = 'ESRI:54012') # World Eckert IV
if(constrain == TRUE){
  glottopols <- glottopols[glottopols$glottocode %in% intersub, ]
}
st_write(obj = glottopols, dsn = paste0("data/glottopols/glottopols_lang", intersublab, ".gpkg"), append = F)
# tm_shape(glottopols) + 
#   tm_fill(col = "glottocode", legend.show = F, alpha = 0.6, palette = "Paired") 
}

# FEATURES: Obtain points for each feature and join with polygons at the language level
if (!file.exists(paste0("data/glottopols/glottopols_feature", intersublab, ".gpkg"))) {
  glottopols <- st_read("data/worldatlas/cldf/traditional/languages.geojson")
  colnames(glottopols)[4] <- "id"
  phoible_raw <- glottospace::glottoget("phoible_raw") %>% 
    dplyr::group_by(.data$id)  %>% 
    dplyr::sample_n(1)  %>% 
    dplyr::ungroup()
  if(constrain == TRUE){
    phoible_raw <- phoible_raw[phoible_raw$id %in% intersub, ]
  }
  glottopols <- phoible_param_sf_pols(phoible_data = phoible_raw, pols = glottopols)
  glottopols <- sf::st_transform(glottopols, crs = 'ESRI:54012') # World Eckert IV
  st_write(obj = glottopols, dsn = paste0("data/glottopols/glottopols_feature", intersublab, ".gpkg"), append = F)
  # tm_shape(glottopols) +
  #   tm_fill(col = "ID", legend.show = F, alpha = 0.6)
}

# Generate worldgrids
if (!file.exists("data/worldgrid")) {
  dir.create("data/worldgrid")
  # Create grid
  for(i in 1:nrow(csdf)){
    cs <- csdf[i,"cs"]
    csl <- csdf[i,"csl"]
    grid <- worldgrid(cs = cs, csl = csl, wrldu = wrldu)
  }
}

# # Calculate richness and endemism ------------------------------------------------------
# Families
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
  grid <- sf::st_read(dsn = paste0("data/worldgrid/worldgrid_", csl, ".gpkg") )
  glottopols <- st_read(dsn = paste0("data/glottopols/glottopols_fam.gpkg"))
  wrldjoin <- glottodivgrid(grid = grid, glottopols = glottopols)
  saveRDS(object = wrldjoin, file = paste0("data/worldjoin/worldjoin_", csl, "_fampols") )
}

# Languages
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
  grid <- sf::st_read(dsn = paste0("data/worldgrid/worldgrid_", csl, ".gpkg") )
  glottopols <- st_read(dsn = paste0("data/glottopols/glottopols_lang", intersublab, ".gpkg"))
  wrldjoin <- glottodivgrid(grid = grid, glottopols = glottopols)
  saveRDS(object = wrldjoin, file = paste0("data/worldjoin/worldjoin_", csl, "_langpols", intersublab) )
}

# Features
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
  grid <- sf::st_read(dsn = paste0("data/worldgrid/worldgrid_", csl, ".gpkg") )
  glottopols <- st_read(dsn = paste0("data/glottopols/glottopols_feature", intersublab, ".gpkg"))
  wrldjoin <- glottodivgrid(grid = grid, glottopols = glottopols, id = "ID")
  saveRDS(object = wrldjoin, file = paste0("data/worldjoin/worldjoin_", csl, "_featurepols", intersublab) )
}

# Hotspot analysis --------------------------------------------------------

# Families (endemism and richness)
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
wrldjoin <- readRDS(file = paste0("data/worldjoin/worldjoin_", csl, "_fampols") )
div = "cwe"
wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_fampols.gpkg"), append = F)
div = "richness"
wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_fampols.gpkg"), append = F)
}

# Languages (endemism and richness)
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
  wrldjoin <- readRDS(file = paste0("data/worldjoin/worldjoin_", csl, "_langpols", intersublab) )
  div = "cwe"
  wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
  st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_langpols", intersublab, ".gpkg"), append = F)
  div = "richness"
  wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
  st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_langpols", intersublab, ".gpkg"), append = F)
}

# Features (endemism and richness)
for(i in 1:nrow(csdf)){
  csl <- csdf[i,"csl"]
  wrldjoin <- readRDS(file = paste0("data/worldjoin/worldjoin_", csl, "_featurepols", intersublab) )
  div = "cwe"
  wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
  st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_featurepols", intersublab, ".gpkg"), append = F)
  div = "richness"
  wrldspot <- glottohotspot(wrldjoin = wrldjoin, div = div)
  st_write(obj = wrldspot, dsn = paste0("output/worldspot_", div, "_", csl, "_featurepols", intersublab, ".gpkg"), append = F)
}


