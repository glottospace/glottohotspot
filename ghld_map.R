
# Script associated with the following publication:
# Norder and van Gijn. 2026. Global hotspots of linguistic diversity: mapping richness and endemism. 
# Journal: Linguistic Typology

library(glottospace)
library(magrittr)
library(rnaturalearth)
library(tmap)
library(dplyr)
library(sf)
sf_use_s2(FALSE)
tmap_mode(mode = "plot")

intersublab <- getisl(FALSE)
csdf <- data.frame(cs = c(500000, 250000, 100000, 50000, 25000), csl = c("500k", "250k", "100k", "50k", "25k"))

for(i in 1:nrow(csdf)){
csl <- getcsl(i)



# BIVARIATE HOTSPOT MAP
# Combining richness and endemism hotspots (for lang, fam, feat) ---------------------------------------------------------

lingunit <- "lang"
wrldspot_rich <- st_read(dsn = paste0("output/worldspot_richness_", csl, "_", lingunit, "pols", intersublab, ".gpkg"))
wrldspot_end <- st_read(dsn = paste0("output/worldspot_cwe_", csl, "_", lingunit, "pols", intersublab, ".gpkg"))
divcorr(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end)
bivmaplang <- bivcolmap(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end, bivlegend = TRUE)
tmap_save(bivmaplang, filename = paste0("maps/bivariate_", csl, "_", lingunit, intersublab, ".png") )

lingunit <- "fam"
wrldspot_rich <- st_read(dsn = paste0("output/worldspot_richness_", csl, "_", lingunit, "pols.gpkg"))
wrldspot_end <- st_read(dsn = paste0("output/worldspot_cwe_", csl, "_", lingunit, "pols.gpkg"))
divcorr(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end)
bivmapfam <- bivcolmap(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end, bivlegend = TRUE)
tmap_save(bivmapfam, filename = paste0("maps/bivariate_", csl, "_", lingunit, ".png") )

lingunit <- "feature"
wrldspot_rich <- st_read(dsn = paste0("output/worldspot_richness_", csl, "_", lingunit, "pols", intersublab, ".gpkg"))
wrldspot_end <- st_read(dsn = paste0("output/worldspot_cwe_", csl, "_", lingunit, "pols", intersublab, ".gpkg"))
divcorr(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end)
bivmapfeat <- bivcolmap(lingunit = lingunit, worldspot_rich = wrldspot_rich, worldspot_end = wrldspot_end, bivlegend = TRUE)
tmap_save(bivmapfeat, filename = paste0("maps/bivariate_", csl, "_", lingunit, intersublab, ".png") )

# HOTSPOT CONGRUENCE: VENN

# Endemism (lang, fam, feat) and Richness (lang, fam, feat)  --------------

# Richness
div <- "richness"
congrumap_rich <- congrumap(div = div, csl = csl)
tmap_save(congrumap_rich, filename = paste0("maps/congrumap_", csl, "_", div, ".png"))

# Endemism
div <- "cwe"
congrumap_end <- congrumap(div = div, csl = csl)
tmap_save(congrumap_end, filename = paste0("maps/congrumap_", csl, "_", div, ".png"))
}



