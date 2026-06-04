
# Script associated with the following publication:
# Norder and van Gijn. 2026. Global hotspots of linguistic diversity: mapping richness and endemism. 
# Journal: Linguistic Typology

library(glottospace)
# library(geospace)
library(magrittr)
library(rnaturalearth)
library(tmap)
library(dplyr)
library(sf)
library(lwgeom)
sf_use_s2(FALSE)

# Functions ---------------------------------------------------------------



basemap <- function(union = TRUE){
  
  if(!file.exists("data/ne_50m_admin_0_countries.shp")){
    # wrld <- rnaturalearth::ne_download(scale = 50, type = "countries",
    #                                     category = "cultural", returnclass = "sf")
    utils::download.file(file.path("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/50m/cultural/ne_50m_admin_0_countries.zip"), zip_file <- tempfile())
    utils::unzip(zip_file, exdir = "data/")
  }
  wrld <- sf::st_read(dsn = "data/ne_50m_admin_0_countries.shp")
  wrld <- wrld %>% dplyr::filter(GEOUNIT != "Antarctica")
  wrld <- sf::st_transform(wrld, crs = 'ESRI:54012') # World Eckert IV
  if(union == TRUE){
    wrldu <- wrld %>% sf::st_geometry() %>% sf::st_union()
    # wrldu <- sf::st_transform(wrldu, crs = 'ESRI:54012') # World Eckert IV
    
    # Remove small islands
    at <- units::set_units(20000, km^2)
    wrldu <- smoothr::drop_crumbs(x = wrldu, threshold = at, drop_empty = TRUE)
    wrldu <- sfheaders::sf_remove_holes(wrldu) # remove holes
    out <- wrldu
  } else {
    out <- wrld
  }
  return(out)
}

getcsl <- function(x){
  # get cellsize label
  csdf <- data.frame(cs = c(500000, 250000, 100000, 50000, 25000), csl = c("500k", "250k", "100k", "50k", "25k"))
  # cs <- csdf[x, "cs"] # cell size in km2
  csl <- csdf[x, "csl"] # cell size label
  csl
}

getisl <- function(constrain = FALSE){
  # Get intersub label to constrain analysis to glottocodes present in both world atlas language polygons (4382) and phoible (2186)
if(constrain == TRUE){
  intersublab <- "intersub" # to label output
} else{
  intersublab <- ""
}
  intersublab
}
  
worldgrid <- function(cs = NULL, csl = NULL, wrldu = NULL){
  # create hexagon grid specifying cell area in km2 # https://github.com/r-spatial/sf/issues/1505
  cellarea <- cs * (1e+6)
  cellsize <- 2 * sqrt(cellarea/((3*sqrt(3)/2))) * sqrt(3)/2
  
  # wrldux <- st_buffer(wrldu, dist = (1000 * (1e+6))  )
  
  grid <- wrldu %>%
    sf::st_make_grid(cellsize = cellsize, square = FALSE) %>%
    st_sf()
  # hex_area <- units::set_units(st_area(grid), km^2) # Check area
  
  # Keep only overlapping (st_intersects)
  grid <- grid[st_intersects(x = grid, wrldu, sparse = FALSE),]
  grid$id <- as.character(1:nrow(grid))
  sf::st_write(obj = grid, dsn = paste0("data/worldgrid/worldgrid_", csl, ".gpkg"), append = F)
  return(grid)
}

glottodivgrid <- function(grid = NULL, glottopols = NULL, id = "glottocode"){
  # Determine which linguistic units fall within a grid cell (intersection matrix)
  
  intdf <- sf::st_intersects(grid, glottopols, sparse = F) %>%
    as.data.frame()
  colnames(intdf) <- glottopols[[id]]
  
  # Endemism by family
  intmat <- as.matrix(intdf)
  aoo <- colSums(intmat) # Area of occupancy (in case of point data: extent of occupancy for languages wouldn't make sense)
  ir <- 1/aoo # inverse range
  
  for(i in 1:nrow(intmat)){
    intmat[i, ] <- ifelse(intmat[i,], ir, 0)
  }
  endemism <- rowSums(intmat)
  intdf$we <- endemism

  # Richness
  richness <- rowSums(intdf)
  intdf$richness <- richness
  intdf <- tibble::rownames_to_column(intdf, "id")
  
  # A high CWE indicates that a large proportion of the languages/families/features in that cell are spatially restricted.
  
  intdf$cwe <- intdf$we/intdf$richness
  intdf$cwe <- ifelse(is.na(intdf$cwe), 0, intdf$cwe) # assign 0 to cells without languages (zero values division)
  
  wrldjoin <- dplyr::left_join(x = intdf, y = grid, by = "id") %>% sf::st_sf()
  
  # Add rank columns
  wrldjoin <- wrldjoin %>% mutate(cwe_rank = rank(desc(cwe), ties.method = 'first')) # alternative could be dense_rank to rank ties equally but this is not what we want, because we're interested in percentages
  wrldjoin <- wrldjoin %>% mutate(cwe_perc = cwe_rank/nrow(wrldjoin)*100)
  
  wrldjoin <- wrldjoin %>% mutate(we_rank = rank(desc(we), ties.method = 'first'))
  wrldjoin <- wrldjoin %>% mutate(we_perc = we_rank/nrow(wrldjoin)*100)
  
  wrldjoin <- wrldjoin %>% mutate(richness_rank = rank(desc(richness), ties.method = 'first')) 
  wrldjoin <- wrldjoin %>% mutate(richness_perc = richness_rank/nrow(wrldjoin)*100)
  
  return(wrldjoin)
}

glottohotspot <- function(wrldjoin, div){
  colnames(wrldjoin)[colnames(wrldjoin) == div] <- "div"
  wrldjoin_nb <- wrldjoin %>% select(id, div)
  
  wrldjoin_nb <- wrldjoin_nb %>% 
    mutate(nb = sfdep::st_contiguity(.),
           wt = st_weights(nb),
           div_lag = st_lag(div, nb, wt)
    )
  
  # Calculate Gi*
  wrldjoin_gi <- wrldjoin_nb %>%  
    mutate(
      Gi = local_gstar_perm(div, nb, wt, nsim = 9999, alternative = "greater")
    ) %>% 
    unnest(Gi)
  
  wrldjoin_hotspots <- wrldjoin_gi %>% 
    select(gi_star, p_value) %>% 
    mutate(
      siglevel = case_when(
        gi_star > 0 & p_value <= 0.001 ~ "sig_veryhigh",
        gi_star > 0 & p_value <= 0.005 ~ "sig_high",
        gi_star > 0 & p_value <= 0.05 ~ "sig_low",
        TRUE ~ "notsig_or_gineg"
      ),
      siglevel = factor(
        siglevel,
        levels = c("sig_veryhigh", "sig_high", "sig_low",
                   "notsig_or_gineg")
      )
    )
 
  # Test for global clustering
  print(global_g_test(wrldjoin_nb$div, wrldjoin_nb$nb, wrldjoin_nb$wt))
  return(wrldjoin_hotspots)
}

# phoible parameter sf polygon object (modified from glottospace package)
phoible_param_sf_pols <- function(phoible_data, pols){
  param_idx <- colnames(phoible_data) |>
    sapply(
      FUN = function(x){
        nchar(x) == 32
      }
    )
  param_ids <- colnames(phoible_data)[param_idx]
  
  phoible_join <- left_join(x = pols, y = phoible_data, by = "id")
  
  data <- phoible_join[, param_ids] %>% stats::na.omit()
  
  # Character vector of parameter IDs
  paramidvec <- st_drop_geometry(data)  %>% 
    apply(MARGIN = 2,
          FUN = function(x){
            !all(x == "absent")
          })  %>% 
    unlist() %>% 
    which() %>% 
    names()
  
  data <- data[, paramidvec]
  
  # data <- data[1:100, 1:100]
  
  result_list <- lapply(paramidvec, function(id) {
    subset_sf <- data[which(data[[id]] != "absent"), ]
    combined_geom <- st_combine(subset_sf)
    st_sf(
      data.frame(ID = id),
      geometry = st_sfc(combined_geom, crs = st_crs(data))
    )
  })
  
  do.call(rbind, result_list)
 
}

bivcolmap <- function(lingunit, worldspot_rich, worldspot_end, bivlegend = TRUE){
  countries <- ne_countries() %>% filter(geounit != "Antarctica") %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  ocean <- sf::st_read(dsn = paste0("data/worldbase/ne_10m_ocean/ne_10m_ocean.shp") ) %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  bivcolmat <- matrix(c("white", "#6c7b99", "#3d5479", "#00305a",
                        "#ffe585", "#c97765", "#c97765", "#c97765",
                        "#ffdc58", "#c97765", "#b04937", "#b04937",
                        "#ffd412", "#c97765", "#b04937", "#94090d"),
                      byrow = TRUE, ncol = 4)
  
  if(lingunit == "lang"){
    maptitle <- "Language diversity hotspots"
  } else if(lingunit == "fam"){
    maptitle <- "Genealogical diversity hotspots"
  } else if(lingunit == "feature"){
    maptitle <- "Structural diversity hotpots"
  }
  
  
  richsig <- wrldspot_rich %>% rename(richsig = siglevel)
  richsig$richsig <- factor(richsig$richsig, levels = c("notsig_or_gineg", "sig_low", "sig_high", "sig_veryhigh")) # Order legend for full bivariate map, in case all categories are included

  endsig <- wrldspot_end %>% rename(endsig = siglevel)
  endsig$endsig <- factor(endsig$endsig, levels = c("notsig_or_gineg", "sig_low", "sig_high", "sig_veryhigh")) # Order legend for full bivariate map, in case all categories are included
  
  # Bivariate map
  wrldspot_comb <- cbind(richsig, endsig) %>%
    select(c("richsig", "endsig")) %>%
    filter(richsig != "notsig_or_gineg" | endsig != "notsig_or_gineg") # Only show cell for which at least richness or endemsim is a hotspot
  
  if(bivlegend == FALSE){
    bivmap <- tm_shape(countries) +
      tm_polygons(fill = "white") +
      tm_shape(wrldspot_comb) +
      tm_fill(fill = tm_vars(c("richsig", "endsig"), multivariate = TRUE),
                  fill.scale = tm_scale_bivariate(
                    scale1 = tm_scale_categorical(labels = c("", "`", "*", "**")),
                    scale2 = tm_scale_categorical(labels = c("", "`", "*", "**")),
                    values = bivcolmat),
                  fill_alpha = 0.7,
                  fill.legend = tm_legend_hide()
      ) +
      tm_shape(ocean) +
      tm_polygons() +
      tm_graticules() +
      tm_title_out(text = maptitle) 
  } else {
    bivmap <- tm_shape(countries) +
      tm_polygons(fill = "white") +
      tm_shape(wrldspot_comb) +
      tm_fill(fill = tm_vars(c("richsig", "endsig"), multivariate = TRUE),
                  fill.scale = tm_scale_bivariate(
                    scale1 = tm_scale_categorical(labels = c("", "`", "*", "**")),
                    scale2 = tm_scale_categorical(labels = c("", "`", "*", "**")),
                    values = bivcolmat),
                  fill_alpha = 0.7,
                  fill.legend = tm_legend_bivariate(
                    xlab = "endemism",
                    ylab = "richness",
                  )
      ) +
      tm_shape(ocean) +
      tm_polygons() +
      tm_graticules()  +
      tm_title_out(text = maptitle) + 
      tm_layout(
        legend.position = c("left","bottom")
      )
  }
  bivmap
}

congruvennlegend <- function(){
  
  world_eck4 <- st_transform(World, "+proj=eck4")
  bb <- st_bbox(world_eck4)
  
  
  y_shift <- (bb["ymax"] - bb["ymin"]) * 0.05
  
  # Choose radius relative to map size
  r <- (bb["xmax"] - bb["xmin"]) * 0.06
  
  circles_sf <- st_sf(
    region = c("lang", "fam", "feat"),
    geometry = st_sfc(
      st_buffer(
        st_point(c(bb["xmin"] + r*1.5, bb["ymin"] + r*1.5 + y_shift)), r
      ),
      st_buffer(
        st_point(c(bb["xmin"] + r*3.0, bb["ymin"] + r*1.5 + y_shift)), r
      ),
      st_buffer(
        st_point(c(bb["xmin"] + r*2.25, bb["ymin"] + r*3.0 + y_shift)), r
      ),
      crs = st_crs(world_eck4)
    )
  )
  
  lang <- circles_sf[1, ]
  fam <- circles_sf[2, ]
  feat <- circles_sf[3, ]
  
  # Triple intersection
  langfamfeat <- st_intersection(lang, fam, feat)
  
  # Exclusive pairwise intersections
  langfam_only <- st_difference(st_intersection(lang, fam), feat)
  langfeat_only <- st_difference(st_intersection(lang, feat), fam)
  famfeat_only <- st_difference(st_intersection(fam, feat), lang)
  
  # Exclusive single regions
  lang_only <- st_difference(lang, st_union(fam, feat))
  fam_only <- st_difference(fam, st_union(lang, feat))
  feat_only <- st_difference(feat, st_union(lang, fam))
  
  # Triple intersection 
  langfamfeat <- st_difference(langfamfeat, langfam_only)
  
  st_geometry_type(langfam_only)
  st_geometry_type(langfamfeat)
  
  
  make_region <- function(x, name) {
    x$region <- name
    st_cast(x[, "region"], "MULTIPOLYGON")
  }
  
  venn_sf <- rbind(
    make_region(lang_only,  "lang"),
    make_region(fam_only,  "fam"),
    make_region(feat_only,  "feat"),
    make_region(langfam_only, "langfam"),
    make_region(langfeat_only, "langfeat"),
    make_region(famfeat_only, "famfeat"),
    make_region(langfamfeat,     "langfamfeat")
  )
  
  labels_sf <- st_centroid(
    circles_sf[, "region"]
  )
  labels_sf[, 1] <- c("Language", "Genealogical", "Structural")
  
  congruvennlegend <- tm_shape(venn_sf) +
    tm_polygons(
      col = "region",
      palette = c(
        "lang" = "blue",
        "fam" = "red",
        "feat" = "yellow",
        "langfam" = "purple",
        "langfeat" = "green",
        "famfeat" = "orange",
        "langfamfeat" = "black"
      ),
      alpha = 0.7,
      border.col = "black"
    ) + tm_shape(labels_sf) +
    tm_text("region", size = 0.45, fontface = "bold") +
    tm_layout(legend.show = FALSE)
  
  congruvennlegend
}

congrumap <- function(div = NULL, csl = NULL){
  countries <- ne_countries() %>% filter(geounit != "Antarctica") %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  ocean <- sf::st_read(dsn = paste0("data/worldbase/ne_10m_ocean/ne_10m_ocean.shp") ) %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  
  if(div == "richness"){
    maptitle <- "Richness hotspot congruence"
  } else if(div == "cwe"){
    maptitle <- "Endemism hotspot congruence"
  } 
  
  wrldspot_ls <- list()
  for(lingunit in c("lang", "fam", "feature")){
    wrldspot_ls[[lingunit]] <- st_read(dsn = paste0("output/worldspot_", div, "_", csl, "_", lingunit, "pols.gpkg"))
  }
  
  langsig <- wrldspot_ls[["lang"]]$siglevel
  famsig <- wrldspot_ls[["fam"]]$siglevel
  featsig <- wrldspot_ls[["feature"]]$siglevel
  
  wrldspot_comb <- cbind(wrldspot_ls[["lang"]], famsig, featsig) %>% 
    select(c("siglevel", "famsig", "featsig")) 
  colnames(wrldspot_comb)[1] <- "langsig"
  
  # Should all be false!
  any(unique(wrldspot_comb$langsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  any(unique(wrldspot_comb$famsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  any(unique(wrldspot_comb$featsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  
  wrldspot_comb <- wrldspot_comb %>% mutate(overlap = case_when(
    langsig != "notsig_or_gineg" & famsig != "notsig_or_gineg"  & featsig != "notsig_or_gineg"  ~ "langfamfeat",
    langsig != "notsig_or_gineg" & famsig != "notsig_or_gineg"  & featsig == "notsig_or_gineg"  ~ "langfam",
    langsig != "notsig_or_gineg" & featsig != "notsig_or_gineg" & famsig == "notsig_or_gineg" ~ "langfeat",
    famsig != "notsig_or_gineg" & featsig != "notsig_or_gineg" & langsig == "notsig_or_gineg" ~ "famfeat",
    langsig != "notsig_or_gineg" & famsig == "notsig_or_gineg"  & featsig == "notsig_or_gineg"  ~ "lang",
    langsig == "notsig_or_gineg" & featsig != "notsig_or_gineg" & famsig == "notsig_or_gineg" ~ "feat",
    famsig != "notsig_or_gineg" & featsig == "notsig_or_gineg" & langsig == "notsig_or_gineg" ~ "fam",
    TRUE ~ "notsig_or_gineg"),
    overlap = factor(
      overlap,
      levels = c("lang", "fam", "feat",
                 "langfam", "langfeat", "famfeat",
                 "langfamfeat",
                 "notsig_or_gineg")
    )
  )
  wrldspot_comb <- wrldspot_comb %>% filter(overlap != "notsig_or_gineg")
  wrldspot_comb$overlap <- droplevels(wrldspot_comb$overlap)
  
  vennlegend <- congruvennlegend()
  
  congrumap <- tm_shape(countries) +
    tm_polygons(fill = "white") +
    tm_shape(wrldspot_comb) +
    tm_fill(fill = "overlap",
            fill.scale = tm_scale_categorical(values = c("blue", "red", "yellow","purple", "green", "orange", "black") ), 
            fill_alpha = 0.7
    ) +
    tm_shape(ocean) +
    tm_polygons() +
    tm_graticules() +
    tm_title_out(text = maptitle) +
    tm_layout(legend.show = FALSE) +
    vennlegend
  congrumap
}

divcorr <- function(lingunit, worldspot_rich, worldspot_end){

  if(lingunit == "lang"){
    maptitle <- "Language diversity hotspots"
  } else if(lingunit == "fam"){
    maptitle <- "Genealogical diversity hotspots"
  } else if(lingunit == "feature"){
    maptitle <- "Structural diversity hotpots"
  }
  
  
  richsig <- wrldspot_rich %>% rename(richsig = siglevel)
  endsig <- wrldspot_end %>% rename(endsig = siglevel)

  
  wrldspot_comb <- cbind(richsig, endsig) %>%
    select(c("richsig", "endsig")) %>%
    mutate(
      across(
        c(richsig, endsig),
        ~case_when(
          .x %in% c("sig_low", "sig_high", "sig_veryhigh") ~ 1,
          .x %in% "notsig_or_gineg" ~ 0
        )
      )
    )
  
  # Create a contingency table
  contingency_table <- table(wrldspot_comb$richsig, wrldspot_comb$endsig)

  # Chi-square test
  reschi <- chisq.test(contingency_table)

  # Pearson correlation
  rescor <- cor(as.numeric(wrldspot_comb$richsig), as.numeric(wrldspot_comb$endsig), method = "pearson")
  
  paste0("Pearson correlation ", round(rescor, 2), " p-value: ", reschi$p.value)
}

congru_explore <- function(div = NULL, csl = NULL){
  countries <- ne_countries() %>% filter(geounit != "Antarctica") %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  ocean <- sf::st_read(dsn = paste0("data/worldbase/ne_10m_ocean/ne_10m_ocean.shp") ) %>% sf::st_transform(wrld, crs = 'ESRI:54012')
  
  if(div == "richness"){
    maptitle <- "Richness hotspot congruence"
  } else if(div == "cwe"){
    maptitle <- "Endemism hotspot congruence"
  } 
  
  wrldspot_ls <- list()
  for(lingunit in c("lang", "fam", "feature")){
    wrldspot_ls[[lingunit]] <- st_read(dsn = paste0("output/worldspot_", div, "_", csl, "_", lingunit, "pols.gpkg"))
  }
  
  langsig <- wrldspot_ls[["lang"]]$siglevel
  famsig <- wrldspot_ls[["fam"]]$siglevel
  featsig <- wrldspot_ls[["feature"]]$siglevel
  
  wrldspot_comb <- cbind(wrldspot_ls[["lang"]], famsig, featsig) %>% 
    select(c("siglevel", "famsig", "featsig")) 
  colnames(wrldspot_comb)[1] <- "langsig"
  
  # Should all be false!
  any(unique(wrldspot_comb$langsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  any(unique(wrldspot_comb$famsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  any(unique(wrldspot_comb$featsig) %nin% c("notsig_or_gineg", "sig_high", "sig_low", "sig_veryhigh"))
  
  wrldspot_comb <- wrldspot_comb %>% mutate(overlap = case_when(
    langsig != "notsig_or_gineg" & famsig != "notsig_or_gineg"  & featsig != "notsig_or_gineg"  ~ "langfamfeat",
    langsig != "notsig_or_gineg" & famsig != "notsig_or_gineg"  & featsig == "notsig_or_gineg"  ~ "langfam",
    langsig != "notsig_or_gineg" & featsig != "notsig_or_gineg" & famsig == "notsig_or_gineg" ~ "langfeat",
    famsig != "notsig_or_gineg" & featsig != "notsig_or_gineg" & langsig == "notsig_or_gineg" ~ "famfeat",
    langsig != "notsig_or_gineg" & famsig == "notsig_or_gineg"  & featsig == "notsig_or_gineg"  ~ "lang",
    langsig == "notsig_or_gineg" & featsig != "notsig_or_gineg" & famsig == "notsig_or_gineg" ~ "feat",
    famsig != "notsig_or_gineg" & featsig == "notsig_or_gineg" & langsig == "notsig_or_gineg" ~ "fam",
    TRUE ~ "notsig_or_gineg"),
    overlap = factor(
      overlap,
      levels = c("lang", "fam", "feat",
                 "langfam", "langfeat", "famfeat",
                 "langfamfeat",
                 "notsig_or_gineg")
    )
  )
  wrldspot_comb <- wrldspot_comb %>% filter(overlap != "notsig_or_gineg")
  wrldspot_comb$overlap <- droplevels(wrldspot_comb$overlap)
  
  vennlegend <- congruvennlegend()
  
  congrumap <- tm_shape(wrldspot_comb) +
    tm_fill(fill = "overlap",
            fill.scale = tm_scale_categorical(values = c("blue", "red", "yellow","purple", "green", "orange", "black") ), 
            fill_alpha = 0.7
    ) +
    tm_shape(ocean) +
    tm_polygons() +
    tm_graticules() +
    tm_title_out(text = maptitle) +
    tm_layout(legend.show = FALSE) 
  congrumap
}
