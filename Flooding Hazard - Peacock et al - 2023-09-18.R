# title:  "Census Flooding Hazard Exposure Presentation (Sep-2023)"
# author: "Alexander Abuabara"
# date:   "Sep-18-2023"

library(tidyverse)
library(sf)
library(tigris)
library(tmap); tmap_mode("plot")
library(units)
library(smoothr)
library(janitor)
library(terra)
library(exactextractr)
library(tidycensus)

options(digits = 4, scipen = 9999, na.strings = "NA",
        stringsAsFactors = FALSE, sf_use_s2 = FALSE,
        tigris_use_cache = TRUE, tigris_year = 2020)

`%nin%` = Negate(`%in%`)

setwd("/Users/alexander/Documents/Research/Census/Hazard exposure")
# setwd("/Users/abuabara/Library/Mobile Documents/com~apple~CloudDocs/")

#### Census ####
tx_counties     <- counties("Texas", cb = FALSE) %>% st_transform(3081)
tx_state        <- tx_counties %>% group_by() %>% summarise()
tx_galveston    <- tx_counties[tx_counties$NAME == "Galveston", ]
tx_counties_cb  <- counties("Texas", cb = TRUE) %>% st_transform(3081)
tx_galveston_cb <- tx_counties_cb[tx_counties_cb$NAME == "Galveston", ]

(w1 <- tm_shape(st_buffer(tx_galveston, set_units(40, km)), unit = "mi") +
    tm_borders(lwd = 0) +
    tm_shape(tx_counties) +
    tm_polygons(col = "grey90", lwd = .5) +
    tm_shape(tx_galveston) +
    tm_fill(col = "gold") +
    tm_borders(col = "brown4", lwd = 1) +
    tm_shape(tx_galveston) +
    tm_text("NAME", scale = 0.9, col = "brown4") + 
    tm_layout(main.title = "Case Study",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              frame = FALSE))

(w2 <- tm_shape(st_buffer(tx_galveston, set_units(40, km)), unit = "mi") +
    tm_borders(lwd = 0) +
    tm_shape(tx_counties) +
    tm_polygons(col = "grey90", lwd = .5) +
    tm_shape(tx_galveston) +
    tm_fill(col = "gold") +
    tm_borders(col = "brown4", lwd = 1) +
    tm_shape(tx_galveston_cb) +
    tm_fill(col = "brown4") +
    tm_borders(col = "brown4", lwd = 1) +
    tm_shape(tx_galveston) +
    tm_text("NAME", scale = 1, col = "brown4", xmod = 3.5, ymod = -3.5) + 
    tm_layout(main.title = "Cartographic Boundary",
              main.title.size = 1.1,
              fontfamily = "Georgia",
              frame = FALSE) +
    tm_compass(color.dark = "grey10",
               size = 1,
               position = c(.9,.1)) +
    tm_scale_bar(color.dark = "grey10",
                 text.size = .75,
                 breaks = c(0, 5, 10),
                 position = c(.80,.02)))

tmap_arrange(w1, w2)
tmap_save(tmap_arrange(w1, w2), "w1 w2.png", width = 3000, height = 1400, asp = 0)

#### Water ####
tx_galveston_water <- area_water(county = "Galveston", state = "TX") %>% st_transform(3081)
area_thresh <- set_units(0.01, km2)
tx_galveston_water_simpl <- tx_galveston_water
tx_galveston_water_simpl <- drop_crumbs(tx_galveston_water_simpl, threshold = area_thresh)
tx_galveston_water_simpl <- fill_holes(tx_galveston_water_simpl, threshold = 10 * area_thresh)

tx_galveston_water_simpl <- tx_galveston_water_simpl %>%
  group_by() %>% summarise() %>% transmute(Water = "1") %>%
  st_crop(tx_galveston) %>% st_make_valid()

as.numeric(object.size(tx_galveston_water_simpl)) / as.numeric(object.size(tx_galveston_water))

(w3 <- tm_shape(tx_galveston) +
    tm_polygons(col = "grey") +
    tm_shape(tx_galveston_water) +
    tm_polygons(col = "orange", lwd = .5) +
    tm_shape(tx_galveston) +
    tm_borders(col = "brown4", lwd = 1) +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "brown4", lwd = 1) +
    tm_layout(main.title = "Original Census County Water",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              legend.position = c("right","bottom"),
              frame = FALSE))

(w4 <- tm_shape(tx_galveston) +
    tm_polygons(col = "grey") +
    tm_shape(tx_galveston_water_simpl) +
    tm_fill(col = "blue", alpha = 1) +
    tm_shape(tx_galveston) +
    tm_borders(col = "brown4", lwd = 1) +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "brown4", lwd = 1) +
    tm_layout(main.title = "Simplified County Water",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              legend.position = c("right","bottom"),
              frame = FALSE))

tmap_arrange(w3, w4)
tmap_save(tmap_arrange(w3, w4), "w3 w4.png", width = 3000, height = 1400, asp = 0)

# (tm_shape(tx_galveston) +
#     tm_polygons(col = "gold") +
#     tm_shape(tx_galveston) +
#     tm_borders(col = "red", lwd = 1) +
#     tm_shape(tx_galveston_water_simpl) +
#     tm_fill(col = "lightblue") +
#     tm_shape(tx_galveston) +
#     tm_borders(col = "red", lwd = 1) +
#     tm_layout(main.title = "County Water",
#               main.title.size = .9,
#               fontfamily = "Georgia",
#               legend.position = c("right","bottom"),
#               frame = FALSE))

#### Flooplain ####
unzip(zipfile = "./Data/FEMA NFHL files/48167C_20221120.zip",
      files = c("S_FLD_HAZ_AR.dbf", "S_FLD_HAZ_AR.prj", "S_FLD_HAZ_AR.shp", "S_FLD_HAZ_AR.shx"),
      exdir = "./Data/Tempdir",
      overwrite = TRUE,
      junkpaths = TRUE)

shp0 <- list.files("./Data/Tempdir", pattern = "\\.shp$")

floodplain_ <- st_read(paste("./Data/Tempdir/", shp0[1], sep = "")) %>%
  st_transform(3081) %>%
  mutate(area = set_units(st_area(.), mi2))

if("ZONE" %in% colnames(floodplain_)){names(floodplain_)[names(floodplain_) == "ZONE"] <- "fld_zone"}
if("ZONE_" %in% colnames(floodplain_)){names(floodplain_)[names(floodplain_) == "ZONE_"] <- "fld_zone"}
if("ZONE_SUBTY" %nin% colnames(floodplain_)){floodplain_$zone_subty <- ""}

floodplain_ <- floodplain_ %>%
  filter(!st_is_empty(.)) %>%
  clean_names() %>%
  mutate_if(is.character, ~replace_na(.,"")) %>%
  mutate(fld_zone_and_subty = if_else(zone_subty != "", paste(fld_zone, zone_subty, sep = "-"), paste(fld_zone)))

st_drop_geometry(floodplain_) %>% count(fld_zone_and_subty) %>% arrange(fld_zone_and_subty)

floodplain <- floodplain_ %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>%
  st_make_valid() %>%
  mutate(recode = plyr::revalue(fld_zone_and_subty,
                                c("A" = "100Yr_X",
                                  "AE" = "100Yr_X",
                                  "AE-FLOODWAY" = "100Yr_V",
                                  "AO" = "100Yr_X",
                                  "OPEN WATER" = "100Yr_V",
                                  "VE" = "100Yr_V",
                                  "VE-RIVERINE FLOODWAY SHOWN IN COASTAL ZONE" = "100Yr_V",
                                  "X-0.2 PCT ANNUAL CHANCE FLOOD HAZARD" = "500Yr",
                                  "X-AREA OF MINIMAL FLOOD HAZARD" = "Out",
                                  "X-AREA WITH REDUCED FLOOD RISK DUE TO LEVEE" = "Levee"))) %>%
  filter(recode != "Out") %>% select(fld_zone, zone_subty, fld_zone_and_subty, recode, geometry)

floodplain_diss <- floodplain %>%
  group_by(recode) %>%
  summarise() %>%
  st_intersection(., tx_galveston_cb) %>%
  st_difference(., tx_galveston_water_simpl) %>%
  select(c("GEOID", "recode"))

floodplain_diss <- rbind(floodplain_diss,
                         tx_galveston_water_simpl %>%
                           transmute(
                             GEOID = "48167",
                             recode = "Water",
                             geometry = geometry)) %>%
  mutate(recode = ordered(recode, levels = c("100Yr_V", "100Yr_X", "500Yr", "Levee", "Water")),
         area = set_units(st_area(.), mi2)) %>% select(GEOID, recode, area)

(w5 <- tm_shape(tx_galveston_cb) +
    tm_borders(lwd = .0) +
    tm_shape(floodplain_) +
    tm_fill(col = "fld_zone_and_subty", title = "Flood Zones and Zones Sub-type", alpha = 1) +
    tm_layout(scale = .8, legend.position = c("right","bottom")) +
    tm_shape(tx_counties_cb) +
    tm_borders(col = "black") +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "blue") +
    tm_layout(main.title = "Floodplain in Galveston County (Original)",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              frame = FALSE))

(w6 <- tm_shape(tx_galveston_cb) +
    tm_borders(lwd = .0) +
    tm_shape(floodplain_diss) +
    tm_fill(col = "recode", title = "Floodplain (Recode)",
            alpha = 1,
            label = c("100-year with velocity", "100-year", "500-year", "Levee protected", "Water"),
            palette = c("maroon", "red1", "pink", "yellow", "lightblue")) +
    tm_layout(scale = .8, legend.position = c("right","bottom")) +
    tm_shape(tx_counties_cb) +
    tm_borders(col = "black") +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "blue") +
    tm_layout(main.title = "Floodplain in Galveston County (Reclassified)",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              frame = FALSE))

tmap_arrange(w5, w6)
tmap_save(tmap_arrange(w5, w6), "w5 w6.png", width = 3000, height = 1400, asp = 0)

#### LandScan ####
landscan <- rast("/Users/alexander/Documents/Research/Census/Hazard exposure/Data/landscan-usa-2019-night/landscan-usa-2019-conus-night.tif")
landscan_tx_galveston <- crop(landscan, st_transform(st_buffer(tx_galveston, set_units(10, km)), crs(landscan)))
landscan_tx_galveston <- project(landscan_tx_galveston, "EPSG:3081")
# landscan_tx_galveston = mask(landscan_tx_galveston, tx_galveston)
landscan_tx_galveston[landscan_tx_galveston < 1 ] = NA
# writeRaster(landscan_tx_galveston, "./Data/landscan-usa-2019-night/landscan_tx_galveston.tif", overwrite = TRUE)
# landscan_tx_galveston = rast("/Users/alexander/Documents/Research/Census/Hazard exposure/Data/landscan-usa-2019-night/landscan_tx_galveston.tif")

(w7 <- tm_shape(tx_galveston) +
    tm_borders(lwd = 0) +
    tm_shape(tx_counties_cb) +
    tm_polygons(col = "grey90",
                title = "Galveston County") +
    tm_shape(landscan_tx_galveston) +
    tm_raster("conus_night",
              title = "Estimate",
              saturation = 4,
              palette = c("green2", "yellow2", "blue2"),
              breaks = c(1, 10, 50, Inf)) + 
    tm_layout(scale = .9, legend.position = c("right","bottom")) +
    tm_shape(tx_counties_cb) +
    tm_borders() +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "brown4") +
    tm_layout(main.title = "ORNL LandScan 2019 Population Counts (Conus Night)",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              frame = FALSE) +
    tm_add_legend("line",
                  col = "brown4",
                  size = 2,
                  labels = c("Galveston County (Cartographic Boundary)")))
tmap_save(w7, "w7.png", width = 1800, height = 1400, asp = 0)

tx_galveston_acs_tract_pop <-
  get_acs(geography = "tract",
          variables = c("B02001_001"),
          county = c("Galveston"),
          state = "TX", 
          year = 2021,
          geometry = TRUE) %>%
  filter(!st_is_empty(.)) %>%
  st_transform(3081) %>%
  mutate(area = set_units(st_area(.), mi2),
         density = estimate/area) %>%
  select(-NAME)

paste("estimate county population", round(exact_extract(landscan_tx_galveston, tx_galveston_cb, "sum", progress = TRUE), digits = 0))
# estimate county population 354,112

tx_galveston_acs_tract_pop %>%
  st_drop_geometry() %>% group_by() %>% summarise(total = sum(estimate),
                                                  moe = moe_sum(estimate, moe))
#   total    moe
# 347,084 38,817

tx_galveston_acs_tract_pop <- tx_galveston_acs_tract_pop %>%
  mutate(landscan_population = exact_extract(landscan_tx_galveston, ., "sum", progress = TRUE),
         diff = estimate - landscan_population) %>%
  select(GEOID, estimate, moe, landscan_population, diff)

mean(tx_galveston_acs_tract_pop$estimate)
sd(tx_galveston_acs_tract_pop$estimate)
# 3403
# 1796

mean(tx_galveston_acs_tract_pop$landscan_population)
sd(tx_galveston_acs_tract_pop$landscan_population)
# 3472
# 1793

par(family = "Georgia")
hist(tx_galveston_acs_tract_pop$diff,
     col = "skyblue",
     main = "Histogram for Census Estimate - Landscan Counts")

t.test(tx_galveston_acs_tract_pop$estimate,
       tx_galveston_acs_tract_pop$landscan_population,
       alternative = c("two.sided"), paired = TRUE, var.equal = FALSE, conf.level = 0.95)
# Paired t-test
# data:  tx_galveston_acs_tract_pop$estimate and tx_galveston_acs_tract_pop$landscan_population
# t = -1, df = 101, p-value = 0.3
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval: -203.28 65.55
# sample estimates: mean difference -68.87 

# Population at risk
(exact_extract(landscan_tx_galveston, filter(floodplain_diss, recode == "Water"), "sum", progress = TRUE))
# 5,141

(exact_extract(landscan_tx_galveston, filter(floodplain_diss, recode == "100Yr_V"), "sum", progress = TRUE))
# 7,712

(exact_extract(landscan_tx_galveston, filter(floodplain_diss, recode == "100Yr_X"), "sum", progress = TRUE))
# 99,136

(exact_extract(landscan_tx_galveston, filter(floodplain_diss, recode == "500Yr"), "sum", progress = TRUE))
# 81,637

(exact_extract(landscan_tx_galveston, filter(floodplain_diss, recode == "Levee"), "sum", progress = TRUE))
# 40,243

354112 - 40243 - 81637 - 99136 - 7712 - 5141
# 120,243

st_drop_geometry(floodplain_diss) %>% adorn_totals()
# GEOID  recode         area
# 48167 100Yr_V 102.14 [mi2]
# 48167 100Yr_X  93.84 [mi2]
# 48167   500Yr  61.07 [mi2]
# 48167   Levee  24.81 [mi2]
# 48167   Water 492.70 [mi2]
# Total       - 774.55 [mi2]

#### Parcels ####
# Source data: https://galvestoncad.org/gis-data/
parcels_ <- st_read("/Users/alexander/Downloads/parcels/parcels.shp") %>%
  filter(!st_is_empty(.)) %>%
  filter(st_geometry_type(.) %in% c("POLYGON", "MULTIPOLYGON")) %>%
  st_make_valid()

parcels <- parcels_ %>%
  filter(grepl('R', LANDUSE)) %>%
  filter(VAL23IMP > 0) %>%
  st_make_valid() %>%
  st_transform(3081)

downtown_bb <- st_bbox(st_buffer(tx_galveston_acs_tract_pop %>%
                                   filter(GEOID %in% c("48167725400")), set_units(5, km)))

(w8 <- tm_shape(tx_galveston_cb, bbox = downtown_bb) +
    tm_polygons(col = "grey90", lwd = 0) +
    tm_shape(floodplain_diss %>%
               filter(recode != "Levee") %>%
               mutate(recode = factor(recode,
                                      levels = c("100Yr_V", "100Yr_X", "500Yr", "Water")))) +
    tm_fill(col = "recode", title = "Floodplain (Recode)",
            alpha = 1,
            label = c("100-year with velocity", "100-year", "500-year", "Water"),
            palette = c("maroon", "red1", "pink", "lightblue")) +
    tm_layout(scale = .9, legend.position = c("right","bottom")) +
    tm_shape(tx_counties_cb) +
    tm_borders(col = "black") +
    tm_shape(tx_galveston_cb) +
    tm_borders(col = "blue") +
    tm_shape(parcels_, bbox = downtown_bb) +
    tm_fill(col = "grey40", alpha = .3) +
    tm_shape(parcels, bbox = downtown_bb) +
    tm_fill(col = "green2", alpha = .6) +
    tm_add_legend("fill",
                  col = "grey40",
                  border.col = NA,
                  alpha = .6,
                  size = 2,
                  labels = c("Parcels")) +
    tm_add_legend("fill",
                  col = "green2",
                  border.col = NA,
                  alpha = .6,
                  size = 2,
                  labels = c("Resid. parcels with improvement")) +
    tm_layout(main.title = "View of Galveston Island Downtown Area",
              main.title.size = 1.2,
              fontfamily = "Georgia",
              legend.bg.color = "white",
              frame = FALSE))
tmap_save(w8, "w8.png", width = 1800, height = 1400, asp = 0)

(tot <- length(st_intersects(tx_galveston, parcels)[[1]]))
(c1 <- length(st_intersects(filter(floodplain_diss, recode == "Water"), st_centroid(parcels))[[1]]))
(c2 <- length(st_intersects(filter(floodplain_diss, recode == "100Yr_V"), st_centroid(parcels))[[1]]))
(c3 <- length(st_intersects(filter(floodplain_diss, recode == "100Yr_X"), st_centroid(parcels))[[1]]))
(c4 <- length(st_intersects(filter(floodplain_diss, recode == "500Yr"), st_centroid(parcels))[[1]]))
(c5 <- length(st_intersects(filter(floodplain_diss, recode == "Levee"), st_centroid(parcels))[[1]]))
c1 + c2 + c3 + c4 + c5
tot - (c1 + c2 + c3 + c4 + c5)

# (a1 <- length(st_intersects(filter(floodplain_diss, recode == "Water"), parcels)[[1]]))
# (a2 <- length(st_intersects(filter(floodplain_diss, recode == "100Yr_V"), parcels)[[1]]))
# (a3 <- length(st_intersects(filter(floodplain_diss, recode == "100Yr_X"), parcels)[[1]]))
# (a4 <- length(st_intersects(filter(floodplain_diss, recode == "500Yr"), parcels)[[1]]))
# (a5 <- length(st_intersects(filter(floodplain_diss, recode == "Levee"), parcels)[[1]]))
# tot - (a1 + a2 + a3 + a4 + a5)
# a1 + a2 + a3 + a4 + a5
#### END ####