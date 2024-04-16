# explore new lake data from KM, make spatial, and join relevant lakes 
# with catchments
# Jenny Hansen
# 06 December 2023

# working in FreshRestore project

# Load required libraries -------------------------------------------------

library(readxl)
library(dplyr)
library(janitor)
library(sf)
library(purrr)
library(tidyr)
library(stringr)

# Import data -------------------------------------------------------------

# new data & fish data from Kim Magnus
new_data <- read_excel("data/new/Oversikt-lokaliteter-lagt inn i basen.xls") %>% 
  clean_names()

fish_data <- read.delim("data/new/innlandsfisk.txt", header = TRUE, sep = ";") %>% 
  clean_names()

# innsjo & regine shapefiles came from GeoNorge.no
lakes <- st_read("vector/Innsjo_Innsjo.shp") %>% 
  st_transform(st_crs(25833)) %>% 
  st_zm() %>% 
  st_make_valid()

catch <- st_read("vector/Nedborfelt_RegineEnhet.shp") %>% 
  st_transform(st_crs(25833)) %>% 
  st_zm() %>% 
  st_make_valid

# Clean up new data -------------------------------------------------------

# function to obtain a unique row for each year
# this is v. complicated because the year column has so many 
# weird formats

parse_years <- function(year_str) {
  # split string by commas
  parts <- str_split(year_str, ",\\s*")[[1]]
  years <- c()
  
  # function to determine the correct century for a year
  determine_century <- function(year, last_full_year) {
    if (year < 100) {
      # if year is two digits, determine the correct century
      century <- (last_full_year %/% 100) * 100
      if (year + century > current_year) {
        # if the result is in future, subtract 100 years
        century <- century - 100
      }
      return(year + century)
    } else {
      return(year)
    }
  }
  
  # get the current year
  current_year <- as.integer(format(Sys.Date(), "%Y"))
  
  # loop through parts and handle ranges
  for (part in parts) {
    if (str_detect(part, "-")) {
      # handle range
      range <- str_split(part, "-")[[1]]
      start_year <- as.integer(range[1])
      end_year <- as.integer(range[2])
      
      # adjust for short year formats (e.g., 2000-01)
      if (nchar(range[2]) < 4) {
        end_year <- determine_century(end_year, start_year)
      }
      
      years <- c(years, seq(start_year, end_year))
    } else {
      # handle single year
      year <- as.integer(part)
      year <- determine_century(year, 
                                if (length(years) > 0) tail(years, 1) else current_year)
      years <- c(years, year)
    }
  }
  
  return(years)
}

# apply function to new data
tidy_data <- new_data %>%
  mutate(year = map(ar_lagt_inn, parse_years)) %>%
  unnest(year)

# QAQC for year separation (ran through 10 samplings; looks good)
sampled_tidy_data <- tidy_data %>% 
  select(ar_lagt_inn, year) %>% 
  sample_n(20)
sampled_tidy_data
range(tidy_data$year) # 1975-2013
table(tidy_data$year)
table(tidy_data$loknavn)

# Find lakes common to new_data & fish_data -------------------------------

fish_data %>% 
  filter(fk_innsjo_nr %in% tidy_data$innsjo_nr) %>% 
  group_by(fk_innsjo_nr) %>% 
  slice(1) # 334

tidy_data %>% 
  filter(innsjo_nr %in% fish_data$fk_innsjo_nr) %>% 
  group_by(innsjo_nr) # 334

# filter for Swedish lakes 
tidy_data %>% 
  filter(str_detect(innsjo_nr, "^Sverige")) %>% 
  pull(loknavn) # n = 4

# drop the lakes beginning with 'Sverige_'
# I previously tried to match names with the Swedish
# lake GIS, but could not find matches
shared_lakes <- tidy_data %>%
  filter(innsjo_nr %in% fish_data$fk_innsjo_nr) %>% 
  mutate(innsjo_nr = as.numeric(innsjo_nr)) %>% 
  filter(!is.na(innsjo_nr))

# Make spatial (join with lakes) ------------------------------------------

lakes_sf <- lakes %>% 
  left_join(shared_lakes, by = c("vatnLnr" = "innsjo_nr"), keep = T) %>% 
  filter(!is.na(innsjo_nr)) %>% 
  select(-c(hoyde_moh:ekspType),
         -objType)
names(lakes_sf)
st_crs(lakes_sf) # 25833
mapview::mapview(lakes_sf)

# Join with catchments ----------------------------------------------------

filtered_catch <- catch %>% 
  st_filter(lakes_sf, .predicate = st_intersects)

joined_catch <- st_join(filtered_catch, lakes_sf[, c("innsjo_nr", "year")], 
                        join = st_intersects) %>% 
  select(vassdragNr, innsjo_nr, year)

mapview::mapview(joined_catch)

# Remove lake portion from catchment --------------------------------------

# needs to be done rowwise to avoid multiple 'holes' in shared catchments

# set 'geometry' as the active geometry column
st_geometry(lakes_sf) <- "geometry"
st_geometry(joined_catch) <- "geometry"

# remove lake area from catchments
catch_minus_lake <- joined_catch %>% 
  rowwise()  %>%  
  mutate(geometry = 
           st_difference(geometry, 
                         lakes_sf[lakes_sf$innsjo_nr == innsjo_nr, ]$geometry[1]))  %>%  
  ungroup()  %>%  
  st_as_sf()

mapview::mapview(catch_minus_lake) +
  mapview::mapview(lakes_sf, col.regions = "orange")

# Check/remove duplicates -------------------------------------------------

catch_minus_lake %>% 
  unite("check_1", c("vassdragNr", "innsjo_nr"), remove = F) %>% 
  unite("check_column", c("check_1", "year")) %>% 
  filter(duplicated(check_column))
# no duplicates

# Transform for work in GEE -----------------------------------------------

catch_4326 <- catch_minus_lake %>% 
  st_transform(4326)

lakes_4326 <- lakes_sf %>% 
  st_transform(4326)

# Write to file -----------------------------------------------------------

st_write(lakes_4326, "vector/time_lakes_4326.shp")
st_write(catch_4326, "vector/time_catch_minus_lake_4326.shp")