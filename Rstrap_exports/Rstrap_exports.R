
library(Rstrap)
library(dplyr)

# Rstrap::setdet %>%
#     filter(data.series == "Campelen", rec == 6, NAFOdiv %in% c("3L", "3N", "3O")) %>%
#     group_by(spec, common.name) %>%
#     summarise(n_obs = n(), total_n = sum(number), total_weight = sum(number)) %>%
#     arrange(-total_weight) %>% View()

focal_sp <- c("Redfish spp." = 794,
              "Greenland Halibut" = 892,
              "Atlantic Cod" = 438,
              "Thorny Skate" = 90,
              "Witch Flounder" = 890,
              "American Plaice" = 889,
              "Yellowtail Flounder" = 891,
              "Atlantic Halibut" = 893,
              "Haddock" = 441,
              "Striped Wolffish" = 700,
              "Spotted Wolffish" = 701,
              "Broadhead Wolffish" = 699,
              "Common Grenadier" = 478,
              "Roughhead Grenadier" = 474,
              "Roundnose Grenadier" = 481)

## Use Rstrap to process raw set details, largely to create zeros
out <- strat.fun(setdet = Rstrap::setdet, program = "strat2",
                 data.series = "Campelen", species = focal_sp,
                 season = "fall", NAFOdiv = c("3L", "3N", "3O"),
                 export = NULL, plot.results = FALSE)

tidy_setdet <- out$raw.data$set.details
focal_sp_names <- names(focal_sp)
names(focal_sp_names) <- focal_sp
tidy_setdet$common.name <- focal_sp_names[as.character(tidy_setdet$spec)]
tidy_setdet <- tidy_setdet[, c("survey.year", "season", "vessel", "trip", "set", "NAFOdiv", "strat", "spec",
                               "common.name", "month", "day", "year", "set.dur", "dist.towed", "set.depth.mean",
                               "bot.temp", "lat.start", "long.start", "number", "weight")]
names(tidy_setdet) <- gsub("\\.", "_", names(tidy_setdet))
tidy_setdet <- tidy_setdet %>%
    rename(nafo_division = NAFOdiv, set_duration = set_dur, distance_towed = dist_towed,
           mean_set_depth = set_depth_mean, bottom_temperature = bot_temp,
           lat = lat_start, long = long_start, number = number, weight = weight,
           species_code = spec)

saveRDS(tidy_setdet, "Rstrap_exports/set_details_fall_3LNO.rds")


