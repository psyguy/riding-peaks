# ESM data, previously cleaned
d_leuven_esm <- readRDS(here("data",
                         "d_leuven_long.rds"))

# Baseline covariates
d_leuven_cov <- read.csv(
  here("data",
       "data_leuven_complete.csv")) %>%
  select(UUID, contains("_BL")) %>%
  group_by(UUID) %>%
  mutate(id = cur_group_id(),
         .before = 1) %>%
  ungroup() %>%
  distinct() %>%
  select(-UUID)

names(d_leuven_cov) <-
  c("id",
    "age",
    "gender",
    "pa_Enthusiastic",
    "na_Anxious",
    "pa_Interested",
    "na_Angry",
    "pa_Determined",
    "na_Stressed",
    "pa_Irritated",
    "na_Jittery",
    "pa_Cheerful",
    "pa_Active",
    "pa_Strong",
    "na_Hostile",
    "na_Embarrased",
    "na_pa_Proud",
    "na_Guilty",
    "pa_Attentive",
    "pa_Inspired",
    "na_Nervous",
    "na_Sad",
    "pa_Alert")

d_leuven_cov <- d_leuven_cov %>%
  group_by(id) %>%
  mutate(bl_pa = mean(c_across(contains("pa")), na.rm = TRUE),
         bl_na = mean(c_across(contains("na")), na.rm = TRUE))

# Making a whole dataset of ESM + covariates

d_leuven_joined <- d_leuven_esm %>%
  full_join(d_leuven_cov, by = "id") %>%
  mutate(t = hour_in_day,
         co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t),
         co2 = cos(4*pi/24 * t),
         si2 = sin(4*pi/24 * t),
         .after = y)

saveRDS(d_leuven_joined,
        here("data",
             "d_leuven_joined.rds"))
