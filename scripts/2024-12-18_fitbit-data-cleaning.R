library(readr)
d_fitbit_raw <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/hourly_fitbit_sema_df_unprocessed.csv") %>%
  select(-1)


## There are fewer missingness
## from 2021-06-14 to 2021-07-25
d_fitbit <- d_fitbit_raw %>%
  select(id,
         date,
         hour,
         temperature,
         bpm,
         age,
         gender,
         bmi) %>%
  mutate(date = lubridate::as_date(date),
         t = hour,
         y = bpm) %>%
  # Selecting three weeks of data
  filter(date >= "2021-07-01",
         date <= "2021-07-21") %>%
  # Selecting 40 users that have exactly 24*21 measurements
  filter(!(id %in% c("621e2ef567b776a24099f889",
                   "621e323667b776a240f19134",
                   "621e2f9167b776a240011ccb",
                   "621e2f1b67b776a240b3d87c",
                   "621e301e67b776a240608a72"))) %>%
  mutate(co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t),
         co2 = cos(4*pi/24 * t),
         si2 = sin(4*pi/24 * t))


