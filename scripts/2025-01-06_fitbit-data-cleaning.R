library(readr)
d_fitbit_raw <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/hourly_fitbit_sema_df_unprocessed.csv") %>%
  select(-1) %>%
  rename(id_hashed = "id")

d_fitbit_personality <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/personality.csv") %>%
  select(-1) %>%
  rename(id_hashed = "user_id") %>%
  select(-type) %>%
  group_by(id_hashed) %>%
  slice_max(order_by = submitdate,
            n = 1) %>%
  select(-submitdate)

d_fitbit_panas <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/panas.csv") %>%
  select(-1) %>%
  rename(id_hashed = "user_id") %>%
  select(-type) %>%
  rename(pa = "positive_affect_score",
         na = "negative_affect_score") %>%
  group_by(id_hashed) %>%
  slice_max(order_by = submitdate,
            n = 1) %>%
  select(-submitdate)

d_fitbit_joined <- d_fitbit_raw %>%
  select(id_hashed,
         date,
         hour,
         bpm,
         age,
         gender,
         bmi
         ) %>%
  mutate(date = lubridate::as_date(date),
         t = hour,
         y = bpm,
         co = cos(2*pi/24 * t),
         si = sin(2*pi/24 * t),
         co2 = cos(4*pi/24 * t),
         si2 = sin(4*pi/24 * t),
         .after = bpm) %>%
  group_by(id_hashed) %>%
  mutate(num_measurements = n()) %>%
  filter(num_measurements > 24) %>%
  ungroup() %>%
  mutate(id = as.integer(factor(id_hashed)),
         .before = 1) %>%
  full_join(d_fitbit_panas,
            by = "id_hashed") %>%
  full_join(d_fitbit_personality,
            by = c("id_hashed", "gender"))

d_fitbit <- d_fitbit_joined %>%
  select(id:si2) %>%
  na.omit()
d_fitbit_cov <- d_fitbit_joined %>%
  select(age) %>%
  na.omit()

item_ <- "fitbit"

