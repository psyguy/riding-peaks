library(readr)
d_fitbit_raw <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/hourly_fitbit_sema_df_unprocessed.csv") %>%
  select(-1) %>%
  rename(id_hashed = "id")

d_fitbit_personality <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/personality.csv") %>%
  select(-1) %>%
  rename(id_hashed = "user_id") %>%
  select(-type) %>%
  group_by(id_hashed) %>%
  slice_min(order_by = submitdate,
            n = 1) %>%
  select(-submitdate)

d_fitbit_panas <- read_csv("data/yfantidou_2022_LifeSnaps4monthMultimodal/panas.csv") %>%
  select(-1) %>%
  rename(id_hashed = "user_id") %>%
  select(-type) %>%
  rename(bl_pa = "positive_affect_score",
         bl_na = "negative_affect_score") %>%
  group_by(id_hashed) %>%
  slice_min(order_by = submitdate,
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
  mutate(num_measurements = sum(!is.na(y))) %>%
  filter(num_measurements > 24) %>%
  # Selecting up to 3 weeks of data per person
  mutate(foo = paste(id_hashed, date)) %>%
  slice_min(order_by = foo,
            n = 21*24,
            na_rm = TRUE) %>%
  select(-foo) %>%
  ungroup() %>%
  left_join(d_fitbit_panas,
            by = "id_hashed") %>%
  left_join(d_fitbit_personality,
            by = c("id_hashed", "gender")) %>%
  group_by(id_hashed) %>%
  mutate(bar = all(if_any(age:ipip_intellect_category,
                      is.na))
  ) %>%
  ungroup() %>%
  mutate(bar = paste(bar, id_hashed)) %>%
  arrange(bar) %>%
  mutate(id = dense_rank(bar),
         .before = 1) %>%
  select(-bar)

saveRDS(d_fitbit_joined,
        here("data",
             "d_fitbit_joined.rds"))
