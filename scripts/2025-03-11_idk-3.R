ee %>%
  # filter(measure == "cor",
  #        par1 %in% c("mesor", "amp", "logsigma"),
  #        par2 %in% c("mesor", "amp", "logsigma"),
  #        correlation_type %in% c("pearson", "spearman")) %>%
  mutate(dataset = item,
         value = case_when(
           correlation_type %in% c("JWM", "Mardia") ~
             value %>% sqrt(),
           TRUE ~ value
         ),
         par_pairs = paste0("(", par1, ",", par2, ")")) %>%
  group_by(item, variable, par1, correlation_type) %>%
  summarize(point = median(value) %>% round(2),
            ci_l = quantile(value, 0.025) %>% round(2),
            ci_u = quantile(value, 0.975) %>% round(2)) %>%
  mutate(est =
           glue::glue("{point}, [{ci_l}, {ci_u}]"),
         .after = 3) %>%
  View()


   item   par1  par2     est               correlation_type point  ci_l  ci_u
   <chr>  <chr> <chr>    <glue>            <chr>            <dbl> <dbl> <dbl>
 1 fitbit amp   logsigma 0.18, [0.08, 0.2… pearson           0.18  0.08  0.27
 2 fitbit amp   logsigma 0.13, [0.05, 0.2… spearman          0.13  0.05  0.21
 3 fitbit mesor amp      -0.32, [-0.38, -… pearson          -0.32 -0.38 -0.25
 4 fitbit mesor amp      -0.32, [-0.4, -0… spearman         -0.32 -0.4  -0.25
 5 fitbit mesor logsigma 0.21, [0.16, 0.2… pearson           0.21  0.16  0.26
 6 fitbit mesor logsigma 0.28, [0.22, 0.3… spearman          0.28  0.22  0.33
 7 pa     amp   logsigma 0.19, [0, 0.43]   pearson           0.19  0     0.43
 8 pa     amp   logsigma 0.2, [0.01, 0.44] spearman          0.2   0.01  0.44
 9 pa     mesor amp      -0.01, [-0.29, 0… pearson          -0.01 -0.29  0.18
10 pa     mesor amp      -0.01, [-0.28, 0… spearman         -0.01 -0.28  0.18
11 pa     mesor logsigma -0.24, [-0.31, -… pearson          -0.24 -0.31 -0.17
12 pa     mesor logsigma -0.25, [-0.32, -… spearman         -0.25 -0.32 -0.17


item   par_pairs        correlation_type est

Dataset 2 (amp,logsigma)   pearson          0.18, [0.18, 0.18]
Dataset 2 (amp,logsigma)   spearman         0.1, [0.1, 0.1]
Dataset 2 (mesor,amp)      pearson          -0.3, [-0.3, -0.3]
Dataset 2 (mesor,amp)      spearman         -0.33, [-0.33, -0.33]
Dataset 2 (mesor,logsigma) pearson          0.23, [0.23, 0.23]
Dataset 2 (mesor,logsigma) spearman         0.26, [0.26, 0.26]
Dataset 1 (amp,logsigma)   pearson          0.42, [0.42, 0.42]
Dataset 1 (amp,logsigma)   spearman         0.43, [0.43, 0.43]
Dataset 1 (mesor,amp)      pearson          0.04, [0.04, 0.04]
Dataset 1 (mesor,amp)      spearman         0.03, [0.03, 0.03]
Dataset 1 (mesor,logsigma) pearson          -0.18, [-0.18, -0.18]
Dataset 1 (mesor,logsigma) spearman         -0.16, [-0.16, -0.16]


Dataset     Startinng time   variable         estimate
______________________________________________________________
Dataset 1  0:00   Amplitude   -0.13, [-0.34, 0.25]
Dataset 1  0:00   Amplitude   -0.11, [-0.3, 0.19]
Dataset 1  0:00   Amplitude   0.23, [0.05, 0.58]
Dataset 1  0:00   Amplitude   0.06, [0, 0.33]
Dataset 1  0:00   MESOR       0.03, [-0.17, 0.31]
Dataset 1  0:00   MESOR       0.02, [-0.16, 0.23]
Dataset 1  0:00   MESOR       0.14, [0.03, 0.42]
Dataset 1  0:00   MESOR       0.02, [0, 0.14]
Dataset 1  0:00   log(SD)     0.05, [-0.16, 0.28]
Dataset 1  0:00   log(SD)     0.03, [-0.17, 0.23]
Dataset 1  0:00   log(SD)     0.28, [0.1, 0.45]
Dataset 1  0:00   log(SD)     0.06, [0.01, 0.16]
Dataset 1  6:00   Amplitude   0.2, [-0.15, 0.38]
Dataset 1  6:00   Amplitude   0.18, [-0.22, 0.38]
Dataset 1  6:00   Amplitude   0.23, [0.05, 0.58]
Dataset 1  6:00   Amplitude   0.06, [0, 0.33]
Dataset 1  6:00   MESOR       -0.03, [-0.31, 0.2]
Dataset 1  6:00   MESOR       -0.04, [-0.39, 0.19]
Dataset 1  6:00   MESOR       0.14, [0.03, 0.42]
Dataset 1  6:00   MESOR       0.02, [0, 0.14]
Dataset 1  6:00   log(SD)     0.05, [-0.21, 0.27]
Dataset 1  6:00   log(SD)     -0.01, [-0.28, 0.23]
Dataset 1  6:00   log(SD)     0.28, [0.1, 0.45]
Dataset 1  6:00   log(SD)     0.06, [0.01, 0.16]
Dataset 1  12:00  Amplitude   0.08, [-0.27, 0.29]
Dataset 1  12:00  Amplitude   0.1, [-0.27, 0.32]
Dataset 1  12:00  Amplitude   0.23, [0.05, 0.58]
Dataset 1  12:00  Amplitude   0.06, [0, 0.33]
Dataset 1  12:00  MESOR       -0.04, [-0.34, 0.18]
Dataset 1  12:00  MESOR       -0.04, [-0.38, 0.18]
Dataset 1  12:00  MESOR       0.14, [0.03, 0.42]
Dataset 1  12:00  MESOR       0.02, [0, 0.14]
Dataset 1  12:00  log(SD)     -0.02, [-0.29, 0.2]
Dataset 1  12:00  log(SD)     -0.05, [-0.32, 0.18]
Dataset 1  12:00  log(SD)     0.28, [0.1, 0.45]
Dataset 1  12:00  log(SD)     0.06, [0.01, 0.16]
Dataset 2  0:00   Amplitude   0.14, [-0.04, 0.4]
Dataset 2  0:00   Amplitude   0.09, [-0.01, 0.24]
Dataset 2  0:00   Amplitude   0.4, [0.27, 0.52]
Dataset 2  0:00   Amplitude   0.04, [0.01, 0.09]
Dataset 2  0:00   MESOR       -0.14, [-0.26, -0.05]
Dataset 2  0:00   MESOR       0.01, [-0.09, 0.1]
Dataset 2  0:00   MESOR       0.3, [0.24, 0.36]
Dataset 2  0:00   MESOR       0.01, [0, 0.03]
Dataset 2  0:00   log(SD)     0.1, [-0.06, 0.17]
Dataset 2  0:00   log(SD)     0.15, [0.05, 0.23]
Dataset 2  0:00   log(SD)     0.25, [0.19, 0.31]
Dataset 2  0:00   log(SD)     0.05, [0.02, 0.1]
Dataset 2  6:00   Amplitude   -0.12, [-0.25, 0.05]
Dataset 2  6:00   Amplitude   -0.01, [-0.09, 0.09]
Dataset 2  6:00   Amplitude   0.4, [0.27, 0.52]
Dataset 2  6:00   Amplitude   0.04, [0.01, 0.09]
Dataset 2  6:00   MESOR       0.2, [0.12, 0.26]
Dataset 2  6:00   MESOR       0.13, [0.05, 0.2]
Dataset 2  6:00   MESOR       0.3, [0.24, 0.36]
Dataset 2  6:00   MESOR       0.01, [0, 0.03]
Dataset 2  6:00   log(SD)     0.19, [0.12, 0.25]
Dataset 2  6:00   log(SD)     0.2, [0.14, 0.27]
Dataset 2  6:00   log(SD)     0.25, [0.19, 0.31]
Dataset 2  6:00   log(SD)     0.05, [0.02, 0.1]
Dataset 2  12:00  Amplitude   -0.29, [-0.46, -0.16]
Dataset 2  12:00  Amplitude   -0.1, [-0.21, 0]
Dataset 2  12:00  Amplitude   0.4, [0.27, 0.52]
Dataset 2  12:00  Amplitude   0.04, [0.01, 0.09]
Dataset 2  12:00  MESOR       0.19, [0.12, 0.31]
Dataset 2  12:00  MESOR       0.14, [0.06, 0.23]
Dataset 2  12:00  MESOR       0.3, [0.24, 0.36]
Dataset 2  12:00  MESOR       0.01, [0, 0.03]
Dataset 2  12:00  log(SD)     0.27, [0.21, 0.33]
Dataset 2  12:00  log(SD)     0.28, [0.21, 0.36]
Dataset 2  12:00  log(SD)     0.25, [0.19, 0.31]
Dataset 2  12:00  log(SD)     0.05, [0.02, 0.1]


Dataset 1   Starting at 0:00    Amplitude   -0.13, [-0.34, 0.25]    Pearson
Dataset 1   Starting at 0:00    MESOR       0.03, [-0.17, 0.31]     Pearson
Dataset 1   Starting at 0:00    log(SD)     0.05, [-0.16, 0.28]     Pearson
Dataset 1   Starting at 6:00    Amplitude   0.2, [-0.15, 0.38]      Pearson
Dataset 1   Starting at 6:00    MESOR       -0.03, [-0.31, 0.2]     Pearson
Dataset 1   Starting at 6:00    log(SD)     0.05, [-0.21, 0.27]     Pearson
Dataset 1   Starting at 12:00   Amplitude   0.08, [-0.27, 0.29]     Pearson
Dataset 1   Starting at 12:00   MESOR       -0.04, [-0.34, 0.18]    Pearson
Dataset 1   Starting at 12:00   log(SD)     -0.02, [-0.29, 0.2]     Pearson
Dataset 1   Starting at 0:00    Amplitude   -0.11, [-0.3, 0.19]     Spearman
Dataset 1   Starting at 0:00    MESOR       0.02, [-0.16, 0.23]     Spearman
Dataset 1   Starting at 0:00    log(SD)     0.03, [-0.17, 0.23]     Spearman
Dataset 1   Starting at 6:00    Amplitude   0.18, [-0.22, 0.38]     Spearman
Dataset 1   Starting at 6:00    MESOR       -0.04, [-0.39, 0.19]    Spearman
Dataset 1   Starting at 6:00    log(SD)     -0.01, [-0.28, 0.23]    Spearman
Dataset 1   Starting at 12:00   Amplitude   0.1, [-0.27, 0.32]      Spearman
Dataset 1   Starting at 12:00   MESOR       -0.04, [-0.38, 0.18]    Spearman
Dataset 1   Starting at 12:00   log(SD)     -0.05, [-0.32, 0.18]    Spearman
Dataset 1   Starting at 0:00    Amplitude   0.48, [0.23, 0.76]      JWM
Dataset 1   Starting at 0:00    MESOR       0.38, [0.16, 0.65]      JWM
Dataset 1   Starting at 0:00    log(SD)     0.53, [0.32, 0.67]      JWM
Dataset 1   Starting at 6:00    Amplitude   0.48, [0.23, 0.76]      JWM
Dataset 1   Starting at 6:00    MESOR       0.38, [0.16, 0.65]      JWM
Dataset 1   Starting at 6:00    log(SD)     0.53, [0.32, 0.67]      JWM
Dataset 1   Starting at 12:00   Amplitude   0.48, [0.23, 0.76]      JWM
Dataset 1   Starting at 12:00   MESOR       0.38, [0.16, 0.65]      JWM
Dataset 1   Starting at 12:00   log(SD)     0.53, [0.32, 0.67]      JWM
Dataset 1   Starting at 0:00    Amplitude   0.25, [0.05, 0.57]      Mardia
Dataset 1   Starting at 0:00    MESOR       0.13, [0.02, 0.38]      Mardia
Dataset 1   Starting at 0:00    log(SD)     0.24, [0.08, 0.41]      Mardia
Dataset 1   Starting at 6:00    Amplitude   0.25, [0.05, 0.57]      Mardia
Dataset 1   Starting at 6:00    MESOR       0.13, [0.02, 0.38]      Mardia
Dataset 1   Starting at 6:00    log(SD)     0.24, [0.08, 0.41]      Mardia
Dataset 1   Starting at 12:00   Amplitude   0.25, [0.05, 0.57]      Mardia
Dataset 1   Starting at 12:00   MESOR       0.13, [0.02, 0.38]      Mardia
Dataset 1   Starting at 12:00   log(SD)     0.24, [0.08, 0.41]      Mardia
Dataset 2   Starting at 0:00    Amplitude   0.14, [-0.04, 0.4]      Pearson
Dataset 2   Starting at 0:00    MESOR       -0.14, [-0.26, -0.05]   Pearson
Dataset 2   Starting at 0:00    log(SD)     0.1, [-0.06, 0.17]      Pearson
Dataset 2   Starting at 6:00    Amplitude   -0.12, [-0.25, 0.05]    Pearson
Dataset 2   Starting at 6:00    MESOR       0.2, [0.12, 0.26]       Pearson
Dataset 2   Starting at 6:00    log(SD)     0.19, [0.12, 0.25]      Pearson
Dataset 2   Starting at 12:00   Amplitude   -0.29, [-0.46, -0.16]   Pearson
Dataset 2   Starting at 12:00   MESOR       0.19, [0.12, 0.31]      Pearson
Dataset 2   Starting at 12:00   log(SD)     0.27, [0.21, 0.33]      Pearson
Dataset 2   Starting at 0:00    Amplitude   0.09, [-0.01, 0.24]     Spearman
Dataset 2   Starting at 0:00    MESOR       0.01, [-0.09, 0.1]      Spearman
Dataset 2   Starting at 0:00    log(SD)     0.15, [0.05, 0.23]      Spearman
Dataset 2   Starting at 6:00    Amplitude   -0.01, [-0.09, 0.09]    Spearman
Dataset 2   Starting at 6:00    MESOR       0.13, [0.05, 0.2]       Spearman
Dataset 2   Starting at 6:00    log(SD)     0.2, [0.14, 0.27]       Spearman
Dataset 2   Starting at 12:00   Amplitude   -0.1, [-0.21, 0]        Spearman
Dataset 2   Starting at 12:00   MESOR       0.14, [0.06, 0.23]      Spearman
Dataset 2   Starting at 12:00   log(SD)     0.28, [0.21, 0.36]      Spearman
Dataset 2   Starting at 0:00    Amplitude   0.63, [0.52, 0.72]      JWM
Dataset 2   Starting at 0:00    MESOR       0.55, [0.49, 0.6]       JWM
Dataset 2   Starting at 0:00    log(SD)     0.5, [0.44, 0.56]       JWM
Dataset 2   Starting at 6:00    Amplitude   0.63, [0.52, 0.72]      JWM
Dataset 2   Starting at 6:00    MESOR       0.55, [0.49, 0.6]       JWM
Dataset 2   Starting at 6:00    log(SD)     0.5, [0.44, 0.56]       JWM
Dataset 2   Starting at 12:00   Amplitude   0.63, [0.52, 0.72]      JWM
Dataset 2   Starting at 12:00   MESOR       0.55, [0.49, 0.6]       JWM
Dataset 2   Starting at 12:00   log(SD)     0.5, [0.44, 0.56]       JWM
Dataset 2   Starting at 0:00    Amplitude   0.19, [0.08, 0.3]       Mardia
Dataset 2   Starting at 0:00    MESOR       0.08, [0.01, 0.17]      Mardia
Dataset 2   Starting at 0:00    log(SD)     0.23, [0.15, 0.32]      Mardia
Dataset 2   Starting at 6:00    Amplitude   0.19, [0.08, 0.3]       Mardia
Dataset 2   Starting at 6:00    MESOR       0.08, [0.01, 0.17]      Mardia
Dataset 2   Starting at 6:00    log(SD)     0.23, [0.15, 0.32]      Mardia
Dataset 2   Starting at 12:00   Amplitude   0.19, [0.08, 0.3]       Mardia
Dataset 2   Starting at 12:00   MESOR       0.08, [0.01, 0.17]      Mardia
Dataset 2   Starting at 12:00   log(SD)     0.23, [0.15, 0.32]      Mardia




Dataset 1 Pearson          MESOR       Starting at 0:00    0.03, [-0.17, 0.31]
Dataset 1 Pearson          MESOR       Starting at 6:00    -0.03, [-0.31, 0.2]
Dataset 1 Pearson          MESOR       Starting at 12:00   -0.04, [-0.34, 0.18]
Dataset 1 Pearson          Amplitude   Starting at 0:00    -0.13, [-0.34, 0.25]
Dataset 1 Pearson          Amplitude   Starting at 6:00    0.2, [-0.15, 0.38]
Dataset 1 Pearson          Amplitude   Starting at 12:00   0.08, [-0.27, 0.29]
Dataset 1 Pearson          log(SD)     Starting at 0:00    0.05, [-0.16, 0.28]
Dataset 1 Pearson          log(SD)     Starting at 6:00    0.05, [-0.21, 0.27]
Dataset 1 Pearson          log(SD)     Starting at 12:00   -0.02, [-0.29, 0.2]
Dataset 1 Spearman         MESOR       Starting at 0:00    0.02, [-0.16, 0.23]
Dataset 1 Spearman         MESOR       Starting at 6:00    -0.04, [-0.39, 0.19]
Dataset 1 Spearman         MESOR       Starting at 12:00   -0.04, [-0.38, 0.18]
Dataset 1 Spearman         Amplitude   Starting at 0:00    -0.11, [-0.3, 0.19]
Dataset 1 Spearman         Amplitude   Starting at 6:00    0.18, [-0.22, 0.38]
Dataset 1 Spearman         Amplitude   Starting at 12:00   0.1, [-0.27, 0.32]
Dataset 1 Spearman         log(SD)     Starting at 0:00    0.03, [-0.17, 0.23]
Dataset 1 Spearman         log(SD)     Starting at 6:00    -0.01, [-0.28, 0.23]
Dataset 1 Spearman         log(SD)     Starting at 12:00   -0.05, [-0.32, 0.18]
Dataset 1 JWM              MESOR       Starting at 0:00    0.14, [0.03, 0.42]
Dataset 1 JWM              MESOR       Starting at 6:00    0.14, [0.03, 0.42]
Dataset 1 JWM              MESOR       Starting at 12:00   0.14, [0.03, 0.42]
Dataset 1 JWM              Amplitude   Starting at 0:00    0.23, [0.05, 0.58]
Dataset 1 JWM              Amplitude   Starting at 6:00    0.23, [0.05, 0.58]
Dataset 1 JWM              Amplitude   Starting at 12:00   0.23, [0.05, 0.58]
Dataset 1 JWM              log(SD)     Starting at 0:00    0.28, [0.1, 0.45]
Dataset 1 JWM              log(SD)     Starting at 6:00    0.28, [0.1, 0.45]
Dataset 1 JWM              log(SD)     Starting at 12:00   0.28, [0.1, 0.45]
Dataset 1 Mardia           MESOR       Starting at 0:00    0.13, [0.02, 0.38]
Dataset 1 Mardia           MESOR       Starting at 6:00    0.13, [0.02, 0.38]
Dataset 1 Mardia           MESOR       Starting at 12:00   0.13, [0.02, 0.38]
Dataset 1 Mardia           Amplitude   Starting at 0:00    0.25, [0.05, 0.57]
Dataset 1 Mardia           Amplitude   Starting at 6:00    0.25, [0.05, 0.57]
Dataset 1 Mardia           Amplitude   Starting at 12:00   0.25, [0.05, 0.57]
Dataset 1 Mardia           log(SD)     Starting at 0:00    0.24, [0.08, 0.41]
Dataset 1 Mardia           log(SD)     Starting at 6:00    0.24, [0.08, 0.41]
Dataset 1 Mardia           log(SD)     Starting at 12:00   0.24, [0.08, 0.41]
Dataset 2 Pearson          MESOR       Starting at 0:00    -0.14, [-0.26, -0.05]
Dataset 2 Pearson          MESOR       Starting at 6:00    0.2, [0.12, 0.26]
Dataset 2 Pearson          MESOR       Starting at 12:00   0.19, [0.12, 0.31]
Dataset 2 Pearson          Amplitude   Starting at 0:00    0.14, [-0.04, 0.4]
Dataset 2 Pearson          Amplitude   Starting at 6:00    -0.12, [-0.25, 0.05]
Dataset 2 Pearson          Amplitude   Starting at 12:00   -0.29, [-0.46, -0.16]
Dataset 2 Pearson          log(SD)     Starting at 0:00    0.1, [-0.06, 0.17]
Dataset 2 Pearson          log(SD)     Starting at 6:00    0.19, [0.12, 0.25]
Dataset 2 Pearson          log(SD)     Starting at 12:00   0.27, [0.21, 0.33]
Dataset 2 Spearman         MESOR       Starting at 0:00    0.01, [-0.09, 0.1]
Dataset 2 Spearman         MESOR       Starting at 6:00    0.13, [0.05, 0.2]
Dataset 2 Spearman         MESOR       Starting at 12:00   0.14, [0.06, 0.23]
Dataset 2 Spearman         Amplitude   Starting at 0:00    0.09, [-0.01, 0.24]
Dataset 2 Spearman         Amplitude   Starting at 6:00    -0.01, [-0.09, 0.09]
Dataset 2 Spearman         Amplitude   Starting at 12:00   -0.1, [-0.21, 0]
Dataset 2 Spearman         log(SD)     Starting at 0:00    0.15, [0.05, 0.23]
Dataset 2 Spearman         log(SD)     Starting at 6:00    0.2, [0.14, 0.27]
Dataset 2 Spearman         log(SD)     Starting at 12:00   0.28, [0.21, 0.36]
Dataset 2 JWM              MESOR       Starting at 0:00    0.3, [0.24, 0.36]
Dataset 2 JWM              MESOR       Starting at 6:00    0.3, [0.24, 0.36]
Dataset 2 JWM              MESOR       Starting at 12:00   0.3, [0.24, 0.36]
Dataset 2 JWM              Amplitude   Starting at 0:00    0.4, [0.27, 0.52]
Dataset 2 JWM              Amplitude   Starting at 6:00    0.4, [0.27, 0.52]
Dataset 2 JWM              Amplitude   Starting at 12:00   0.4, [0.27, 0.52]
Dataset 2 JWM              log(SD)     Starting at 0:00    0.25, [0.19, 0.31]
Dataset 2 JWM              log(SD)     Starting at 6:00    0.25, [0.19, 0.31]
Dataset 2 JWM              log(SD)     Starting at 12:00   0.25, [0.19, 0.31]
Dataset 2 Mardia           MESOR       Starting at 0:00    0.08, [0.01, 0.17]
Dataset 2 Mardia           MESOR       Starting at 6:00    0.08, [0.01, 0.17]
Dataset 2 Mardia           MESOR       Starting at 12:00   0.08, [0.01, 0.17]
Dataset 2 Mardia           Amplitude   Starting at 0:00    0.19, [0.08, 0.3]
Dataset 2 Mardia           Amplitude   Starting at 6:00    0.19, [0.08, 0.3]
Dataset 2 Mardia           Amplitude   Starting at 12:00   0.19, [0.08, 0.3]
Dataset 2 Mardia           log(SD)     Starting at 0:00    0.23, [0.15, 0.32]
Dataset 2 Mardia           log(SD)     Starting at 6:00    0.23, [0.15, 0.32]
Dataset 2 Mardia           log(SD)     Starting at 12:00   0.23, [0.15, 0.32]
