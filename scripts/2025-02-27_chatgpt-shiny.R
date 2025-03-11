library(shiny)
library(dplyr)
library(ggplot2)
library(ggforce)       # for geom_circle
library(circular)      # for circular::mean.circular

# We'll assume compute_correlations is defined externally.
# e.g. compute_correlations <- function(phi, z) { ... }

# ---------------------------------------------------------
#                 MARDIA PLOT FUNCTION
# ---------------------------------------------------------
f_plot_mard <- function(ed, point_size = 0.7) {
  # Return an empty circle if no data
  if (nrow(ed) == 0) {
    return(
      ggplot() +
        xlim(-1, 1) + ylim(-1, 1) +
        coord_fixed(1) +
        theme_minimal() +
        theme(
          panel.grid = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()
        )
    )
  }

  # Build histogram from rank-transform
  r_p  <- rank(ed$phi)
  r_z  <- rank(ed$z)
  phi_vec <- rep(2*pi * r_p / nrow(ed), times = r_z)

  bins  <- min(length(r_p), 24 * 60)
  arc   <- 2*pi / bins
  brks  <- seq(0, 2*pi, length.out = bins + 1)
  h     <- hist.default(phi_vec, breaks = brks, plot = FALSE, right = TRUE)
  mids  <- seq(arc/2, 2*pi - arc/2, length.out = bins)

  bins_count <- h$counts
  inc        <- point_size / max(bins_count, na.rm = TRUE)  # spacing outward

  d <- do.call(rbind, lapply(seq_len(bins), function(i) {
    count <- bins_count[i]
    if(count == 0) return(NULL)
    j <- seq(0, count - 1)
    r <- 1 + j * inc
    data.frame(r = r, x = r * cos(mids[i]), y = r * sin(mids[i]))
  }))

  disk_radius <- sqrt(compute_correlations(ed$phi, ed$z)$mard_cor)

  ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.3, color = "gray60") +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "gray60") +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1),
                color = "gray70", inherit.aes = FALSE) +
    geom_point(data = d, aes(x = x, y = y),
               size = point_size, color = "azure4") +
    geom_circle(aes(x0 = 0, y0 = 0, r = disk_radius),
                fill = "chartreuse3", color = NA,
                alpha = 0.7, inherit.aes = FALSE) +
    coord_fixed(1) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      axis.title      = element_blank()
    )
}

ui <- fluidPage(
  # ---------------- Row 1: Data load
  fluidRow(
    column(12,
           textAreaInput("dplyrCode",
                         "dplyr pipeline (including dataframe name and pipeline)",
                         placeholder = "e.g. myData %>% filter(...) %>% mutate(...)",
                         rows = 3, width = "90%"),
           actionButton("loadDF", "Load DF from Workspace")
    )
  ),
  tags$hr(),

  # ---------------- Row 2: Controls
  fluidRow(
    column(2, radioButtons("mode", "Mode:", choices = c("Add", "Remove"), inline = TRUE)),
    column(2, actionButton("reset", "Reset Points")),

    # Sliders for shifting phi and z, with extended widths
    column(4, sliderInput("phiShift", "Shift (phi):",
                          min = -3.14, max = 3.14, value = 0, step = 0.01,
                          width = "100%")),
    column(4, sliderInput("zShift",   "Shift (z):",
                          min = -100, max = 100, value = 0, step = 1,
                          width = "100%"))
  ),
  tags$hr(),

  # ---------------- Row 3: Main 2D plot
  fluidRow(
    column(12,
           plotOutput("plot", click = "plot_click", height = "400px")
    )
  ),
  tags$hr(),

  # ---------------- Row 4: Circular and Mardia, side by side
  fluidRow(
    column(6, plotOutput("circularPlot", click = "circularPlot_click", height = "400px")),
    column(6, plotOutput("mardiaPlot", height = "400px"))
  ),

  # ---------------- Regression summary + Download
  fluidRow(
    column(12, verbatimTextOutput("reg_summary"))
  ),
  fluidRow(
    column(12, br(), downloadButton("downloadData", "Download Data as CSV"))
  )
)

server <- function(input, output, session) {

  # Store original data (phi, z)
  rv <- reactiveValues(df = data.frame(phi = numeric(), z = numeric()))

  # -------------- Load Data --------------
  observeEvent(input$loadDF, {
    req(input$dplyrCode)
    pipeline <- input$dplyrCode
    dfName   <- sub(" %>%.*$", "", pipeline)

    if(!exists(dfName, envir = .GlobalEnv)) {
      showNotification(paste("No object named", dfName, "in global environment."), type = "error")
      return()
    }

    dfEnv <- get(dfName, envir = .GlobalEnv)
    if(!is.data.frame(dfEnv)) {
      showNotification(paste(dfName, "is not a data.frame."), type = "error")
      return()
    }

    if(nzchar(pipeline)) {
      expr <- parse(text = pipeline)
      tryCatch({
        dfRes <- eval(expr)
        rv$df <- dfRes %>% select(phi, z)
        showNotification("Data loaded successfully!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in pipeline:", e$message), type = "error")
      })
    } else {
      rv$df <- dfEnv %>% select(phi, z)
      showNotification("Data loaded successfully!", type = "message")
    }
  })

  # -------------- Effective Data --------------
  # We apply shifts to phi and z
  effectiveData <- reactive({
    df <- rv$df
    if(nrow(df) > 0) {
      # Wrap phi into [0, 2π) after shift
      df$phi_eff <- (df$phi + input$phiShift) %% (2*pi)
      df$z_eff   <- df$z + input$zShift
    }
    df
  })

  # -------------- Reactive for dynamic Y-limits on main plot --------------
  lim_min <- reactive({
    ed <- effectiveData()
    if(nrow(ed) == 0) return(-1)
    # take min(z_eff) - 5, but at least -5
    dynamic_min <- min(ed$z_eff, na.rm = TRUE) - 5
    min(-1, dynamic_min)
  })

  lim_max <- reactive({
    ed <- effectiveData()
    if(nrow(ed) == 0) return(10)
    # take max(z_eff) + 5, but at least 20
    dynamic_max <- max(ed$z_eff, na.rm = TRUE) + 5
    max(10, dynamic_max)
  })

  # -------------- Add/Remove Points: Main Plot --------------
  observeEvent(input$plot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      rv$df <- rbind(rv$df, data.frame(
        phi = (input$plot_click$x - input$phiShift) %% (2*pi),
        z   = input$plot_click$y - input$zShift
      ))
    } else {
      # For removing, we compare to the displayed (phi_eff, z_eff)
      pt <- nearPoints(ed, input$plot_click,
                       xvar = "phi_eff", yvar = "z_eff",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(
          round((rv$df$phi + input$phiShift) %% (2*pi), 5) == round(pt$phi_eff, 5) &
            round(rv$df$z + input$zShift, 5)               == round(pt$z_eff,   5)
        )[1]
        if(!is.na(idx)) rv$df <- rv$df[-idx, , drop = FALSE]
      }
    }
  })

  # -------------- Add/Remove Points: Circular Plot --------------
  observeEvent(input$circularPlot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      X_click <- input$circularPlot_click$x
      Y_click <- input$circularPlot_click$y
      new_r   <- sqrt(X_click^2 + Y_click^2)

      # We'll compare to lim_max() for allowable radius
      if(new_r <= lim_max()) {
        new_phi <- atan2(Y_click, X_click)
        if(new_phi < 0) new_phi <- new_phi + 2*pi
        # Inverse shift to store the "original" phi, z
        orig_phi <- (new_phi - input$phiShift) %% (2*pi)
        orig_z   <- new_r - input$zShift
        rv$df <- rbind(rv$df, data.frame(phi = orig_phi, z = orig_z))
      }
    } else {
      # We have the displayed (phi_eff, z_eff)
      df_cart <- ed %>% mutate(
        x_eff = z_eff * cos(phi_eff),
        y_eff = z_eff * sin(phi_eff)
      )
      pt <- nearPoints(df_cart, input$circularPlot_click,
                       xvar = "x_eff", yvar = "y_eff",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(
          round((rv$df$z + input$zShift) * cos((rv$df$phi + input$phiShift) %% (2*pi)), 5) ==
            round(pt$x_eff, 5) &
            round((rv$df$z + input$zShift) * sin((rv$df$phi + input$phiShift) %% (2*pi)), 5) ==
            round(pt$y_eff, 5)
        )[1]
        if(!is.na(idx)) rv$df <- rv$df[-idx, , drop = FALSE]
      }
    }
  })

  # -------------- Reset --------------
  observeEvent(input$reset, {
    rv$df <- data.frame(phi = numeric(), z = numeric())
  })

  # -------------- Main 2D Plot --------------
  output$plot <- renderPlot({
    ed <- effectiveData()
    if(nrow(ed) == 0) {
      ed <- data.frame(phi_eff = numeric(), z_eff = numeric())
    }

    p <- ggplot(ed, aes(x = phi_eff, y = z_eff)) +
      labs(x = "phi (shifted)", y = "z (shifted)", title = "Main 2D Plot") +
      ylim(lim_min(), lim_max())

    if(nrow(ed) == 0) {
      return(p)
    }

    p <- p + geom_point(color = "blue", size = 2)

    if(nrow(ed) < 2) {
      return(p)
    }

    # Linear: z_eff ~ phi_eff
    fit_lin <- lm(z_eff ~ phi_eff, data = ed)
    phi_grid <- seq(0, 2*pi, length.out = 200)
    df_lin   <- data.frame(phi_eff = phi_grid)
    df_lin$z_eff <- predict(fit_lin, newdata = df_lin)

    p <- p + geom_line(data = df_lin, aes(phi_eff, z_eff),
                       color = "red", size = 1)

    # Cosinor: z_eff ~ cos(phi_eff) + sin(phi_eff)
    fit_cos <- lm(z_eff ~ cos(phi_eff) + sin(phi_eff), data = ed)
    df_cos  <- data.frame(phi_eff = phi_grid)
    df_cos$z_eff <- predict(fit_cos, newdata = df_cos)
    p <- p + geom_line(data = df_cos, aes(phi_eff, z_eff),
                       color = "green", size = 1)

    # MESOR line
    coefs <- coef(fit_cos)
    M <- coefs[1]
    p <- p + geom_hline(yintercept = M, color = "purple", linetype = "dashed")

    # Correlations if ≥ 3 points
    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi_eff, ed$z_eff))
      lbl <- paste(
        sprintf("Pearson:  %.3f (p=%.3f)", corr$pearson_cor,  corr$pearson_pval),
        sprintf("Spearman: %.3f (p=%.3f)", corr$spearman_cor, corr$spearman_pval),
        sprintf("JWM:      %.3f (p=%.3f)", corr$jwm_cor,      corr$jwm_pval),
        sprintf("Mardia:   %.3f (p=%.3f)", corr$mard_cor,     corr$mard_pval),
        sep = "\n"
      )
      p <- p + annotate("label",
                        x = 0.2, y = lim_max() * 0.95,
                        label = lbl, fill = "white", alpha = 0.5,
                        hjust = 0, vjust = 1, size = 3
      )
    }

    p
  })

  # -------------- Circular Plot --------------
  output$circularPlot <- renderPlot({
    ed <- effectiveData()
    if(nrow(ed) == 0) {
      ed <- data.frame(phi_eff = numeric(), z_eff = numeric())
    }

    df_cart <- ed %>% mutate(
      x_eff = z_eff * cos(phi_eff),
      y_eff = z_eff * sin(phi_eff)
    )

    p <- ggplot(df_cart, aes(x = x_eff, y = y_eff)) +
      xlim(-lim_max(), lim_max()) +
      ylim(-lim_max(), lim_max()) +
      coord_fixed() +
      labs(title = "Circular Plot") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray")

    if(nrow(df_cart) == 0) {
      return(p)
    }

    p <- p + geom_point(size = 2, color = "blue")

    # median circle
    r_median  <- median(df_cart$z_eff)
    theta_seq <- seq(0, 2*pi, length.out = 200)
    df_median <- data.frame(
      x_eff = r_median * cos(theta_seq),
      y_eff = r_median * sin(theta_seq)
    )
    p <- p + geom_path(data = df_median, aes(x_eff, y_eff),
                       color = "red", linetype = "dashed", size = 1)

    # boundary circle
    df_bound <- data.frame(
      x_eff = lim_max() * cos(theta_seq),
      y_eff = lim_max() * sin(theta_seq)
    )
    p <- p + geom_path(data = df_bound, aes(x_eff, y_eff),
                       color = "gray80", size = 1)

    if(nrow(df_cart) < 2) {
      return(p)
    }

    # Linear fit
    fit_lin <- lm(z_eff ~ phi_eff, data = ed)
    phi_grid <- seq(0, 2*pi, length.out = 200)
    z_lin    <- predict(fit_lin, newdata = data.frame(phi_eff = phi_grid))
    df_lin   <- data.frame(
      x_eff = z_lin * cos(phi_grid),
      y_eff = z_lin * sin(phi_grid)
    )
    p <- p + geom_path(data = df_lin, aes(x_eff, y_eff),
                       color = "red", size = 1)

    # Cosinor fit
    fit_cos <- lm(z_eff ~ cos(phi_eff) + sin(phi_eff), data = ed)
    z_cos   <- predict(fit_cos, newdata = data.frame(phi_eff = phi_grid))
    df_cos  <- data.frame(
      x_eff = z_cos * cos(phi_grid),
      y_eff = z_cos * sin(phi_grid)
    )
    p <- p + geom_path(data = df_cos, aes(x_eff, y_eff),
                       color = "green", size = 1)

    # Mean angle arrow
    phase_mean <- suppressWarnings(as.numeric(circular::mean.circular(ed$phi_eff)))
    if(phase_mean < 0) phase_mean <- phase_mean + 2*pi
    p <- p + geom_segment(
      x = 0, y = 0,
      xend = r_median * cos(phase_mean),
      yend = r_median * sin(phase_mean),
      arrow = arrow(length = unit(0.2, "cm")),
      color = "blue", size = 1
    )

    # Correlations if ≥ 3 points
    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi_eff, ed$z_eff))
      lbl <- paste(
        sprintf("Pearson:  %.3f (p=%.3f)", corr$pearson_cor,  corr$pearson_pval),
        sprintf("Spearman: %.3f (p=%.3f)", corr$spearman_cor, corr$spearman_pval),
        sprintf("JWM:      %.3f (p=%.3f)", corr$jwm_cor,      corr$jwm_pval),
        sprintf("Mardia:   %.3f (p=%.3f)", corr$mard_cor,     corr$mard_pval),
        sep = "\n"
      )
      p <- p + annotate("label",
                        x = -0.95 * lim_max(), y = 0.95 * lim_max(),
                        label = lbl, fill = "white", alpha = 0.5,
                        hjust = 0, vjust = 1, size = 3
      )
    }

    p
  })

  # -------------- Mardia Plot --------------
  output$mardiaPlot <- renderPlot({
    ed <- effectiveData()
    # We'll do the Mardia plot on the shifted columns
    ed_for_mard <- ed
    ed_for_mard$phi <- ed_for_mard$phi_eff
    ed_for_mard$z   <- ed_for_mard$z_eff
    f_plot_mard(ed_for_mard[complete.cases(ed_for_mard), ])
  })

  # -------------- Regression Summary --------------
  output$reg_summary <- renderPrint({
    ed <- effectiveData()
    if(nrow(ed) < 2) {
      cat("Need at least 2 points for regression.\n")
      return()
    }

    # Linear: z_eff ~ phi_eff
    fit_lin   <- lm(z_eff ~ phi_eff, data = ed)
    lin_coefs <- coef(fit_lin)
    lin_sum   <- summary(fit_lin)
    lin_sd    <- lin_sum$sigma
    std_slope <- lin_coefs[2] * (sd(ed$phi_eff) / sd(ed$z_eff))

    # Cosinor: z_eff ~ cos(phi_eff) + sin(phi_eff)
    fit_cos   <- lm(z_eff ~ cos(phi_eff) + sin(phi_eff), data = ed)
    cos_coefs <- coef(fit_cos)
    cos_sum   <- summary(fit_cos)
    cos_sd    <- cos_sum$sigma

    M <- cos_coefs[1]
    C <- cos_coefs[2]
    S <- cos_coefs[3]
    A <- sqrt(C^2 + S^2)
    cos_phase <- atan2(S, C) %% (2*pi)

    std_C <- C * (sd(cos(ed$phi_eff)) / sd(ed$z_eff))
    std_S <- S * (sd(sin(ed$phi_eff)) / sd(ed$z_eff))

    cat("Linear Regression Parameter Estimates:\n")
    cat(sprintf("  Intercept: %.3f\n", lin_coefs[1]))
    cat(sprintf("  Slope: %.3f (standardized: %.3f)\n", lin_coefs[2], std_slope))
    cat(sprintf("  Error SD: %.3f, Variance: %.3f\n\n", lin_sd, lin_sd^2))

    cat("Cosinor Regression Parameter Estimates:\n")
    cat(sprintf("  M (MESOR): %.3f\n", M))
    cat(sprintf("  C (cos):   %.3f (std: %.3f)\n", C, std_C))
    cat(sprintf("  S (sin):   %.3f (std: %.3f)\n", S, std_S))
    cat(sprintf("  A (Amplitude): %.3f\n", A))
    cat(sprintf("  Phase: %.3f radians\n", cos_phase))
    cat(sprintf("  Error SD: %.3f, Variance: %.3f\n\n", cos_sd, cos_sd^2))

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi_eff, ed$z_eff))
      cat("Correlations (phi_eff vs. z_eff):\n")
      cat(sprintf("  Pearson:  %.3f (p=%.3f)\n",  corr$pearson_cor,  corr$pearson_pval))
      cat(sprintf("  Spearman: %.3f (p=%.3f)\n",  corr$spearman_cor, corr$spearman_pval))
      cat(sprintf("  JWM:      %.3f (p=%.3f)\n",  corr$jwm_cor,      corr$jwm_pval))
      cat(sprintf("  Mardia:   %.3f (p=%.3f)\n",  corr$mard_cor,     corr$mard_pval))
    } else {
      cat("Not enough points (need >= 3) for correlation.\n")
    }
  })

  # -------------- Download Data --------------
  output$downloadData <- downloadHandler(
    filename = function() { "my_data.csv" },
    content  = function(file) {
      write.csv(rv$df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
