library(shiny)
library(dplyr)
library(tibble)
library(Rfast)      # For Rfast::circlin.cor (if used in compute_correlations)
library(circular)   # For circular::mean.circular

options(warn = -1)  # Silence warnings

compute_correlations <- function(theta, x) {
  valid <- complete.cases(theta, x)
  theta <- theta[valid]
  x <- x[valid]
  n <- length(theta)
  if(n < 3) {
    return(tibble(
      pearson_cor = NA_real_, pearson_pval = NA_real_,
      spearman_cor = NA_real_, spearman_pval = NA_real_,
      jwm_cor = NA_real_, jwm_pval = NA_real_,
      mard_cor = NA_real_, mard_pval = NA_real_
    ))
  }
  # Pearson & Spearman:
  pearson <- cor(theta, x, method = "pearson")
  p_pearson <- cor.test(theta, x, method = "pearson")$p.value
  spearman <- cor(theta, x, method = "spearman")
  p_spearman <- cor.test(theta, x, method = "spearman")$p.value

  # JWM using Rfast::circlin.cor:
  jwm <- tryCatch({
    out <- Rfast::circlin.cor(theta, x)
    list(estimate = out[1] %>% sqrt(),
         p.value = out[2])
  }, error = function(e) list(estimate = NA_real_, p.value = NA_real_))

  # Mardia's circular-linear rank correlation:
  r_theta <- rank(theta, ties.method = "average")
  r_x <- rank(x, ties.method = "average")
  r_theta_star <- r_theta * 2 * pi / n
  T_c <- sum(r_x * cos(r_theta_star))
  T_s <- sum(r_x * sin(r_theta_star))
  a <- ifelse(n %% 2 == 0,
              1 / (1 + 5 / tan(pi/n)^2 + 4 / tan(pi/n)^4),
              2 * sin(pi/n)^4 / (1 + cos(pi/n))^3)
  D <- a * (T_c^2 + T_s^2)
  U_n <- 24 * (T_c^2 + T_s^2) / (n^2 * (n + 1))
  p_val <- 1 - pchisq(U_n, df = 2)

  tibble(
    pearson_cor = pearson,
    pearson_pval = p_pearson,
    spearman_cor = spearman,
    spearman_pval = p_spearman,
    jwm_cor = jwm$estimate,
    jwm_pval = jwm$p.value,
    mard_cor = D,
    mard_pval = p_val
  )
}


# Helper function to draw text with a semi-transparent white background "box"
drawTextBox <- function(x, y, text, cex = 0.8, col = "black", alpha.bg = 0.5) {
  lines <- strsplit(text, "\n")[[1]]
  lineH <- strheight("M", cex = cex)
  boxH  <- lineH * length(lines) * 1.2
  boxW  <- max(sapply(lines, strwidth, cex = cex)) * 1.1
  rect(x, y, x + boxW, y - boxH,
       col = adjustcolor("white", alpha.f = alpha.bg),
       border = NA)
  text(x, y, labels = text, adj = c(0, 1), cex = cex, col = col)
}

ui <- fluidPage(
  # Row 1: Controls for loading data from workspace
  fluidRow(
    column(width = 12,
           textAreaInput("dplyrCode", "dplyr pipeline (including dataframe name and pipeline)",
                         placeholder = "e.g. dfName %>% filter(x > 0) %>% mutate(z = x * 2)",
                         rows = 3, width = "90%")
    ),˜
    column(width = 12,
           actionButton("loadDF", "Load DF from Workspace")
    )
  ),
  tags$hr(),
  # Row 2: Controls for mode, reset, rotation slider
  fluidRow(
    column(width = 3,
           radioButtons("mode", "Mode:", choices = c("Add", "Remove"), inline = TRUE)
    ),
    column(width = 3,
           actionButton("reset", "Reset Points")
    ),
    column(width = 3,
           sliderInput("rotation", "Rotation (hours):",
                       min = -12, max = 12, value = 0, step = 0.1, width = "100%")
    ),
    column(width = 3,
           sliderInput("zshift", "Vertical Shift (z):",
                       min = -100, max = 100, value = 0, step = 1, width = "100%")
    )
  ),
  tags$hr(),
  # Row 3: Two side-by-side 2D plots (Main and Circular)
  fluidRow(
    column(width = 6,
           plotOutput("plot", click = "plot_click", height = "400px")
    ),
    column(width = 6,
           plotOutput("circularPlot", click = "circularPlot_click", height = "400px")
    )
  ),
  # Row 4: Model output text
  fluidRow(
    column(12, verbatimTextOutput("reg_summary"))
  ),
  # Row 5: Download button
  fluidRow(
    column(12, br(),
           downloadButton("downloadData", "Download Data as CSV")
    )
  )
)

server <- function(input, output, session) {
  # Reactive storage for data (must have columns: x and z)
  rv <- reactiveValues(df = data.frame(x = numeric(), z = numeric()))

  # Load dataframe from workspace when button is clicked
  observeEvent(input$loadDF, {
    req(input$dplyrCode)
    pipeline <- input$dplyrCode
    dfName <- sub(" %>%.*$", "", pipeline)  # Extract dataframe name from the input pipeline
    if (!exists(dfName, envir = .GlobalEnv)) {
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
        rv$df <- dfRes %>% select(x, z)  # Key fix: Ensure only x/z columns
        showNotification("Data loaded successfully!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in pipeline:", e$message), type = "error")
      })
    } else {
      rv$df <- dfEnv %>% select(x, z)  # Key fix: Ensure only x/z columns
      showNotification("Data loaded successfully!", type = "message")
    }
  })

  # Compute effective data, adding x_eff and phi, and apply vertical shift to z.
  effectiveData <- reactive({
    df <- rv$df
    if(nrow(df) > 0) {
      df$x_eff <- (df$x + input$rotation) %% 24
      df$z_eff <- df$z + input$zshift
      df$phi <- df$x_eff * pi/12
    }
    df
  })

  # Main Plot: Plot x_eff vs. z_eff with limits x: 0-24, z: -100 to 100.
  observeEvent(input$plot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      newPoint <- data.frame(x = input$plot_click$x, z = input$plot_click$y)
      rv$df <- rbind(rv$df, newPoint)
    } else {
      pt <- nearPoints(ed, input$plot_click, xvar = "x_eff", yvar = "z_eff",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(round(ed$x_eff, 5) == round(pt$x_eff[1], 5) &
                       round(ed$z_eff, 5) == round(pt$z_eff[1], 5))[1]
        if(!is.na(idx)) {
          rv$df <- rv$df[-idx, , drop = FALSE]
        }
      }
    }
  })

  # Circular Plot: Points at (z_eff*cos(phi), z_eff*sin(phi)); fixed limits -100 to 100.
  observeEvent(input$circularPlot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      X_click <- input$circularPlot_click$x
      Y_click <- input$circularPlot_click$y
      new_z <- sqrt(X_click^2 + Y_click^2)
      if(new_z > 100) return()  # ignore points outside 100
      new_phase <- atan2(Y_click, X_click)
      new_x_eff <- (new_phase * 12/pi) %% 24
      new_x <- (new_x_eff - input$rotation) %% 24
      rv$df <- rbind(rv$df, data.frame(x = new_x, z = new_z))
    } else {
      pt <- nearPoints(ed, input$circularPlot_click,
                       xvar = "z_eff * cos(phi)", yvar = "z_eff * sin(phi)",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(round(ed$z_eff*cos(ed$phi), 5) == round(pt[1, "x"], 5) &
                       round(ed$z_eff*sin(ed$phi), 5) == round(pt[1, "y"], 5))[1]
        if(!is.na(idx)) {
          rv$df <- rv$df[-idx, , drop = FALSE]
        }
      }
    }
  })

  # Reset points
  observeEvent(input$reset, {
    rv$df <- data.frame(x = numeric(), z = numeric())
  })

  # Main 2D Plot (Left)
  output$plot <- renderPlot({
    ed <- effectiveData()
    plot(NA, xlim = c(0, 24), ylim = c(-100, 100),
         xlab = "X (hours)", ylab = "Z", main = "Main 2D Plot")
    if(nrow(ed) == 0) return()
    points(ed$x_eff, ed$z_eff, col = "blue", pch = 19)
    if(nrow(ed) < 2) return()

    fit_lin <- lm(z_eff ~ x_eff, data = ed)
    abline(fit_lin, col = "red", lwd = 2)

    fit_cos <- lm(z_eff ~ cos(2*pi*x_eff/24) + sin(2*pi*x_eff/24), data = ed)
    coefs <- coef(fit_cos)
    M <- coefs[1]
    C <- coefs[2]
    S <- coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi

    x_grid <- seq(0,24, length.out = 200)
    z_cos <- M + C*cos(2*pi*x_grid/24) + S*sin(2*pi*x_grid/24)
    lines(x_grid, z_cos, col = "green", lwd = 2)

    abline(h = M, lty = 2, col = "purple")
    segments(x0 = psi, y0 = M, x1 = psi, y1 = M + A, col = "black", lwd = 2)

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
      txt <- paste(
        sprintf("Pearson: %.3f (p=%.3f)", corr$pearson_cor, corr$pearson_pval),
        sprintf("Spearman: %.3f (p=%.3f)", corr$spearman_cor, corr$spearman_pval),
        sprintf("JWM: %.3f (p=%.3f)", corr$jwm_cor, corr$jwm_pval),
        sprintf("Mardia: %.3f (p=%.3f)", corr$mard_cor, corr$mard_pval),
        sep="\n"
      )
      drawTextBox(1, 95, txt, cex = 0.8, col = "black", alpha.bg = 0.5)
    }
  }, bg = "white")

  # Circular 2D Plot (Right)
  output$circularPlot <- renderPlot({
    ed <- effectiveData()
    plot(NA, xlim = c(-100, 100), ylim = c(-100, 100),
         xlab = "x", ylab = "y", main = "Circular Plot", asp = 1)
    abline(v = 0, lty = 2, col = "gray")
    abline(h = 0, lty = 2, col = "gray")
    if(nrow(ed) == 0) return()

    phase <- ed$phi
    points(ed$z_eff * cos(phase), ed$z_eff * sin(phase), col = "blue", pch = 19)

    theta_seq <- seq(0, 2*pi, length.out = 100)
    r_dashed <- median(ed$z_eff)
    lines(r_dashed * cos(theta_seq), r_dashed * sin(theta_seq),
          lty = 2, col = "red", lwd = 2)

    lines(100 * cos(theta_seq), 100 * sin(theta_seq),
          lty = 1, col = "gray90", lwd = 2)

    if(nrow(ed) < 2) return()

    phase_mean <- suppressWarnings(as.numeric(circular::mean.circular(phase)))
    arrows(0, 0, r_dashed * cos(phase_mean), r_dashed * sin(phase_mean),
           col = "blue", lwd = 2, length = 0.1)

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(phase, ed$z_eff))
      txt <- paste(
        sprintf("Pearson: %.3f (p=%.3f)", corr$pearson_cor, corr$pearson_pval),
        sprintf("Spearman: %.3f (p=%.3f)", corr$spearman_cor, corr$spearman_pval),
        sprintf("JWM: %.3f (p=%.3f)", corr$jwm_cor, corr$jwm_pval),
        sprintf("Mardia: %.3f (p=%.3f)", corr$mard_cor, corr$mard_pval),
        sep="\n"
      )
      drawTextBox(-95, 95, txt, cex = 0.8, col = "black", alpha.bg = 0.5)
    }
  }, width = 400, height = 400, bg = "white")

  # Regression and correlation summary
  output$reg_summary <- renderPrint({
    ed <- effectiveData()
    if(nrow(ed) < 2) {
      cat("Need at least 2 points for regression.\n")
      return()
    }

    fit_lin <- lm(z_eff ~ x_eff, data = ed)
    lin_coefs <- coef(fit_lin)
    lin_summary <- summary(fit_lin)
    lin_sd <- lin_summary$sigma
    lin_var <- lin_sd^2
    sd_x_eff <- sd(ed$x_eff)
    sd_z <- sd(ed$z_eff)
    std_slope <- lin_coefs[2] * (sd_x_eff / sd_z)

    fit_cos <- lm(z_eff ~ cos(2*pi*x_eff/24) + sin(2*pi*x_eff/24), data = ed)
    cos_coefs <- coef(fit_cos)
    M <- cos_coefs[1]
    C <- cos_coefs[2]
    S <- cos_coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi
    cos_summary <- summary(fit_cos)
    cos_sd <- cos_summary$sigma
    cos_var <- cos_sd^2

    x_cos <- cos(2*pi*ed$x_eff/24)
    x_sin <- sin(2*pi*ed$x_eff/24)
    sd_x_cos <- sd(x_cos)
    sd_x_sin <- sd(x_sin)
    std_C <- C * (sd_x_cos / sd_z)
    std_S <- S * (sd_x_sin / sd_z)

    cat("Linear Regression Parameter Estimates:\n")
    cat(sprintf("  Intercept: %.3f\n", lin_coefs[1]))
    cat(sprintf("  Slope: %.3f (standardized: %.3f)\n", lin_coefs[2], std_slope))
    cat(sprintf("  Error SD: %.3f, Error Variance: %.3f\n\n", lin_sd, lin_var))

    cat("Cosinor Regression Parameter Estimates:\n")
    cat(sprintf("  M (MESOR): %.3f\n", M))
    cat(sprintf("  C (cosine coefficient): %.3f (standardized: %.3f)\n", C, std_C))
    cat(sprintf("  S (sine coefficient): %.3f (standardized: %.3f)\n", S, std_S))
    cat(sprintf("  A (Amplitude): %.3f\n", A))
    cat(sprintf("  ψ (Peak Offset, hours): %.3f\n", psi))
    cat(sprintf("  Error SD: %.3f, Error Variance: %.3f\n\n", cos_sd, cos_var))

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
      cat("Correlations between phase (x_eff * pi/12) and z_eff:\n")
      cat(sprintf("  Pearson: %.3f (p=%.3f)\n", corr$pearson_cor, corr$pearson_pval))
      cat(sprintf("  Spearman: %.3f (p=%.3f)\n", corr$spearman_cor, corr$spearman_pval))
      cat(sprintf("  JWM: %.3f (p=%.3f)\n", corr$jwm_cor, corr$jwm_pval))
      cat(sprintf("  Mardia: %.3f (p=%.3f)\n", corr$mard_cor, corr$mard_pval))
    } else {
      cat("Not enough points (need >= 3) for correlation.\n")
    }
  })

  # Download handler for CSV data
  output$downloadData <- downloadHandler(
    filename = function() {
      "my_data.csv"
    },
    content = function(file) {
      write.csv(rv$df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
