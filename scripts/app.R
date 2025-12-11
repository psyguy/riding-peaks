library(shiny)
library(dplyr)
library(tibble)
library(DT)
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
  D_sqrt <- sqrt(a * (T_c^2 + T_s^2))
  U_n <- 24 * (T_c^2 + T_s^2) / (n^2 * (n + 1))
  p_val <- 1 - pchisq(U_n, df = 2)

  tibble(
    pearson_cor = pearson,
    pearson_pval = p_pearson,
    spearman_cor = spearman,
    spearman_pval = p_spearman,
    jwm_cor = jwm$estimate,
    jwm_pval = jwm$p.value,
    mard_cor = D_sqrt,
    mard_pval = p_val
  )
}

# Mardia correlation function
cor_mardia <- function(theta, x) {
  n <- length(theta)
  r_theta <- rank(theta, ties.method = "average")
  r_x <- rank(x, ties.method = "average")
  r_theta_star <- r_theta * 2 * pi / n
  T_c <- sum(r_x * cos(r_theta_star))
  T_s <- sum(r_x * sin(r_theta_star))
  a <- ifelse(n %% 2 == 0,
              1 / (1 + 5 / tan(pi/n)^2 + 4 / tan(pi/n)^4),
              2 * sin(pi/n)^4 / (1 + cos(pi/n))^3)
  D_sqrt <- sqrt(a * (T_c^2 + T_s^2))
  resultant_length <- sqrt(T_c^2 + T_s^2) / sum(r_x)

  list(estimate = resultant_length, statistic = D_sqrt)
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
  titlePanel("Enhanced Circular Regression Analysis Tool"),

  # Improved styling
  tags$head(
    tags$style(HTML("
      .content-wrapper { margin: 10px; }
      .well { background-color: #f8f9fa; border: 1px solid #dee2e6; }
      .btn-primary { background-color: #007bff; border-color: #007bff; }
    "))
  ),

  div(class = "content-wrapper",
      # Data Input Section
      wellPanel(
        h4("Data Input"),
        fluidRow(
          column(width = 6,
                 h5("Option 1: Load CSV File"),
                 fileInput("csvFile", "Choose CSV File",
                           accept = c(".csv"),
                           placeholder = "Select a CSV file..."),
                 helpText("CSV should contain columns 'psi' (peak time/hours) and 'z' (response variable)")
          ),
          column(width = 6,
                 h5("Option 2: Load from R Workspace"),
                 textAreaInput("dplyrCode", "dplyr pipeline (including dataframe name)",
                               placeholder = "e.g. dfName %>% filter(psi > 0) %>% mutate(z = psi * 2)",
                               rows = 3, width = "100%"),
                 actionButton("loadDF", "Load DF from Workspace", class = "btn-primary")
          )
        ),

        # Column mapping for CSV
        conditionalPanel(
          condition = "output.showColumnMapping",
          fluidRow(
            column(width = 6,
                   selectInput("psiColumn", "Select Psi column (peak time/hours):", choices = NULL)
            ),
            column(width = 6,
                   selectInput("zColumn", "Select Z column (response):", choices = NULL)
            )
          ),
          actionButton("applyMapping", "Apply Column Mapping", class = "btn-primary")
        )
      ),

      # Controls Section
      wellPanel(
        h4("Analysis Controls"),
        fluidRow(
          column(width = 3,
                 radioButtons("mode", "Interaction Mode:",
                              choices = c("Add Points" = "Add", "Remove Points" = "Remove"),
                              inline = TRUE)
          ),
          column(width = 3,
                 actionButton("reset", "Reset All Points", class = "btn btn-warning"),
                 br(), br(),
                 checkboxInput("showConfBands", "Show Confidence Bands", value = FALSE)
          ),
          column(width = 3,
                 sliderInput("rotation", "Time Rotation (hours):",
                             min = -12, max = 12, value = 0, step = 0.1, width = "100%")
          ),
          column(width = 3,
                 sliderInput("zshift", "Vertical Shift (z):",
                             min = -10, max = 10, value = 0, step = 0.01, width = "100%")
          )
        )
      ),

      # Plots Section
      wellPanel(
        h4("Visualization"),
        fluidRow(
          column(width = 8,
                 h5("Time Series Plot"),
                 plotOutput("plot", click = "plot_click", height = "400px")
          ),
          column(width = 4,
                 h5("Mardia's Rank Correlation Plot"),
                 plotOutput("mardiaPlot", height = "400px")
          )
        )
      ),

      # Results Section
      wellPanel(
        h4("Analysis Results"),
        tabsetPanel(
          tabPanel("Model Summary",
                   verbatimTextOutput("reg_summary")),
          tabPanel("Data Table",
                   DT::dataTableOutput("dataTable")),
          tabPanel("Correlation Matrix",
                   DT::dataTableOutput("corrTable"))
        )
      ),

      # Export Section
      wellPanel(
        h4("Export Results"),
        fluidRow(
          column(width = 4,
                 downloadButton("downloadData", "Download Data (CSV)", class = "btn-success")
          ),
          column(width = 4,
                 downloadButton("downloadResults", "Download Results (CSV)", class = "btn-success")
          ),
          column(width = 4,
                 downloadButton("downloadPlots", "Download Plots (PDF)", class = "btn-success")
          )
        )
      )
  )
)

server <- function(input, output, session) {
  # Reactive storage for data
  rv <- reactiveValues(
    df = data.frame(psi = numeric(), z = numeric()),
    csvData = NULL,
    columnChoices = NULL
  )

  # Show column mapping panel conditionally
  output$showColumnMapping <- reactive({
    !is.null(rv$csvData) && is.null(rv$columnChoices)
  })
  outputOptions(output, "showColumnMapping", suspendWhenHidden = FALSE)

  # Compute effective data
  effectiveData <- reactive({
    df <- rv$df
    if (nrow(df) > 0) {
      df$psi_eff <- (df$psi + input$rotation) %% 24
      df$z_eff   <-  df$z  + input$zshift
      df$phi     <-  df$psi_eff * pi/12
    }
    df
  })

  # Handle CSV file upload
  observeEvent(input$csvFile, {
    req(input$csvFile)

    tryCatch({
      df <- read.csv(input$csvFile$datapath, stringsAsFactors = FALSE)
      rv$csvData <- df

      # Auto-detect psi and z columns if they exist
      cols <- names(df)
      if("psi" %in% cols && "z" %in% cols) {
        rv$df <- df %>% select(psi, z)
        rv$columnChoices <- c("psi", "z")
        showNotification("CSV loaded successfully! Auto-detected psi and z columns.", type = "message")
      } else {
        # Show column mapping interface
        updateSelectInput(session, "psiColumn", choices = cols,
                          selected = ifelse("psi" %in% cols, "psi", cols[1]))
        updateSelectInput(session, "zColumn", choices = cols,
                          selected = ifelse("z" %in% cols, "z",
                                            ifelse(length(cols) > 1, cols[2], cols[1])))
        showNotification("Please map columns to psi (peak time) and z (response) variables.", type = "warning")
      }
    }, error = function(e) {
      showNotification(paste("Error reading CSV:", e$message), type = "error")
    })
  })

  # Apply column mapping
  observeEvent(input$applyMapping, {
    req(rv$csvData, input$psiColumn, input$zColumn)

    tryCatch({
      df <- rv$csvData
      rv$df <- data.frame(
        psi = as.numeric(df[[input$psiColumn]]),
        z = as.numeric(df[[input$zColumn]])
      )
      rv$columnChoices <- c(input$psiColumn, input$zColumn)
      showNotification("Column mapping applied successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error applying mapping:", e$message), type = "error")
    })
  })

  # Load dataframe from workspace
  observeEvent(input$loadDF, {
    req(input$dplyrCode)
    pipeline <- input$dplyrCode
    dfName <- sub(" %>%.*$", "", pipeline)

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
        if(!all(c("psi", "z") %in% names(dfRes))) {
          showNotification("Dataframe must contain columns 'psi' and 'z'", type = "error")
          return()
        }
        rv$df <- dfRes %>% select(psi, z)
        showNotification("Data loaded successfully from workspace!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error in pipeline:", e$message), type = "error")
      })
    }
  })

  # Main Plot interaction
  observeEvent(input$plot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      newPoint <- data.frame(psi = input$plot_click$x, z = input$plot_click$y)
      rv$df <- rbind(rv$df, newPoint)
    } else {
      pt <- nearPoints(ed, input$plot_click, xvar = "psi_eff", yvar = "z_eff",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(round(ed$psi_eff, 5) == round(pt$psi_eff[1], 5) &
                       round(ed$z_eff, 5) == round(pt$z_eff[1], 5))[1]
        if(!is.na(idx)) {
          rv$df <- rv$df[-idx, , drop = FALSE]
        }
      }
    }
  })

  # Reset points
  observeEvent(input$reset, {
    rv$df <- data.frame(psi = numeric(), z = numeric())
  })

  # Enhanced Main Plot
  output$plot <- renderPlot({
    ed <- effectiveData()
    plot(NA, xlim = c(0, 24), ylim = c(-10, 10),
         xlab = "Psi (hours)", ylab = "Z", main = "Time Series Plot")
    grid(col = "lightgray", lty = "dotted")

    if(nrow(ed) == 0) return()
    points(ed$psi_eff, ed$z_eff, col = "blue", pch = 19, cex = 1.2)

    if(nrow(ed) < 2) return()

    # Linear regression
    fit_lin <- lm(z_eff ~ psi_eff, data = ed)
    abline(fit_lin, col = "red", lwd = 2)

    # Confidence bands for linear regression
    if(input$showConfBands && nrow(ed) > 2) {
      psi_seq <- seq(0, 24, length.out = 100)
      pred_lin <- predict(fit_lin, newdata = data.frame(psi_eff = psi_seq),
                          interval = "confidence")
      lines(psi_seq, pred_lin[,2], col = "red", lty = 2)
      lines(psi_seq, pred_lin[,3], col = "red", lty = 2)
    }

    # Cosinor regression
    fit_cos <- lm(z_eff ~ cos(2*pi*psi_eff/24) + sin(2*pi*psi_eff/24), data = ed)
    coefs <- coef(fit_cos)
    M <- coefs[1]
    C <- coefs[2]
    S <- coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi

    psi_grid <- seq(0, 24, length.out = 200)
    z_cos <- M + C*cos(2*pi*psi_grid/24) + S*sin(2*pi*psi_grid/24)
    lines(psi_grid, z_cos, col = "green", lwd = 2)

    # Confidence bands for cosinor
    if(input$showConfBands && nrow(ed) > 3) {
      pred_cos <- predict(fit_cos,
                          newdata = data.frame(psi_eff = psi_grid),
                          interval = "confidence")
      lines(psi_grid, pred_cos[,2], col = "green", lty = 2)
      lines(psi_grid, pred_cos[,3], col = "green", lty = 2)
    }

    abline(h = M, lty = 2, col = "purple")
    segments(x0 = psi, y0 = M, x1 = psi, y1 = M + A, col = "black", lwd = 2)

    # Enhanced legend
    legend("topright",
           legend = c("Data Points", "Linear Fit", "Cosinor Fit", "MESOR", "Amplitude"),
           col = c("blue", "red", "green", "purple", "black"),
           pch = c(19, NA, NA, NA, NA),
           lty = c(NA, 1, 1, 2, 1),
           lwd = c(NA, 2, 2, 1, 2),
           bg = "white")

    # Correlation info box
    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
      txt <- paste(
        sprintf("n = %d points", nrow(ed)),
        sprintf("Pearson: r=%.3f, p=%.3f", corr$pearson_cor, corr$pearson_pval),
        sprintf("Spearman: ρ=%.3f, p=%.3f", corr$spearman_cor, corr$spearman_pval),
        sprintf("JWM: r=%.3f, p=%.3f", corr$jwm_cor, corr$jwm_pval),
        sep="\n"
      )
      drawTextBox(1, 95, txt, cex = 0.7, col = "black", alpha.bg = 0.8)
    }
  }, bg = "white")

  # Mardia's Rank Correlation Plot
  output$mardiaPlot <- renderPlot({
    ed <- effectiveData()

    if(nrow(ed) < 3) {
      plot(NA, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "", main = "Mardia's Rank Correlation Plot", asp = 1)
      text(0, 0, "Need at least 3 points\nfor Mardia correlation",
           cex = 1.2, col = "gray50")
      return()
    }

    # Build histogram from rank-transform
    r_phi <- rank(ed$phi)
    r_z <- rank(ed$z_eff)
    phi_vec <- rep(2 * pi * r_phi / nrow(ed), times = r_z)

    bins <- min(length(r_phi), 100)
    arc <- 2 * pi / bins
    brks <- seq(0, 2 * pi, length.out = bins + 1)
    h <- hist.default(phi_vec,
                      breaks = brks,
                      plot = FALSE,
                      right = TRUE)
    mids <- seq(arc / 2, 2 * pi - arc / 2, length.out = bins)

    bins_count <- h$counts
    point_size <- 2

    max_count <- max(bins_count, na.rm = TRUE)
    if(max_count == 0) {
      plot(NA, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "", main = "Mardia's Rank Correlation Plot", asp = 1)
      text(0, 0, "No valid data\nfor histogram", cex = 1.2, col = "gray50")
      return()
    }

    inc <- 0.5 * point_size / max_count

    d <- do.call(rbind, lapply(seq_len(bins), function(i) {
      count <- bins_count[i]
      if (count == 0)
        return(NULL)
      j <- seq(0, count - 1)
      r <- 1 + j * inc
      data.frame(r = r,
                 x = r * cos(mids[i]),
                 y = r * sin(mids[i]))
    }))

    if(is.null(d) || nrow(d) == 0) {
      plot(NA, xlim = c(-1, 1), ylim = c(-1, 1),
           xlab = "", ylab = "", main = "Mardia's Rank Correlation Plot", asp = 1)
      text(0, 0, "No data points\nto display", cex = 1.2, col = "gray50")
      return()
    }

    lim_mardia <- max(abs(c(d$x, d$y)), na.rm = TRUE)
    if(is.na(lim_mardia) || lim_mardia == 0) lim_mardia <- 1.5

    T_c <- sum(r_z * cos(r_phi * 2 * pi / nrow(ed)))
    T_s <- sum(r_z * sin(r_phi * 2 * pi / nrow(ed)))

    mardia_result <- tryCatch({
      cor_mardia(ed$phi, ed$z_eff)
    }, error = function(e) {
      list(estimate = NA, statistic = NA)
    })

    resultant_l <- mardia_result$estimate
    resultant_angle <- atan2(T_s, T_c)

    plot(NA, xlim = c(-lim_mardia, lim_mardia), ylim = c(-lim_mardia, lim_mardia),
         xlab = "", ylab = "", main = "Mardia's Rank Correlation Plot", asp = 1,
         axes = FALSE)

    axis(1, at = pretty(c(-lim_mardia, lim_mardia)), labels = FALSE, tcl = -0.3)
    axis(2, at = pretty(c(-lim_mardia, lim_mardia)), labels = FALSE, tcl = -0.3)

    abline(v = 0, lty = 2, col = "gray60", lwd = 0.8)
    abline(h = 0, lty = 2, col = "gray60", lwd = 0.8)

    theta_seq <- seq(0, 2*pi, length.out = 200)
    lines(cos(theta_seq), sin(theta_seq), col = "gray70", lwd = 1.5)

    points(d$x, d$y, pch = 19, cex = 0.6, col = "azure4")

    if(!is.na(resultant_l) && !is.na(resultant_angle) && resultant_l > 0) {
      arrows(0, 0,
             resultant_l * cos(resultant_angle),
             resultant_l * sin(resultant_angle),
             col = "cornflowerblue", lwd = 3, length = 0.15)

      if(resultant_l <= lim_mardia) {
        lines(resultant_l * cos(theta_seq), resultant_l * sin(theta_seq),
              lty = 2, col = "cornflowerblue", lwd = 2)
      }
    }

    if(!is.na(resultant_l)) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))

      info_lines <- c(
        sprintf("Mardia r = %.3f", resultant_l),
        sprintf("p = %.3f", ifelse(is.na(corr$mard_pval), NA, corr$mard_pval)),
        sprintf("Angle = %.1f°", resultant_angle * 180/pi %% 360),
        sprintf("n = %d", nrow(ed))
      )

      info_lines <- info_lines[!grepl("NA", info_lines)]
      info_text <- paste(info_lines, collapse = "\n")

      text_x <- lim_mardia * 0.95
      text_y <- lim_mardia * 0.95

      text_width <- max(strwidth(info_lines, cex = 0.8)) * 1.1
      text_height <- length(info_lines) * strheight("M", cex = 0.8) * 1.2

      rect(text_x - text_width, text_y - text_height, text_x, text_y,
           col = adjustcolor("white", alpha.f = 0.9), border = "gray80")

      text(text_x - text_width/2, text_y - text_height/2, info_text,
           cex = 0.8, col = "black", adj = c(0.5, 0.5))
    }

  }, bg = "white")

  # Enhanced regression summary
  output$reg_summary <- renderPrint({
    ed <- effectiveData()
    if(nrow(ed) < 2) {
      cat("Need at least 2 points for regression analysis.\n")
      cat("Current data points:", nrow(ed), "\n")
      cat("Click on plots to add points or load data from CSV/workspace.\n")
      return()
    }

    cat("=== CIRCULAR REGRESSION ANALYSIS RESULTS ===\n\n")
    cat("Data Summary:\n")
    cat(sprintf("  Number of observations: %d\n", nrow(ed)))
    cat(sprintf("  Psi range: %.2f - %.2f hours\n", min(ed$psi_eff), max(ed$psi_eff)))
    cat(sprintf("  Response range: %.2f - %.2f\n", min(ed$z_eff), max(ed$z_eff)))
    cat(sprintf("  Applied rotation: %.1f hours\n", input$rotation))
    cat(sprintf("  Applied vertical shift: %.1f\n\n", input$zshift))

    # Linear regression
    fit_lin <- lm(z_eff ~ psi_eff, data = ed)
    lin_summary <- summary(fit_lin)

    cat("LINEAR REGRESSION (z ~ psi):\n")
    cat(sprintf("  Equation: z = %.3f + %.3f * psi\n",
                coef(fit_lin)[1], coef(fit_lin)[2]))
    cat(sprintf("  R-squared: %.4f\n", lin_summary$r.squared))
    cat(sprintf("  Adjusted R-squared: %.4f\n", lin_summary$adj.r.squared))
    cat(sprintf("  RMSE: %.3f\n", lin_summary$sigma))
    cat(sprintf("  F-statistic: %.3f (p = %.4f)\n\n",
                lin_summary$fstatistic[1],
                pf(lin_summary$fstatistic[1], lin_summary$fstatistic[2],
                   lin_summary$fstatistic[3], lower.tail = FALSE)))

    # Cosinor regression
    fit_cos <- lm(z_eff ~ cos(2*pi*psi_eff/24) + sin(2*pi*psi_eff/24), data = ed)
    cos_summary <- summary(fit_cos)
    cos_coefs <- coef(fit_cos)
    M <- cos_coefs[1]
    C <- cos_coefs[2]
    S <- cos_coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi

    cat("COSINOR REGRESSION (z ~ cos(2πpsi/24) + sin(2πpsi/24)):\n")
    cat(sprintf("  MESOR (M): %.3f ± %.3f\n", M, cos_summary$coefficients[1,2]))
    cat(sprintf("  Amplitude (A): %.3f\n", A))
    cat(sprintf("  Acrophase (ψ): %.3f hours (%.1f°)\n", psi, phi * 180/pi))
    cat(sprintf("  R-squared: %.4f\n", cos_summary$r.squared))
    cat(sprintf("  Adjusted R-squared: %.4f\n", cos_summary$adj.r.squared))
    cat(sprintf("  RMSE: %.3f\n", cos_summary$sigma))
    cat(sprintf("  F-statistic: %.3f (p = %.4f)\n\n",
                cos_summary$fstatistic[1],
                pf(cos_summary$fstatistic[1], cos_summary$fstatistic[2],
                   cos_summary$fstatistic[3], lower.tail = FALSE)))


    # Correlations
    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
      cat("CIRCULAR-LINEAR CORRELATIONS:\n")
      cat(sprintf("  Pearson: r = %.3f (p = %.4f)\n", corr$pearson_cor, corr$pearson_pval))
      cat(sprintf("  Spearman: ρ = %.3f (p = %.4f)\n", corr$spearman_cor, corr$spearman_pval))
      cat(sprintf("  Johnson-Wehrly-Mardia: r = %.3f (p = %.4f)\n", corr$jwm_cor, corr$jwm_pval))
      cat(sprintf("  Mardia: D = %.3f (p = %.4f)\n", corr$mard_cor, corr$mard_pval))
    }
  })

  # Data table output
  output$dataTable <- DT::renderDataTable({
    ed <- effectiveData()
    if(nrow(ed) > 0) {
      display_data <- ed %>%
        mutate(
          psi_original = rv$df$psi,
          z_original = rv$df$z,
          psi_effective = round(psi_eff, 3),
          z_effective = round(z_eff, 3),
          phase_radians = round(phi, 3),
          phase_degrees = round(phi * 180/pi, 1)
        ) %>%
        select(psi_original, z_original, psi_effective, z_effective,
               phase_radians, phase_degrees)

      DT::datatable(display_data,
                    options = list(pageLength = 10, scrollX = TRUE),
                    colnames = c("Psi (Original)", "Z (Original)", "Psi (Effective)",
                                 "Z (Effective)", "Phase (rad)", "Phase (°)"))
    }
  })

  # Correlation table
  output$corrTable <- DT::renderDataTable({
    ed <- effectiveData()
    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
      corr_table <- data.frame(
        Method = c("Pearson", "Spearman", "Johnson-Wehrly-Mardia", "Mardia"),
        Correlation = c(corr$pearson_cor, corr$spearman_cor,
                        corr$jwm_cor, corr$mard_cor),
        P_Value = c(corr$pearson_pval, corr$spearman_pval,
                    corr$jwm_pval, corr$mard_pval),
        Significant = c(corr$pearson_pval < 0.05, corr$spearman_pval < 0.05,
                        corr$jwm_pval < 0.05, corr$mard_pval < 0.05)
      )

      DT::datatable(corr_table,
                    options = list(pageLength = 10, dom = 't')) %>%
        DT::formatRound(c("Correlation", "P_Value"), 4)
    }
  })

  # Download handlers
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("circular_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      ed <- effectiveData()
      if(nrow(ed) > 0) {
        write.csv(ed, file, row.names = FALSE)
      }
    }
  )

  output$downloadResults <- downloadHandler(
    filename = function() {
      paste0("circular_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      ed <- effectiveData()
      if(nrow(ed) >= 2) {
        # Compile results
        fit_lin <- lm(z_eff ~ psi_eff, data = ed)
        fit_cos <- lm(z_eff ~ cos(2*pi*psi_eff/24) + sin(2*pi*psi_eff/24), data = ed)

        lin_summary <- summary(fit_lin)
        cos_summary <- summary(fit_cos)
        cos_coefs <- coef(fit_cos)

        M <- cos_coefs[1]
        C <- cos_coefs[2]
        S <- cos_coefs[3]
        A <- sqrt(C^2 + S^2)
        phi <- (atan2(S, C)) %% (2*pi)
        psi <- phi * 12/pi

        results <- data.frame(
          Parameter = c("Sample_Size", "Linear_Intercept", "Linear_Slope",
                        "Linear_R_Squared", "Linear_RMSE", "Linear_P_Value",
                        "Cosinor_MESOR", "Cosinor_Amplitude", "Cosinor_Acrophase_Hours",
                        "Cosinor_Acrophase_Degrees", "Cosinor_R_Squared", "Cosinor_RMSE",
                        "Cosinor_P_Value", "Linear_AIC", "Cosinor_AIC"),
          Value = c(nrow(ed), coef(fit_lin)[1], coef(fit_lin)[2],
                    lin_summary$r.squared, lin_summary$sigma,
                    pf(lin_summary$fstatistic[1], lin_summary$fstatistic[2],
                       lin_summary$fstatistic[3], lower.tail = FALSE),
                    M, A, psi, phi * 180/pi, cos_summary$r.squared, cos_summary$sigma,
                    pf(cos_summary$fstatistic[1], cos_summary$fstatistic[2],
                       cos_summary$fstatistic[3], lower.tail = FALSE),
                    AIC(fit_lin), AIC(fit_cos))
        )

        # Add correlations if available
        if(nrow(ed) >= 3) {
          corr <- suppressWarnings(compute_correlations(ed$phi, ed$z_eff))
          corr_results <- data.frame(
            Parameter = c("Pearson_Correlation", "Pearson_P_Value",
                          "Spearman_Correlation", "Spearman_P_Value",
                          "JWM_Correlation", "JWM_P_Value",
                          "Mardia_Correlation", "Mardia_P_Value"),
            Value = c(corr$pearson_cor, corr$pearson_pval,
                      corr$spearman_cor, corr$spearman_pval,
                      corr$jwm_cor, corr$jwm_pval,
                      corr$mard_cor, corr$mard_pval)
          )
          results <- rbind(results, corr_results)
        }

        write.csv(results, file, row.names = FALSE)
      }
    }
  )

  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste0("circular_plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ed <- effectiveData()
      if(nrow(ed) > 0) {
        pdf(file, width = 12, height = 8)

        # Create a 2x1 layout
        par(mfrow = c(1, 2), mar = c(4, 4, 3, 2))

        # Main plot
        plot(NA, xlim = c(0, 24), ylim = c(-10 , 10),
             xlab = "Psi (hours)", ylab = "Z", main = "Time Series Plot")
        grid(col = "lightgray", lty = "dotted")
        points(ed$psi_eff, ed$z_eff, col = "blue", pch = 19, cex = 1.2)

        if(nrow(ed) >= 2) {
          fit_lin <- lm(z_eff ~ psi_eff, data = ed)
          abline(fit_lin, col = "red", lwd = 2)

          fit_cos <- lm(z_eff ~ cos(2*pi*psi_eff/24) + sin(2*pi*psi_eff/24), data = ed)
          coefs <- coef(fit_cos)
          M <- coefs[1]
          C <- coefs[2]
          S <- coefs[3]
          A <- sqrt(C^2 + S^2)
          phi <- (atan2(S, C)) %% (2*pi)
          psi <- phi * 12/pi

          psi_grid <- seq(0, 24, length.out = 200)
          z_cos <- M + C*cos(2*pi*psi_grid/24) + S*sin(2*pi*psi_grid/24)
          lines(psi_grid, z_cos, col = "green", lwd = 2)

          abline(h = M, lty = 2, col = "purple")
          segments(x0 = psi, y0 = M, x1 = psi, y1 = M + A, col = "black", lwd = 2)

          legend("topright",
                 legend = c("Data", "Linear", "Cosinor", "MESOR", "Amplitude"),
                 col = c("blue", "red", "green", "purple", "black"),
                 pch = c(19, NA, NA, NA, NA),
                 lty = c(NA, 1, 1, 2, 1), lwd = c(NA, 2, 2, 1, 2))
        }

        # Mardia plot
        plot(NA, xlim = c(-10, 10), ylim = c(-10, 10),
             xlab = "X Component", ylab = "Y Component",
             main = "Mardia's Rank Correlation Plot", asp = 1)

        if(nrow(ed) >= 3) {
          # Build histogram from rank-transform (simplified for PDF)
          r_phi <- rank(ed$phi)
          r_z <- rank(ed$z_eff)
          phi_vec <- rep(2 * pi * r_phi / nrow(ed), times = r_z)

          bins <- min(length(r_phi), 100)
          arc <- 2 * pi / bins
          brks <- seq(0, 2 * pi, length.out = bins + 1)
          h <- hist.default(phi_vec, breaks = brks, plot = FALSE, right = TRUE)
          mids <- seq(arc / 2, 2 * pi - arc / 2, length.out = bins)

          bins_count <- h$counts
          point_size <- 1.5
          inc <- 0.5 * point_size / max(bins_count, na.rm = TRUE)

          d <- do.call(rbind, lapply(seq_len(bins), function(i) {
            count <- bins_count[i]
            if (count == 0) return(NULL)
            j <- seq(0, count - 1)
            r <- 1 + j * inc
            data.frame(r = r, x = r * cos(mids[i]), y = r * sin(mids[i]))
          }))

          if(!is.null(d) && nrow(d) > 0) {
            theta_seq <- seq(0, 2*pi, length.out = 100)
            lines(cos(theta_seq), sin(theta_seq), col = "gray70", lwd = 1)
            abline(v = 0, lty = 2, col = "gray60")
            abline(h = 0, lty = 2, col = "gray60")

            points(d$x, d$y, pch = 19, cex = 0.5, col = "azure4")

            T_c <- sum(r_z * cos(r_phi * 2 * pi / nrow(ed)))
            T_s <- sum(r_z * sin(r_phi * 2 * pi / nrow(ed)))
            mardia_result <- cor_mardia(ed$phi, ed$z_eff)
            resultant_l <- mardia_result$estimate
            resultant_angle <- atan2(T_s, T_c)

            if(!is.na(resultant_l) && !is.na(resultant_angle)) {
              arrows(0, 0, resultant_l * cos(resultant_angle),
                     resultant_l * sin(resultant_angle),
                     col = "cornflowerblue", lwd = 2, length = 0.1)
            }
          }
        } else {
          text(0, 0, "Need ≥3 points\nfor Mardia plot", cex = 1, col = "gray50")
        }

        dev.off()
      }
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
