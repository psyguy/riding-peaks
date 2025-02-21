library(shiny)
library(dplyr)
library(tibble)
library(Rfast)     # For Rfast::circlin.cor
library(circular)  # For circular::mean.circular

options(warn = -1)

# Helper function to draw text with a semi-transparent white background "box"
drawTextBox <- function(x, y, text, cex = 0.8, col = "black", alpha.bg = 0.5) {
  lines <- strsplit(text, "\n")[[1]]
  lineH <- strheight("M", cex = cex)
  boxH <- lineH * length(lines) * 1.2
  boxW <- max(sapply(lines, strwidth, cex = cex)) * 1.1
  rect(x, y, x + boxW, y - boxH,
       col = adjustcolor("white", alpha.f = alpha.bg),
       border = NA)
  text(x, y, labels = text, adj = c(0, 1), cex = cex, col = col)
}

ui <- fluidPage(
  # Top row: put all three inputs side by side
  fluidRow(
    # Radio buttons in first column
    column(
      width = 2,
      radioButtons("mode", "Mode:", choices = c("Add", "Remove"), inline = TRUE)
    ),
    # Reset button in second column
    column(
      width = 2,
      actionButton("reset", "Reset Points")
    ),
    # Rotation slider in third column
    column(
      width = 4,
      sliderInput("rotation", "Rotation (hours):",
                  min = -12, max = 12, value = 0, step = 0.1, width = "100%")
    )
  ),

  # Middle row: main plot (left) and circular plot (right)
  fluidRow(
    column(
      width = 8,
      plotOutput("plot", click = "plot_click", height = "400px")
    ),
    column(
      width = 4,
      plotOutput("circularPlot", click = "circularPlot_click",
                 width = "400px", height = "400px")
    )
  ),

  # Console output row
  fluidRow(
    column(12, verbatimTextOutput("reg_summary"))
  ),

  # Bottom row: download button
  fluidRow(
    column(12,
           br(),
           downloadButton("downloadData", "Download Data as CSV")
    )
  )
)

server <- function(input, output, session) {
  # Reactive data
  rv <- reactiveValues(df = data.frame(x = numeric(), y = numeric()))

  # Effective data
  effectiveData <- reactive({
    df <- rv$df
    if(nrow(df) > 0) {
      df$x_eff <- (df$x + input$rotation) %% 24
      df$phase <- df$x_eff * pi/12
      df$x_circ <- df$y * cos(df$phase)
      df$y_circ <- df$y * sin(df$phase)
    }
    df
  })

  # Main plot clicks
  observeEvent(input$plot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      newPoint <- data.frame(x = input$plot_click$x, y = input$plot_click$y)
      rv$df <- rbind(rv$df, newPoint)
    } else {
      pt <- nearPoints(ed, input$plot_click, xvar = "x_eff", yvar = "y",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(
          round(ed$x_eff, 5) == round(pt$x_eff[1], 5) &
            round(ed$y, 5) == round(pt$y[1], 5)
        )[1]
        if(!is.na(idx)) {
          rv$df <- rv$df[-idx, , drop = FALSE]
        }
      }
    }
  })

  # Circular plot clicks
  observeEvent(input$circularPlot_click, {
    ed <- effectiveData()
    if(input$mode == "Add") {
      X_click <- input$circularPlot_click$x
      Y_click <- input$circularPlot_click$y
      new_y <- sqrt(X_click^2 + Y_click^2)
      if(new_y > 100) return()
      new_phase <- atan2(Y_click, X_click)
      new_x_eff <- (new_phase * 12/pi) %% 24
      new_x <- (new_x_eff - input$rotation) %% 24
      rv$df <- rbind(rv$df, data.frame(x = new_x, y = new_y))
    } else {
      pt <- nearPoints(ed, input$circularPlot_click,
                       xvar = "x_circ", yvar = "y_circ",
                       threshold = 10, maxpoints = 1, addDist = FALSE)
      if(nrow(pt) > 0) {
        idx <- which(
          round(ed$x_circ, 5) == round(pt$x_circ[1], 5) &
            round(ed$y_circ, 5) == round(pt$y_circ[1], 5)
        )[1]
        if(!is.na(idx)) {
          rv$df <- rv$df[-idx, , drop = FALSE]
        }
      }
    }
  })

  # Reset
  observeEvent(input$reset, {
    rv$df <- data.frame(x = numeric(), y = numeric())
  })

  # Main plot
  output$plot <- renderPlot({
    ed <- effectiveData()
    plot(NA, xlim = c(0, 24), ylim = c(0, 100),
         xlab = "Peak offset (in hours)", ylab = "Y",
         main = "Main Plot")
    if(nrow(ed) == 0) return()

    points(ed$x_eff, ed$y, col = "blue", pch = 19)
    if(nrow(ed) < 2) return()

    fit_linear <- lm(y ~ x_eff, data = ed)
    abline(fit_linear, col = "red", lwd = 2)

    fit_cosinor <- lm(y ~ cos(2*pi*x_eff/24) + sin(2*pi*x_eff/24), data = ed)
    coefs <- coef(fit_cosinor)
    M <- coefs[1]
    C <- coefs[2]
    S <- coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi

    x_grid <- seq(0, 24, length.out = 200)
    y_cosinor <- M + C*cos(2*pi*x_grid/24) + S*sin(2*pi*x_grid/24)
    lines(x_grid, y_cosinor, col = "green", lwd = 2)

    abline(h = M, lty = 2, col = "purple")
    segments(x0 = psi, y0 = M, x1 = psi, y1 = M + A, col = "black", lwd = 2)

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(ed$phase, ed$y))
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

  # Circular plot (fixed 400x400, aspect ratio 1:1)
  output$circularPlot <- renderPlot({
    ed <- effectiveData()
    plot(NA, xlim = c(-100, 100), ylim = c(-100, 100),
         xlab = "x", ylab = "y", main = "Circular Plot", asp = 1)
    abline(v = 0, lty = 2, col = "gray")
    abline(h = 0, lty = 2, col = "gray")

    if(nrow(ed) == 0) return()

    phase <- ed$phase
    points(ed$y * cos(phase), ed$y * sin(phase), col = "blue", pch = 19)

    theta_seq <- seq(0, 2*pi, length.out = 100)
    r_dashed <- median(ed$y)
    lines(r_dashed * cos(theta_seq), r_dashed * sin(theta_seq), lty = 2, col = "red", lwd = 2)

    lines(100 * cos(theta_seq), 100 * sin(theta_seq),
          lty = 1, col = "gray90", lwd = 2)

    if(nrow(ed) < 2) return()

    phase_mean <- suppressWarnings(as.numeric(circular::mean.circular(phase)))
    arrows(0, 0, r_dashed * cos(phase_mean), r_dashed * sin(phase_mean),
           col = "blue", lwd = 2, length = 0.1)

    if(nrow(ed) >= 3) {
      corr <- suppressWarnings(compute_correlations(phase, ed$y))
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

  # Console output
  output$reg_summary <- renderPrint({
    ed <- effectiveData()
    if(nrow(ed) < 2) {
      cat("Need at least 2 points for regression.\n")
      return()
    }

    fit_linear <- lm(y ~ x_eff, data = ed)
    lin_coefs <- coef(fit_linear)
    lin_summary <- summary(fit_linear)
    lin_sd <- lin_summary$sigma
    lin_var <- lin_sd^2
    sd_x_eff <- sd(ed$x_eff)
    sd_y <- sd(ed$y)
    std_slope <- lin_coefs[2] * (sd_x_eff / sd_y)

    fit_cosinor <- lm(y ~ cos(2*pi*x_eff/24) + sin(2*pi*x_eff/24), data = ed)
    cos_coefs <- coef(fit_cosinor)
    M <- cos_coefs[1]
    C <- cos_coefs[2]
    S <- cos_coefs[3]
    A <- sqrt(C^2 + S^2)
    phi <- (atan2(S, C)) %% (2*pi)
    psi <- phi * 12/pi
    cos_summary <- summary(fit_cosinor)
    cos_sd <- cos_summary$sigma
    cos_var <- cos_sd^2

    x_cos <- cos(2*pi*ed$x_eff/24)
    x_sin <- sin(2*pi*ed$x_eff/24)
    sd_x_cos <- sd(x_cos)
    sd_x_sin <- sd(x_sin)
    std_C <- C * (sd_x_cos / sd_y)
    std_S <- S * (sd_x_sin / sd_y)

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
      corr <- suppressWarnings(compute_correlations(ed$phase, ed$y))
      cat("Correlations between phase (x_eff * pi/12) and y:\n")
      cat(sprintf("  Pearson: %.3f (p=%.3f)\n", corr$pearson_cor, corr$pearson_pval))
      cat(sprintf("  Spearman: %.3f (p=%.3f)\n", corr$spearman_cor, corr$spearman_pval))
      cat(sprintf("  JWM: %.3f (p=%.3f)\n", corr$jwm_cor, corr$jwm_pval))
      cat(sprintf("  Mardia: %.3f (p=%.3f)\n", corr$mard_cor, corr$mard_pval))
    } else {
      cat("Not enough points (need >= 3) for correlation.\n")
    }
  })

  # Download data as CSV
  output$downloadData <- downloadHandler(
    filename = function() {
      "my_data.csv"
    },
    content = function(file) {
      df <- rv$df
      write.csv(df, file, row.names = FALSE)
    }
  )
}

shinyApp(ui, server)
