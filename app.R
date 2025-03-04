library(rstan)
library(ggplot2)
library(parallel)
library(MyStanPkg)
library(shiny)
library(boot)
library(bslib)
#library(styler)

# Set rstan options (optional)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# --------------------------------------------------------------------------
# UI
# --------------------------------------------------------------------------


# Define the bslib theme with navbar customization
ui <- navbarPage(
title = "Bayesian Modeling Suite",

# Tab 1: Home
tabPanel(
  "Home",
  h2("Welcome to the Bayesian Modeling Suite"),
  h3("This application provides three advanced Bayesian modeling tools for drug development analysis:"),
  br(),
  tags$ul(
    tags$li(
      h4(("Bayesian Hierarchical Modeling:"), 
      "Implements a Bayesian hierarchical model for two-group analysis, incorporating priors and group variability to improve decision-making in preclinical research."
    )),
    br(),
    tags$li(
      h4(("Bayesian Dose-Response Modeling:"), 
      "Performs Bayesian dose-response analysis, estimating target effect probabilities and visualizing uncertainty to optimize drug efficacy."
    )),
    br(),
    tags$li(
      h4(("Bayesian Toxicity Analysis Using Logistic Regression:"), 
      "Utilizes Bayesian logistic regression to estimate toxicity probability and determine the Maximum Tolerated Dose (MTD), balancing efficacy and safety in drug trials."
    )),
    br()
  ),
  h3("Use the tabs above to navigate to each tool and explore their functionalities."),
  br(),
  h3(HTML("<u>Getting Help</u>")),
  h4(
    span("All questions can be sent to Christian Dide-Agossou, PhD:", style = "color:black; display: inline-block;"),
    span(htmlOutput("homeEmail"), style = "display: inline-block;")
  )
),

# Tab 2: Hierarchical Modeling
tabPanel(
"Intervention vs Control",
titlePanel("Bayesian: Hierarchical Modeling for Two-Group Comparison"),
sidebarLayout(
sidebarPanel(
h4("1. Input Priors (e.g., percent reduction)"),
numericInput("prior_mu_drug_mean", "Prior Mean for Drug X (%)", value = 30, step = 1),
numericInput("prior_mu_drug_sd", "Prior SD for Drug X (%)", value = 10, step = 1),
numericInput("prior_mu_placebo_mean", "Prior Mean for Placebo (%)", value = 0, step = 1),
numericInput("prior_mu_placebo_sd", "Prior SD for Placebo (%)", value = 10, step = 1),
numericInput("prior_sigma_scale", "Prior Scale for Sigma (Cauchy)", value = 10, step = 1),

h4("2. Enter Observed Data"),
helpText("Enter percent reduction values (in %), comma-separated. Example: 35, 28, 42, ..."),
textAreaInput("drug_data", "Drug X Data (%)",
value = "35, 30, 45, 25, 40, 39, 32, 50, 42, 27, 34, 36, 29, 37, 41, 38, 33, 46, 48, 31",
rows = 5),
textAreaInput("placebo_data", "Placebo Data (%)",
value = "5, 15, -2, 10, 0, 8, 2, 4, 3, 7, -5, 1, 6, 9, 0, 11, -1, 12, 3, 1",
rows = 5),

h4("3. Effectiveness Threshold"),
numericInput("threshold", "Threshold for (mu_drug - mu_placebo) (%)", value = 10, step = 1),

actionButton("runAnalysis1", "Run Bayesian Analysis")
),

mainPanel(
h3(HTML("<u>Results</u>")),
verbatimTextOutput("summaryOutput1"),
plotOutput("posteriorPlot1")
)
)
),

# Tab 3: Dose-Response Modeling
tabPanel(
"Dose-Response",
titlePanel("Bayesian: Dose-Response Modeling"),
sidebarLayout(
sidebarPanel(
h4("1. Input Data"),
helpText("Enter dose levels and observed percent reduction means (comma-separated)."),
textAreaInput("dose_levels", "Dose Levels (e.g., mg/kg)", "0, 5, 10, 20, 40", rows = 2),
textAreaInput("observed_reduction", "Observed Percent Reduction (%)", "5, 15, 30, 45, 50", rows = 2),
numericInput("std_dev", "Known SD of Observations", value = 5, step = 1),

h4("2. Priors for Emax Model"),
fluidRow(
column(6, numericInput("E0_mu",  "E0 prior mean",  value = 5, step = 1)),
column(6, numericInput("E0_sd",  "E0 prior SD",    value = 5, step = 1))
),
fluidRow(
column(6, numericInput("Emax_mu", "Emax prior mean", value = 50, step = 5)),
column(6, numericInput("Emax_sd", "Emax prior SD",   value = 10, step = 1))
),
fluidRow(
column(6, numericInput("EC50_mu", "EC50 prior mean", value = 15, step = 1)),
column(6, numericInput("EC50_sd", "EC50 prior SD",   value = 5, step = 1))
),
fluidRow(
column(6, numericInput("h_mu",    "h prior mean",    value = 1, step = 0.1)),
column(6, numericInput("h_sd",    "h prior SD",      value = 0.5, step = 0.1))
),

h4("3. Target Effect Probability"),
numericInput("target_effect", "Target Reduction (%)", value = 40, step = 5),
helpText("We'll compute P(E(d) > target). For each dose, we can see which dose crosses the threshold."),

actionButton("runAnalysis2", "Run Bayesian Analysis")
),

mainPanel(
h4(HTML("<u>Posterior Summaries</u>")),
verbatimTextOutput("posteriorSummary2"),

br(),

h4(HTML("<u>Dose-Response Plot</u>")),
plotOutput("doseResponsePlot", height = "400px"),

br(),

h4(HTML("<u>Target Effect Probability per Dose</u>")),
tableOutput("targetEffectTable")
)
)
),

# Tab 4: Toxicity Modeling
tabPanel(
"Toxicity Evaluation",
titlePanel("Bayesian: Drug Toxicity Modeling"),
sidebarLayout(
sidebarPanel(
h4(HTML("<u>1. Input Data</u>")),
helpText("Enter vectors of equal length, comma-separated."),
textAreaInput("dose_levels", 
"Dose Levels (e.g., mg/kg)", 
"0, 5, 10, 20, 40", 
rows = 2),
textAreaInput("tox_events",  
"Number of Toxicities at Each Dose", 
"0, 1, 3, 4, 9", 
rows = 2),
textAreaInput("group_size",  
"Total Mice in Each Dose Group",     
"5, 5, 5, 5, 5", 
rows = 2),

h4(HTML("<u>2. Priors</u>")),
fluidRow(
column(6, numericInput("alpha_mu", "Alpha prior mean", 0, step = 1)),
column(6, numericInput("alpha_sd", "Alpha prior SD",   5, step = 1))
),
fluidRow(
column(6, numericInput("beta_mu", "Beta prior mean", 0, step = 0.5)),
column(6, numericInput("beta_sd", "Beta prior SD",   2, step = 0.5))
),

h4(HTML("<u>3. Maximum Tolerated Dose (MTD)</u>")),
numericInput("tox_threshold", 
"Toxicity Probability Threshold (e.g., 0.3 = 30%)", 
value = 0.3, 
step = 0.05),
helpText("We'll estimate the dose at which Probability(Toxicity) <= threshold."),

actionButton("runAnalysis3", "Run Bayesian Analysis")
),

mainPanel(
h4(HTML("<u>Posterior Summaries</u>")),
verbatimTextOutput("posteriorSummary3"),

h4(HTML("<u>Posterior Dose-Toxicity Curve</u>")),
plotOutput("toxCurvePlot", height = "400px"),

h4(HTML("<u>Estimated MTD</u>")),

verbatimTextOutput("mtdText")
)
)
)
)

# --------------------------------------------------------------------------
# SERVER
# --------------------------------------------------------------------------
server <- function(input, output, session) {

# Home Tab Email
output$homeEmail <- renderUI({
  email <- "christian.dideagossou@gmail.com"
  link <- paste0("mailto:", email)
  tags$a(href = link, email)  # Use tags$a to create the hyperlink
})

# Hierarchical Modeling Tab
runModel1 <- eventReactive(input$runAnalysis1, {
  # Convert text input into numeric vectors
  drug_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$drug_data, "[,\\s]+"))))
  placebo_vals <- suppressWarnings(as.numeric(unlist(strsplit(input$placebo_data, "[,\\s]+"))))
  
  # Remove NA values (caused by bad formatting)
  drug_vals <- na.omit(drug_vals)
  placebo_vals <- na.omit(placebo_vals)
  
  # Ensure at least 2 data points exist in each group
  if (length(drug_vals) < 2 || length(placebo_vals) < 2) {
    showNotification("Error: Please enter at least two values for each group.", type = "error")
    return(NULL)
  }
       
       stan_data_list1 <- list(
         N_drug = length(drug_vals),
         drug_x_data = drug_vals,
         N_placebo = length(placebo_vals),
         placebo_data = placebo_vals,
         prior_mu_drug_mean = input$prior_mu_drug_mean,
         prior_mu_drug_sd = input$prior_mu_drug_sd,
         prior_mu_placebo_mean = input$prior_mu_placebo_mean,
         prior_mu_placebo_sd = input$prior_mu_placebo_sd,
         prior_sigma_scale = input$prior_sigma_scale
       )
       
       fit <- fit_my_stan_model1(
         data = stan_data_list1,
         iter = 3000,
         warmup = 1000,
         chains = 4,
         seed = 42
       )
       
       return(fit)
})

        output$summaryOutput1 <- renderPrint({
          req(runModel1())
          fit <- runModel1()
          posterior_samples <- extract(fit)
          drug_minus_placebo <- posterior_samples$mu_drug - posterior_samples$mu_placebo
          threshold <- input$threshold
          
          prob_effective <- mean(drug_minus_placebo > threshold)
          
          cat("Posterior Summary:\n")
          print(summary(fit, probs = c(0.025, 0.5, 0.975))$summary)
          
          cat("\n-----------------------------------------------------------\n")
          cat(sprintf("Probability(Drug X - Placebo > %0.2f) = %1.2f%%\n",
                      threshold, 100 * prob_effective))
        })
        
        output$posteriorPlot1 <- renderPlot({
          req(runModel1())
          fit <- runModel1()
          posterior_samples <- extract(fit)
          
          df_all <- data.frame(
            value = c(
              posterior_samples$mu_drug,
              posterior_samples$mu_placebo,
              posterior_samples$mu_drug - posterior_samples$mu_placebo
            ),
            group = rep(
              c("mu_drug", "mu_placebo", "mu_drug - mu_placebo"),
              each = length(posterior_samples$mu_drug)
            )
          )
          
          ggplot(df_all, aes(x = value, fill = group)) +
            geom_density(alpha = 0.4) +
            facet_wrap(~ group, scales = "free") +
            theme_minimal() +
            xlab("Posterior Estimate") + ylab("Density") +
            ggtitle("Posterior Distributions")
        })
        
        # Dose-Response Modeling Tab
        runModel2 <- eventReactive(input$runAnalysis2, {
          dose_vec <- as.numeric(strsplit(input$dose_levels, "[,\\s]+")[[1]])
          resp_vec <- as.numeric(strsplit(input$observed_reduction, "[,\\s]+")[[1]])
          
          validate(
            need(length(dose_vec) == length(resp_vec), "Dose and response vectors must have the same length!")
          )
          
          stan_data_list2 <- list(
            N = length(dose_vec),
            dose = dose_vec,
            response = resp_vec,
            sigma = input$std_dev,
            prior_E0_mu    = input$E0_mu,
            prior_E0_sd    = input$E0_sd,
            prior_Emax_mu  = input$Emax_mu,
            prior_Emax_sd  = input$Emax_sd,
            prior_EC50_mu  = input$EC50_mu,
            prior_EC50_sd  = input$EC50_sd,
            prior_h_mu     = input$h_mu,
            prior_h_sd     = input$h_sd
          )
          
          fit <- fit_my_stan_model2(
            data = stan_data_list2,
            iter = 3000,
            warmup = 1000,
            chains = 4,
            seed = 42
          )
          return(list(fit = fit, dose_vec = dose_vec, resp_vec = resp_vec))
        })
        
        output$posteriorSummary2 <- renderPrint({
          req(runModel2())
          fit <- runModel2()$fit
          print(summary(fit, probs = c(0.025, 0.5, 0.975))$summary)
        })
        
        output$doseResponsePlot <- renderPlot({
          req(runModel2())
          fit      <- runModel2()$fit
          dose_vec <- runModel2()$dose_vec
          resp_vec <- runModel2()$resp_vec
          
          posterior_draws <- rstan::extract(fit)
          
          dose_grid <- seq(0, max(dose_vec)*1.2, length.out = 100)
          
          E_draws_mat <- sapply(dose_grid, function(d) {
            E0_draws   <- posterior_draws$E0
            Emax_draws <- posterior_draws$Emax
            EC50_draws <- posterior_draws$EC50
            h_draws    <- posterior_draws$h
            
            E0_draws + (Emax_draws * d^h_draws) / (EC50_draws^h_draws + d^h_draws)
          })
          
          E_mean  <- apply(E_draws_mat, 2, mean)
          E_lower <- apply(E_draws_mat, 2, quantile, probs = 0.025)
          E_upper <- apply(E_draws_mat, 2, quantile, probs = 0.975)
          
          df_plot <- data.frame(
            dose     = dose_grid,
            mean_est = E_mean,
            lower    = E_lower,
            upper    = E_upper
          )
          
          ggplot() +
            geom_ribbon(
              data = df_plot,
              aes(x = dose, ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.2
            ) +
            geom_line(
              data = df_plot,
              aes(x = dose, y = mean_est),
              color = "blue", linewidth = 1
            ) +
            geom_point(
              data = data.frame(dose = dose_vec, resp = resp_vec),
              aes(x = dose, y = resp),
              color = "red", size = 4
            ) +
            theme_minimal() +
            labs(
              x = "Dose (e.g., mg/kg)",
              y = "Percent Reduction (%)",
              title = "Posterior Mean Dose-Response Curve (with 95% Credible Interval)"
            )+
            theme_minimal(base_size = 14) +     # Increase the overall base font size
            theme(
              plot.title   = element_text(size = 16, face = "bold"), # Title size
              axis.title.x = element_text(size = 16, face = "bold"), # X-axis label size
              axis.title.y = element_text(size = 16, face = "bold"), # Y-axis label size
              axis.text.x  = element_text(size = 14),                # X tick label size
              axis.text.y  = element_text(size = 14)                 # Y tick label size
            )
        })
        
        output$targetEffectTable <- renderTable({
          req(runModel2())
          fit      <- runModel2()$fit
          dose_vec <- runModel2()$dose_vec
          posterior_draws <- rstan::extract(fit)
          
          target_val <- input$target_effect
          
          results <- lapply(seq_along(dose_vec), function(i){
            d <- dose_vec[i]
            E0_draws   <- posterior_draws$E0
            Emax_draws <- posterior_draws$Emax
            EC50_draws <- posterior_draws$EC50
            h_draws    <- posterior_draws$h
            
            Ed_draws <- E0_draws + (Emax_draws * d^h_draws) / (EC50_draws^h_draws + d^h_draws)
            prob     <- mean(Ed_draws > target_val)
            
            data.frame(
              Dose = d,
              "Probability E(d) > Target" = sprintf("%.1f%%", 100 * prob),
              check.names = FALSE
            )
          })
          
          do.call(rbind, results)
        }, digits = 0)
        
        # Toxicity Modeling Tab
        runModel3 <- eventReactive(input$runAnalysis3, {
          dose_vec <- as.numeric(strsplit(input$dose_levels, "[,\\s]+")[[1]])
          y_vec    <- as.integer(strsplit(input$tox_events,  "[,\\s]+")[[1]])
          n_vec    <- as.integer(strsplit(input$group_size,  "[,\\s]+")[[1]])
          
          validate(
            need(length(dose_vec) == length(y_vec) &&
                   length(y_vec)  == length(n_vec),
                 "Dose, #Toxicities, and Group Size vectors must all have the same length!"),
            need(all(y_vec <= n_vec),
                 "Number of toxicities cannot exceed the group size in any dose!"),
            need(all(y_vec >= 0),
                 "Toxicities must be non-negative."),
            need(all(n_vec > 0),
                 "Group size must be > 0 for each dose.")
          )
          
          K <- length(dose_vec)
          
          stan_data_list3 <- list(
            K = K,
            dose = dose_vec,
            y    = y_vec,
            n    = n_vec,
            prior_alpha_mu = input$alpha_mu,
            prior_alpha_sd = input$alpha_sd,
            prior_beta_mu  = input$beta_mu,
            prior_beta_sd  = input$beta_sd
          )
          
          fit <- fit_my_stan_model3(
            data = stan_data_list3,
            iter   = 3000,
            warmup = 1000,
            chains = 4,
            seed   = 42,
            control = list(adapt_delta = 0.95, max_treedepth = 15)
          )
          
          list(
            fit      = fit, 
            dose_vec = dose_vec,
            y_vec    = y_vec,
            n_vec    = n_vec
          )
        })
        
        output$posteriorSummary3 <- renderPrint({
          req(runModel3())
          fit_obj <- runModel3()$fit
          print(summary(fit_obj, pars = c("alpha", "beta"),
                        probs = c(0.025, 0.5, 0.975))$summary)
        })
        
        output$toxCurvePlot <- renderPlot({
          req(runModel3())
          fit_obj  <- runModel3()$fit
          dose_vec <- runModel3()$dose_vec
          y_vec    <- runModel3()$y_vec
          n_vec    <- runModel3()$n_vec
          
          draws <- rstan::extract(fit_obj)
          
          dose_grid <- seq(0, max(dose_vec)*1.2, length.out = 100)
          
          p_mat <- sapply(dose_grid, function(d){
            alpha_draws <- draws$alpha
            beta_draws  <- draws$beta
            inv.logit(alpha_draws + beta_draws * d)
          })
          
          p_mean  <- apply(p_mat, 2, mean)
          p_lower <- apply(p_mat, 2, quantile, 0.025)
          p_upper <- apply(p_mat, 2, quantile, 0.975)
          
          df_plot <- data.frame(
            dose  = dose_grid,
            p     = p_mean,
            lower = p_lower,
            upper = p_upper
          )
          
          prop_vec <- y_vec / n_vec
          obs_dat  <- data.frame(
            dose = dose_vec,
            prop = prop_vec
          )
          
          ggplot() +
            geom_ribbon(
              data = df_plot,
              aes(x = dose, ymin = lower, ymax = upper),
              fill = "blue", alpha = 0.2
            ) +
            geom_line(
              data = df_plot,
              aes(x = dose, y = p),
              color = "blue", linewidth = 1
            ) +
            geom_point(
              data = obs_dat,
              aes(x = dose, y = prop),
              color = "red", size = 3, alpha = 0.6
            ) +
            labs(
              x = "Dose (e.g., mg/kg)",
              y = "Observed/Posterior Probability of Toxicity",
              title = "Posterior Toxicity Curve (Aggregated Binomial)"
            ) +
            theme_minimal(base_size = 14)
        })
        
        output$mtdText <- renderPrint({
          req(runModel3())
          fit_obj   <- runModel3()$fit
          draws     <- rstan::extract(fit_obj)
          threshold <- input$tox_threshold
          
          alpha_draws <- draws$alpha
          beta_draws  <- draws$beta
          
          dose_draws <- (qlogis(threshold) - alpha_draws) / beta_draws
          
          mtd_median <- median(dose_draws)
          ci_lower   <- quantile(dose_draws, 0.025)
          ci_upper   <- quantile(dose_draws, 0.975)
          
          cat("Estimated MTD (median of posterior) where P(Toxicity) =",
              threshold, ":\n")
          cat(sprintf("  MTD: %.2f (e.g., mg/kg) [95%% CI: %.2f, %.2f]\n",
                      mtd_median, ci_lower, ci_upper))
        })
}

shinyApp(ui = ui, server = server)

