## Load libraries ----
{
  library(rtdists)
  library(sft)
  library(survival)
  library(cmdstanr)
  library(here)
}

## Load Models ----
{
  # mod_cens <- cmdstan_model(
  #   here("SimulationScripts", "StanModels", "GP_LBA_diff2.stan")
  # )
  mod_cens <- cmdstan_model(
    here("SimulationScripts", "StanModels", "GP_LBA.stan")
  )
}


## Set LBA parameters ----
{
  n_iter <- 10 # should be on par for estimating var
  A <- 1 # accumulation starting point
  b <- 1.4 # boundary threshold
  t0 <- .3 # non-decision time
  mean_v <- 1 # drift rate mean
  sd_v <- 1 # drift rate sd
  true_params <- c(A, b, t0, mean_v) # create a vector for the true parameters
  true_mat <- matrix(rep(true_params, each = n_iter), n_iter, 4)
}

## Setup bias & precision values to save ----
{
  bias <- c()
  prec <- c()
  # bias_n <- c()
  # prec_n <- c()

  ## Set simulation parameters
  n_samp_vec <- c(50, 100, 200, 400)
  censor_vec <- seq(0.7, 0.95, by = 0.05)
}

## Define Hazard Functions ----
estimateNAH <- function(RT, CR = NULL) {
  nt <- length(RT)
  if (is.null(CR) | length(CR) != nt) {
    CR <- rep(1, nt)
  }
  RTx <- sort(RT, index.return = TRUE)
  RT <- RTx$x
  CR <- as.logical(CR)[RTx$ix]

  Y <- rep(NA, nt)
  for (i in 1:nt) {
    Y[i] <- sum(RT >= RT[i])
  }

  H <- stepfun(RT[CR], c(0, cumsum(1 / Y[CR])))
  H.v <- stepfun(RT[CR], c(0, cumsum(1 / Y[CR]^2)))
  return(list(H = H, Var = H.v))
}

estimateKMH <- function(RT, CR = NULL) {
  nt <- length(RT)
  if (is.null(CR) | length(CR) != nt) {
    CR <- rep(1, nt)
  }
  surv_object <- Surv(RT, CR)
  km_fit <- survfit(surv_object ~ 1)
  km_summary <- summary(km_fit)
  survival_prob <- km_summary$surv
  time_points <- km_summary$time

  cumulative_hazard <- -log(survival_prob)

  n_risk <- km_summary$n.risk
  n_events <- km_summary$n.event
  hazard_var <- cumsum(n_events / (n_risk * (n_risk - n_events)))

  H <- stepfun(time_points, c(0, cumulative_hazard))
  H_v <- stepfun(time_points, c(0, hazard_var))

  return(list(H = H, Var = H_v))
}

## Simulation loop ----
for (n_samples in n_samp_vec) {
  for (censoring in censor_vec) {
    ## NOTE TO BRY: try changing what will be censored! eg., 0.05 to 0.3
    all_params <- matrix(NA, n_iter, 4)
    for (i in 1:n_iter) {
      ## NOTE TO BRY: change code to censoring early tail instead of late;
      rt <- rlba_norm(n_samples, A, b, t0, mean_v, sd_v, posdrift = TRUE)
      censored <- rt[, "rt"] > quantile(rt[, "rt"], censoring)

      # Replace all censored RTs with the max RT from censoring limit
      rt[censored, ] <- matrix(
        rep(
          c(quantile(rt[, "rt"], censoring), 0),
          each = sum(censored)
        ),
        nrow = sum(censored), ncol = 2
      )

      tvec <- sort(unique(rt[, "rt"]))
      RT <- rt[, "rt"]
      CR <- rt[, "response"]
      # nahData <- estimateNAH(RT=RT, CR=CR)
      # Kaplan-Meier
      kmhData <- estimateKMH(RT = RT, CR = CR)
      stan_dat <- list(x = tvec, y = kmhData$H(tvec), N = length(tvec))
      fit_cens <- mod_cens$optimize(data = stan_dat, iter = 1E6)

      # Pull parameter estimates from the Stan model
      tmp_params <- c(
        fit_cens$summary(
          # c("A", "bMinusA", "tau", "v[1]")
          c("A", "bMinusA", "t0", "mean_v")
        )[, 2]
      )$estimate

      # Create a list of parameter estimates
      all_params[i, ] <- c(
        tmp_params[1], tmp_params[2] + tmp_params[1], tmp_params[3:4]
      )

      # Haven't I just done the censored data? 
      # We also do not have the below Stan model
      # Create a censored data subset with only correct responses
      # RT <- rt[rt[, "response"] == 1, ]
      # Pull the data we will fit with Stan model
      # stan_dat <- list(
      #   RT = RT,
      #   N = length(RT[, 1]),
      #   NUM_CHOICES = 1
      # )

    }

    # Bias: Estimated - True with means
    bias <- rbind(
      bias,
      c(
        n_samples,
        censoring,
        apply(
          all_params - true_mat, 2, mean
        )
      )
    )

    # Precision: Estimated - True with sd
    prec <- rbind(
      prec,
      c(
        n_samples,
        censoring,
        apply(
          all_params - true_mat, 2, sd
        )
      )
    )
  }
}

# Create a file for the results
save(bias, prec,
  file = here("Data", "simstudy.Rdata")
)
