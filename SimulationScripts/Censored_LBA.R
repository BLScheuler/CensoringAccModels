library(rtdists)
library(sft)
library(survival)
library(cmdstanr)
library(here)

# mod_cens <- cmdstan_model(file.path("../analysis/GP_LBA_diff2.stan"))
mod_norm <- cmdstan_model(
  here("SimulationScripts", "StanModels", "GP_LBA.stan")
)



{
  n_iter <- 10 # should be on par for estimating var
  A <- 1 # accumulation starting point
  b <- 1.4 # boundary threshold
  t0 <- .3 # non-decision time
  mean_v <- 1 # drift rate mean
  sd_v <- 1 # drift rate sd
  true_params <- c(A, b, t0, mean_v) # create a vector for the true parameters
  true_mat <- matrix(rep(true_params, each = n_iter), 10, 4)
}

# bias & precision
{
  bias <- c()
  prec <- c()
  bias_n <- c()
  prec_n <- c()

  # sample sizes that will be tested
  n_samp_vec <- c(50, 100, 200, 400)
  censor_vec <- .7 # seq(.7, 95, by = 0.05)
}

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

estimateKM <- function(RT, CR = NULL) {
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

# Outer Loop: different sample sizes
for (n_samples in n_samp_vec) {
  # Middle Loop: different censoring levels
  for (censoring in censor_vec) {
    ## NOTE TO BRY: try changing what will be censored! eg., 0.05 to 0.3
    # Create matrix for parameter estimates
    all_params_n <- matrix(NA, n_iter, 4)
    # Inner Loop: Repeating the actual simulation for each n_iter (estimate variance due to data)
    for (i in 1:n_iter) {
      # Use normal distributions of LBA drift rates to generate RTs
      rt <- rlba_norm(n_samples, A, b, t0, mean_v, sd_v, posdrift = TRUE)
      # Note: posdrift requires df rates to be pos

      # Pull the RTs that would be censored
      censored <- rt[, "rt"] > quantile(rt[, "rt"], censoring)
      # Replace all censored RTs with the max RT from censoring limit
      rt[censored, ] <- matrix(
        rep(
          c(quantile(rt[, "rt"], censoring), 0),
          each = sum(censored)
        ),
        nrow = sum(censored), ncol = 2
      )
      ## NOTE TO BRY: change code to censoring early tail instead of late; change < to > as needed in line 31

      tvec <- sort(unique(rt[, "rt"]))
      RT <- rt[, "rt"]
      CR <- rt[, "response"]
      # nahData <- estimateNAH(RT=RT, CR=CR)
      # Kaplan-Meier
      kmhData <- estimateKM(RT = RT, CR = CR)

      # plot(nahData$H)apply(rbind(sqrt((-apply(bias_n[, 3:6], 2, min))**2), apply(bias_n[, 3:6], 2, max)), 2, max)
      # plot(kmhData$H)

      stan_dat <- list(x = tvec, y = kmhData$H(tvec), N = length(tvec))
      # Create a subset where data is correct/ censored, and all the uncensored disappears (aka: how bad are things?)
      RT <- rt[rt[, "response"] == 1, ]

      # Pull the data we will fit with Stan model
      # stan_dat <- list(
      #   RT = RT,
      #   N = length(RT[, 1]),
      #   NUM_CHOICES = 1
      # )

      # Stan Model (see GP_LBA file)
      fit_optim2 <- mod_norm$optimize(data = stan_dat, iter = 1E6)
      # fit_optim2$print()

      # Pull parameter estimates from the Stan model
      tmp_params <- c(
        fit_optim2$summary(
          # c("A", "bMinusA", "tau", "v[1]")
          c("A", "bMinusA", "t0", "mean_v")
        )[, 2]
      )$estimate

      # Create a list  for the parameter estimates
      ## NOTE: why are we adding A & bMinusA?
      all_params_n[i, ] <- c(
        tmp_params[1], tmp_params[2] + tmp_params[1], tmp_params[3:4]
      )
    }

    # Bias: Estimated - True with means
    bias_n <- rbind(
      bias_n,
      c(
        n_samples,
        censoring,
        apply(
          all_params_n - true_mat, 2, mean
        )
      )
    )

    # Precision: Estimated - True with sd
    prec_n <- rbind(
      prec_n,
      c(
        n_samples,
        censoring,
        apply(
          all_params_n - true_mat, 2, sd
        )
      )
    )
  }
}

# Create a file for the results
save( # bias, prec,
  bias_n, prec_n,
  file = here("Data", "simstudy.Rdata")
)
