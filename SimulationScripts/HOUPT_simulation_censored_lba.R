library(cmdstanr)
library(rtdists)
library(sft)
library(here)

print(here)
paths <- list(
  sim = here("SimulationScripts"),
  fig = here("Figures")
)
print(paths)

# mod <- cmdstan_model(file.path("../analysis/GP_LBA.stan"))
mod_cens <- cmdstan_model(file.path("../analysis/GP_LBA_diff2.stan"))
mod_norm <- cmdstan_model(file.path("../analysis/LBA.stan"))

PLOT <- FALSE

n_iter <- 10

A <- 1
b <- 1.4
t0 <- .3
mean_v <- 1
sd_v <- 1
true_params <- c(A, b, t0, mean_v)

true_mat <- matrix(rep(true_params, each = n_iter), 10, 4)

if (PLOT) {
  rt <- rlba_norm(1000, A, b, t0, mean_v, sd_v, posdrift = TRUE)
  tvec <- sort(unique(rt))
  plot(estimateNAH(RT = rt[, "rt"], CR = rt[, "response"])$H,
    do.p = FALSE, col = "red",
    xlim = c(0, quantile(rt[, "rt"], .95)),
    main = "Nelson-Aalen Estimator of\nCumulative Hazard",
    xlab = "Time", ylab = "H(t)"
  )
}

bias <- c()
bias_n <- c()
prec <- c()
prec_n <- c()

for (n_samples in c(50, 100, 200, 400)) {
  for (censoring in seq(.7, .95, by = .05)) {
    all_params <- matrix(NA, n_iter, 4)
    all_params_n <- matrix(NA, n_iter, 4)
    for (i in 1:n_iter) {
      rt <- rlba_norm(n_samples, A, b, t0, mean_v, sd_v, posdrift = TRUE)

      censored <- rt[, "rt"] > quantile(rt[, "rt"], censoring)
      rt[censored, ] <- matrix(rep(c(quantile(rt[, "rt"], censoring), 0), each = sum(censored)), sum(censored), 2)

      tvec <- sort(unique(rt[, "rt"]))
      hData <- estimateNAH(RT = rt[, "rt"], CR = rt[, "response"])

      stan_dat <- list(x = tvec, y = hData$H(tvec), N = length(tvec))

      fit_optim <- mod_cens$optimize(data = stan_dat, iter = 1E6)
      all_params[i, ] <- c(fit_optim$summary(c("A", "b", "t0", "mean_v"))[, 2])$estimate
      if (PLOT) {
        h <- c(fit_optim$summary("mu")[, 2])$estimate
        lines(tvec, hData$H(tvec), col = "purple")
        lines(tvec, h, col = grey(.9))
      }

      RT <- cbind(rt[rt[, "response"] == 1, "rt"], 1)
      stan_dat <- list(RT = RT, NUM_CHOICES = 1, N = length(RT[, 1]))
      fit_optim2 <- mod_norm$optimize(data = stan_dat, iter = 1E6)
      tmp_params <- c(fit_optim2$summary(c("A", "bMinusA", "tau", "v[1]"))[, 2])$estimate
      all_params_n[i, ] <- c(tmp_params[1], tmp_params[2] + tmp_params[1], tmp_params[3:4])
    }
    bias <- rbind(bias, c(n_samples, censoring, apply(all_params - true_mat, 2, mean)))
    prec <- rbind(prec, c(n_samples, censoring, apply(all_params - true_mat, 2, sd)))
    bias_n <- rbind(bias_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, mean)))
    prec_n <- rbind(prec_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, sd)))
  }
}

save(bias, prec, bias_n, prec_n, file = "simstudy_20240509.Rdata")


n_samp_vec <- c(50, 100, 200, 400)
parmnames <- c("A", "b-A", "t0", "v")
parmmag <- c(.35, .45, .2, .25)
colvec <- viridis(4)
par(mfrow = c(2, 4))
for (parm in 1:4) {
  plot(bias[bias[, 1] == 50, 2], bias[bias[, 1] == 50, 2 + parm],
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm] * c(-1, 1), main = parmnames[parm], ylab = "Bias", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    lines(prec[bias[, 1] == b, 2], bias[bias[, 1] == b, 2 + parm], col = colvec[i], lwd = 2)
  }
  abline(h = 0, col = grey(.7))
}

for (parm in 1:4) {
  plot(bias_n[bias_n[, 1] == 50, 2], bias_n[bias_n[, 1] == 50, 2 + parm],
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm] * c(-1, 1), main = parmnames[parm], ylab = "Bias", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    lines(prec[bias_n[, 1] == b, 2], bias_n[bias_n[, 1] == b, 2 + parm], col = colvec[i], lwd = 2)
  }
  abline(h = 0, col = grey(.7))
}



n_samp_vec <- c(50, 100, 200, 400)
parmnames <- c("A", "b-A", "t0", "v")
parmmag <- c(.7, .7, .25, .42)
colvec <- viridis(4)
par(mfrow = c(2, 4))
for (parm in 1:4) {
  plot(prec[prec[, 1] == 50, 2], prec[prec[, 1] == 50, 2 + parm],
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm] * c(0, 1), main = parmnames[parm], ylab = "SD", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    lines(prec[prec[, 1] == b, 2], prec[prec[, 1] == b, 2 + parm], col = colvec[i], lwd = 2)
  }
  abline(h = 0, col = grey(.7))
}

for (parm in 1:4) {
  plot(prec_n[prec_n[, 1] == 50, 2], prec_n[prec_n[, 1] == 50, 2 + parm],
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm] * c(0, 1), main = parmnames[parm], ylab = "SD", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    lines(prec[prec_n[, 1] == b, 2], prec_n[prec_n[, 1] == b, 2 + parm], col = colvec[i], lwd = 2)
  }
  abline(h = 0, col = grey(.7))
}

tvec <- sort(unique(dat$rt))[-1]
stan_dat_p <- list(x = tvec, y = Hdata$H(tvec), N = length(tvec))
fit_optim_p <- mod_cens$optimize(data = stan_dat_p, iter = 1E6)
c(fit_optim_$summary(c("A", "b", "t0", "mean_v"))[, 2])$estimate


RT <- cbind(dat$rt[dat$rt > 0], 1)
stan_dat_np <- list(RT = RT, NUM_CHOICES = 1, N = dim(RT)[1])
fit_optim_np <- mod_norm$optimize(data = stan_dat_np, iter = 1E6)
