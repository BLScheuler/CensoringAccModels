## Load libraries and data ----
{
  library(here)
  library(viridis)
  library(tidyverse)

  # check Cesnored_LBA, but should load:
  # 'bias', 'prec', and one day, 'bias_n' and 'prec_n'
  load(
    file = here("Data", "simstudy.Rdata")
  )
}

## plotting functions ----
{
  get_parmag <- function(parm_vec) {
    y_max <- apply(parm_vec, 2, max)
    y_min <- sqrt((-apply(parm_vec, 2, min))**2)
    y_range <- apply(rbind(y_max, y_min), 2, max)
    y_range <- round(y_range, 2) + 0.05
    return(y_range)
  }

}

## set variables
{
  parmnames <- c("A", "b-A", "t0", "v")
  parmmag <- get_parmag(bias[, 3:6]) #c(.35, .45, .2, .25)
  n_samp_vec <- unique(bias[, 1])
  censor_vec <- unique(bias[, 2])
  colvec <- viridis(length(n_samp_vec))
}


par(mfrow = c(2, 4))
for (parm in seq_along(parmnames)) {
  parm_vals <- bias[bias[, 1] == 50, 2 + parm]
  plot(
    censor_vec, parm_vals,
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm] * c(-1, 1),
    main = parmnames[parm], ylab = "Bias", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    parm_vals <- bias[bias[, 1] == b, 2 + parm]
    lines(
      censor_vec, parm_vals,
      col = colvec[i], lwd = 2
    )
  }
  abline(h = 0, col = grey(.7))
}

for (parm in 1:seq_along(parmnames)) {
  parm_vals <- bias_n[bias_n[, 1] == 50, 2 + parm]
  plot(
    censor_vec, parm_vals
    type = "l", col = colvec[1], lwd = 2,
    ylim = parmmag[parm]*c(-1,1),
    main = parmnames[parm], ylab = "Bias", xlab = "Censoring Point"
  )
  for (i in 2:4) {
    b <- n_samp_vec[i]
    lines(prec[bias_n[,1]==b,2], bias_n[bias_n[,1]==b,2+parm], col=colvec[i], lwd=2)
  }
  abline(h=0, col=grey(.7))
}

parmnames <- c("A", "b-A", "t0", "v")
parmmag <- c(.7, .7, .25, .42)
colvec <- viridis(4) #will need changed later
par(mfrow=c(2,4))
for (parm in 1:4) { 
  plot(prec[prec[,1]==50,2], prec[prec[,1]==50,2+parm], type='l', col=colvec[1], lwd=2,
       ylim=parmmag[parm]*c(0,1), main=parmnames[parm], ylab="SD", xlab="Censoring Point")
  for (i in 2:4) { 
    b = n_samp_vec[i]
    lines(prec[prec[,1]==b,2], prec[prec[,1]==b,2+parm], col=colvec[i], lwd=2)
  }
  abline(h=0, col=grey(.7))
}

for (parm in 1:4) {
  plot(prec_n[prec_n[,1]==50,2], prec_n[prec_n[,1]==50,2+parm], type='l', col=colvec[1], lwd=2,
       ylim=parmmag[parm]*c(0,1), main=parmnames[parm], ylab="SD", xlab="Censoring Point")
  for (i in 2:4) {
    b = n_samp_vec[i]
    lines(prec[prec_n[,1]==b,2], prec_n[prec_n[,1]==b,2+parm], col=colvec[i], lwd=2)
  }
  abline(h=0, col=grey(.7))
}

tvec <- sort(unique(dat$rt))[-1]
stan_dat_p <- list(x=tvec, y=Hdata$H(tvec), N=length(tvec))
fit_optim_p <- mod_cens$optimize(data=stan_dat_p, iter=1E6)
c(fit_optim_$summary(c("A", "b", "t0", "mean_v"))[,2])$estimate


RT=cbind(dat$rt[dat$rt>0],1)
stan_dat_np <- list(RT=RT, NUM_CHOICES=1, N=dim(RT)[1])
fit_optim_np <- mod_norm$optimize(data=stan_dat_np, iter=1E6)
