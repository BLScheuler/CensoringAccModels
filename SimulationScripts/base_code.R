library(cmdstanr)
library(rtdists)

CENS <- FALSE
NOCENS <- TRUE

mod_norm <- cmdstan_model(file.path("StanModels/LBA_censored.stan"))

PLOT <- FALSE

n_iter <- 10

A <- 1
b <- 1.4
t0 <- .3
mean_v <- 1
sd_v <- 1
true_params <- c(A, b, t0, mean_v)

true_mat <- matrix(rep(true_params, each=n_iter), 10, 4)

if (CENS) { 
  bias_cens <- c()
  prec_cens <- c()
}
if (NOCENS) { 
  bias_nocens <- c()
  prec_nocens <- c()
}

for (n_samples in c(50, 100, 200, 400)) { 
  for (censoring in seq(.7, .95, by=.05)) { 
      
    if (CENS) params_cens <- matrix(NA, n_iter, 4)
    if (NOCENS) params_nocens <- matrix(NA, n_iter, 4)
    for (i in 1:n_iter) { 
      # Generate data directly from the model
      rt <- rlba_norm(n_samples, A=true_mat[i,1], 
		                 b=true_mat[i,2],  
		                 t0=true_mat[i,3],  
		                 mean_v=true_mat[i,4],  
				 sd_v=sd_v, posdrift = TRUE)
      
      stan_data_raw <- list(
          NUM_CHOICES = 1,
	  N = dim(rt)[1],
	  upper_cutoff = max(rt[,1]) + 1,
	  lower_cutoff = 0,
	  RT = rt
      )

      # Censor data
      upper_cutoff <- quantile(rt[,"rt"], censoring)
      censored <- rt[,"rt"] > upper_cutoff
      rt_c <- rt
      rt_c[censored,] <- matrix(rep(c(quantile(rt[,"rt"], censoring), 0), each=sum(censored)), sum(censored), 2)

      stan_data_cens <- list(
          NUM_CHOICES = 1,
	  N = dim(rt)[1],
	  upper_cutoff = upper_cutoff,
	  lower_cutoff = 0,
	  RT = rt_c
      )

      # Fit model that accounts for censoring
      if (CENS) { 
        fit_cens <- mod_norm$sample(data=stan_data_cens, chains=4, parallel_chains=4)
        params_raw <- fit_cens$summary(c("bMinusA", "A", "tau", "v_std[1]"))[,2]$mean 
	names(params_raw) <- c("bMinusA", "A", "tau", "v_std[1]")
	params <- rep(NA, 4); names(params) <- c("A", "b", "t0", "mean_v")
	params["A"] <- params_raw["A"]
	params["b"] <- params_raw["bMinusA"] + params_raw["A"]
	params["t0"] <- params_raw["tau"]
	params["mean_v"] <- 1+.1*params_raw["v_std[1]"]
        params_cens[i,] <- params
      }
    
      # Fit model that ignores censoring NOT RUN!!!!
      if (NOCENS) { 
        RT=cbind(rt[rt[,"response"]==1,"rt"], 1)
        dat_nocens <- list(RT=RT, NUM_CHOICES=1, N=length(RT[,1]))
        fit_nocens <- mod_norm$sample(data=dat_nocens)
        params <- c(fit_optim2$summary(c("A", "bMinusA", "tau", "v[1]"))[,2])$estimate
        params_nocens[i,] <- c(params[1], params[2] +params[1], params[3:4])
      }
    }
    if (CENS) { 
      bias_cens <- rbind(bias_cens, c(n_samples, censoring, apply(params_cens - true_mat, 2, mean)))
      prec_cens <- rbind(prec_cens, c(n_samples, censoring, apply(params_cens - true_mat, 2, sd)))
    }
    if (NOCENS) { 
      bias_nocens <- rbind(bias_nocens, c(n_samples, censoring, apply(params_nocens - true_mat, 2, mean)))
      prec_nocens <- rbind(prec_nocens, c(n_samples, censoring, apply(params_nocens - true_mat, 2, sd)))
    }
  }
}

