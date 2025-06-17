library(cmdstanr)
library(rtdists)

CENS <- FALSE
NOCENS <- TRUE

mod_norm <- cmdstan_model(file.path("StanModels/LBA.stan"))

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
      rt <- rlba_norm(n_samples, true_mat[i,1], 
		                 true_mat[i,2],  
		                 true_mat[i,3],  
		                 true_mat[i,4],  
				 sd_v, posdrift = TRUE)
      
      # Censor data
      censored <- rt[,"rt"] > quantile(rt[,"rt"], censoring)
      rt[censored,] <- matrix(rep(c(quantile(rt[,"rt"], censoring), 0), each=sum(censored)), sum(censored), 2)

      # Fit model that accounts for censoring
      if (CENS) { 
        fit_cens <- mod_cens$optimize(data=dat_cens, iter=1E6)
        params_cens[i,] <- c(fit_optim$summary(c("A", "b", "t0", "mean_v"))[,2])$estimate
      }
    
      # Fit model that ignores censoring
      if (NOCENS) { 
        RT=cbind(rt[rt[,"response"]==1,"rt"], 1)
        dat_nocens <- list(RT=RT, NUM_CHOICES=1, N=length(RT[,1]))
        fit_nocens <- mod_norm$optimize(data=dat_nocens, iter=1E6)
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

