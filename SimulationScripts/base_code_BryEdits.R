### Set Up

## Packages
library(cmdstanr)
library(rtdists)

## Toggle for to run censored or non-censored versions
CENS <- FALSE 
NOCENS <- TRUE

## Pull in Stan Model
mod_norm <- cmdstan_model(file.path("StanModels/LBA_censored.stan"))

## Toggle for plotting
PLOT <- FALSE

### Simulations

## Number of Iterations
n_iter <- 10

## True Parameters
A <- 1  #Start point
b <- 1.4  #Boundary
t0 <- .3  #Non-decision time
mean_v <- 1  #Drift rate mean
sd_v <- 1  #Drift rate sd

## Vector for True Parms
true_params <- c(A, b, t0, mean_v)

## Martrix for True Parms for each iteration
## 10 and 4???
true_mat <- matrix(rep(true_params, each=n_iter), nrow= n_iter , byrow = TRUE)

## Create matrices for bias and precision
if (CENS) { 
  bias_cens <- c()
  prec_cens <- c()
}
if (NOCENS) { 
  bias_nocens <- c()
  prec_nocens <- c()
}

##Outer Loops
## Loop for sample sizes
for (n_samples in c(50, 100, 200, 400)) { 
  
  ## Loop for censoring levels
  ## .7 and .95?
  for (censoring in seq(.7, .95, by=.05)) { 
    
    ## Matrices for each condition
    if (CENS) params_cens <- matrix(NA, n_iter, 4)
    if (NOCENS) params_nocens <- matrix(NA, n_iter, 4)
    
    
    ##Inner Loop: Simulate multiple datasets
    for (i in 1:n_iter) { 
      # Generate data directly from the model
      rt <- rlba_norm(n_samples, A=true_mat[i,1], 
                      b=true_mat[i,2],  
                      t0=true_mat[i,3],  
                      mean_v=true_mat[i,4],  
                      sd_v=sd_v, posdrift = TRUE)
      
      ## RAW STAN DATA
      stan_data_raw <- list(
        NUM_CHOICES = 1,  #number of choices
        N = dim(rt)[1],  #number of trials
        upper_cutoff = max(rt[,1]) + 1, #Upper cutoff (where it's not currently cutting anything off)
        lower_cutoff = 0,  #Lower cutoff
        RT = rt #RT is the variable in the list, rt is what was just generated
      )
      
      
      ## CENSORED MODEL 
      if (CENS) { 
        # Determine Cutoffs
        cens_upper_cutoff <- quantile(rt[,"rt"], censoring)
        cens_lower_cutoff <- 0
        
        # Censor data and replace with cutoffs (I think?)
        censored <- rt[,"rt"] > cens_upper_cutoff  #Flag/ remove RTs beyond cutoff
        rt_c <- rt
        rt_c[censored,] <- matrix(rep(c(quantile(rt[,"rt"], censoring), 0), each=sum(censored)), sum(censored), 2)
        
        #Create list of censored data for Stan
        stan_data_cens <- list(
          NUM_CHOICES = 1,
          N = dim(rt)[1],
          upper_cutoff = cens_upper_cutoff,
          lower_cutoff = cens_lower_cutoff,
          RT = rt_c
        )
        
        #Fit model using sampling not optim
        fit_cens <- mod_norm$sample(data=stan_data_cens, chains=4, parallel_chains=4)
        
        #Extract parameters
        params_raw <- fit_cens$summary(c("bMinusA", "A", "tau", "v_std[1]"))[,2]$mean 
        
        #Convert back to original parameter scale
        names(params_raw) <- c("bMinusA", "A", "tau", "v_std[1]")
        params <- rep(NA, 4); names(params) <- c("A", "b", "t0", "mean_v")
        params["A"] <- params_raw["A"]
        params["b"] <- params_raw["bMinusA"] + params_raw["A"]
        params["t0"] <- params_raw["tau"]
        params["mean_v"] <- 1+.1*params_raw["v_std[1]"]
        
        #Store censored parameter estimates
        params_cens[i,] <- params
      }
      
      
      ## NON-CENSORED MODEL
      if (NOCENS) { 
        # State Cutoffs (Not actually cutting off, but stan wants it)
        nocens_upper_cutoff <- max(rt[,1]) + 1
        nocens_lower_cutoff <- 0
        
        #Create list of noncensored data for Stan (matches raw)
        stan_data_nocens <- list(
          NUM_CHOICES = 1,
          N = dim(rt)[1],
          upper_cutoff = nocens_upper_cutoff,
          lower_cutoff = nocens_lower_cutoff,
          RT = rt
        )
        
        #Fit model using sampling not optim
        fit_nocens <- mod_norm$sample(data=stan_data_nocens, chains=4, parallel_chains=4)
        
        #Extract parameters
        params_raw <- fit_nocens$summary(c("bMinusA", "A", "tau", "v_std[1]"))[,2]$mean 
        
        #Convert back to original parameter scale
        names(params_raw) <- c("bMinusA", "A", "tau", "v_std[1]")
        params <- rep(NA, 4); names(params) <- c("A", "b", "t0", "mean_v")
        params["A"] <- params_raw["A"]
        params["b"] <- params_raw["bMinusA"] + params_raw["A"]
        params["t0"] <- params_raw["tau"]
        params["mean_v"] <- 1+.1*params_raw["v_std[1]"]
        
        #Store censored parameter estimates
        params_nocens[i,] <- params
      }
      
      
    ## Calculate and store bias and precision
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
}

