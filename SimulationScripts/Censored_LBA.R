library(rtdists)
library(sft)
library(cmdstanr)
library(here)


n_iter <- 10  #should be on par for estimating var

A <- 1   #accumulation starting point
b <- 1.4  #boundary threshold
t0 <- .3  #non-decision time
mean_v <- 1  #drift rate mean
sd_v <- 1  #drift rate sd
true_params <- c(A, b, t0, mean_v)  #create a vector for the true parameters

#Create a matrix for the true parameters
true_mat <- matrix(rep(true_params, each=n_iter), 10, 4)

# mod_cens <- cmdstan_model(file.path("../analysis/GP_LBA_diff2.stan"))
mod_norm <- cmdstan_model(here("SimulationScripts", "StanModels", "LBA.stan"))

#bias & precision 
bias_n <- c()
prec_n <- c()

#sample sizes that will be tested
n_samp_vec <- c(50, 100, 200, 400)

#Outer Loop: different sample sizes
for (n_samples in n_samp_vec) {  
  #Middle Loop: different censoring levels
  for (censoring in seq(.7, .95, by=.05)) { 
    ## NOTE TO BRY: try changing what will be censored! eg., 0.05 to 0.3

#Create matrix for parameter estimates 
    all_params_n <- matrix(NA, n_iter, 4)
    #Inner Loop: Repeating the actual simulation for each n_iter (estimate variance due to data)
    for (i in 1:n_iter) { 
      #Use normal distributions of LBA drift rates to generate RTs
      rt <- rlba_norm(n_samples, A, b, t0, mean_v, sd_v, posdrift = TRUE)
      #Note: posdrift requires df rates to be pos  
      
      #Pull the RTs that would be censored
      censored <- rt[,"rt"] > quantile(rt[,"rt"], censoring) #pulling from middle loop in steps of .05; T/F based on percentile (< or > than eg 70th percentile)
      
      #Replace all censored RTs with the max RT from censoring limit
      rt[censored,] <- matrix(rep(c(quantile(rt[,"rt"], censoring), 0), each=sum(censored)), sum(censored), 2)
      ##NOTE TO BRY: change code to censoring early tail instead of late; change < to > as needed in line 31
      
      #Create a subset where data is correct/ censored, and all the uncensored disappears (aka: how bad are things?)
      RT=rt[rt[,"response"]==1,] 
      
      #Pull the data we will fit with Stan model
      stan_dat <- list(RT=RT, NUM_CHOICES=1, N=length(RT[,1])) #rt, num choice; num

      #Stan Model (see GP_LBA file)
      fit_optim2 <- mod_norm$optimize(data=stan_dat, iter=1E6)
      
      #Pull parameter estimates from the Stan model
      tmp_params <- c(fit_optim2$summary(c("A", "bMinusA", "tau", "v[1]"))[,2])$estimate
      
      #Create a list  for the parameter estimates
      ##NOTE: why are we adding A & bMinusA?
      all_params_n[i,] <- c(tmp_params[1], tmp_params[2] +tmp_params[1], tmp_params[3:4])
    }
    
    #Bias: Estimated - True with means
    bias_n <- rbind(bias_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, mean)))
    
    #Precision: Estimated - True with sd
    prec_n <- rbind(prec_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, sd)))
  }
}

#Create a file for the results
save(bias, prec, bias_n, prec_n, file="simstudy_Pract1.Rdata")


##PLOT IT OUT


parmnames <- c("A", "b-A", "t0", "v")
parmmag <- c(.35, .45, .2, .25)
colvec <- viridis(4) #color profile; #4 will need changed later for n_samp_vec
par(mfrow=c(2,4))
for (parm in 1:4) { 
  plot(bias[bias[,1]==50,2], bias[bias[,1]==50,2+parm], type='l', col=colvec[1], lwd=2,
       ylim=parmmag[parm]*c(-1,1), main=parmnames[parm], ylab="Bias", xlab="Censoring Point")
  for (i in 2:4) { 
    b = n_samp_vec[i]
    lines(prec[bias[,1]==b,2], bias[bias[,1]==b,2+parm], col=colvec[i], lwd=2)
  }
  abline(h=0, col=grey(.7))
}

for (parm in 1:4) { 
  plot(bias_n[bias_n[,1]==50,2], bias_n[bias_n[,1]==50,2+parm], type='l', col=colvec[1], lwd=2,
       ylim=parmmag[parm]*c(-1,1), main=parmnames[parm], ylab="Bias", xlab="Censoring Point")
  for (i in 2:4) { 
    b = n_samp_vec[i]
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
