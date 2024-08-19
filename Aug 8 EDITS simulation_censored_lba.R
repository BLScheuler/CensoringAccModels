library(rtdists)
library(sft)


n_iter <- 10  #should be on par for estimating var

A <- 1
b <- 1.4
t0 <- .3
mean_v <- 1
sd_v <- 1
true_params <- c(A, b, t0, mean_v)

true_mat <- matrix(rep(true_params, each=n_iter), 10, 4)



bias <- c()
bias_n <- c()
prec <- c()
prec_n <- c()
n_samp_vec <- c(50, 100, 200, 400)

for (n_samples in n_samp_vec) { #diff size data sets
  for (censoring in seq(.7, .95, by=.05)) { #eg., 0.05 to 0.3
      
    all_params <- matrix(NA, n_iter, 4)
    all_params_n <- matrix(NA, n_iter, 4)
    for (i in 1:n_iter) { #inner loop: repeating sim study; estimate var due to data
      rt <- rlba_norm(n_samples, A, b, t0, mean_v, sd_v, posdrift = TRUE)
  #r random; lba norm < look up in help; n samples from lba model; posdrift requires df rates to be pos  
      censored <- rt[,"rt"] > quantile(rt[,"rt"], censoring) #pulling from middle loop in steps of .05; T/F based on percentile (< or > than eg 70th percentile)
      rt[censored,] <- matrix(rep(c(quantile(rt[,"rt"], censoring), 0), each=sum(censored)), sum(censored), 2)
      #replacing all censored rt with the max rt from censoring limit
      #change code to censoring early tail instead of late; change < to > as needed in line 31
     
      #Phase 2 section:
      # hData <- estimateNAH(RT=rt[,"rt"], CR=rt[,"response"]) #will change to K-M estimator 
      
      #stan_dat <- list(x=tvec, y=hData$H(tvec), N=length(tvec)) # need adapted for K-M
      
    #  fit_optim <- mod_cens$optimize(data=stan_dat, iter=1E6)
     # all_params[i,] <- c(fit_optim$summary(c("A", "b", "t0", "mean_v"))[,2])$estimate
      # Phase 2 end (for now)
    
      RT=rt[rt[,"response"]==1,] #cbind > column bind; subset where data is correct/ uncensored disappeared; #typical approach > how bad are things?
      stan_dat <- list(RT=RT, NUM_CHOICES=1, N=length(RT[,1])) #rt, num choice; num
      #need stan files for optimization component
      fit_optim2 <- mod_norm$optimize(data=stan_dat, iter=1E6)
      tmp_params <- c(fit_optim2$summary(c("A", "bMinusA", "tau", "v[1]"))[,2])$estimate
      all_params_n[i,] <- c(tmp_params[1], tmp_params[2] +tmp_params[1], tmp_params[3:4])
    }
#    bias <- rbind(bias, c(n_samples, censoring, apply(all_params - true_mat, 2, mean)))
#    prec <- rbind(prec, c(n_samples, censoring, apply(all_params - true_mat, 2, sd)))
    bias_n <- rbind(bias_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, mean)))
    prec_n <- rbind(prec_n, c(n_samples, censoring, apply(all_params_n - true_mat, 2, sd)))
  }
}

save(bias, prec, bias_n, prec_n, file="simstudy_20240509.Rdata")



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
