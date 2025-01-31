library(tidyverse) 
library(data.table)
library(survival)
library(progress)
library(e1071)
library(splines)
source("Simulation_livers/Function_definition.R")

# Function generation settings: no treatment strategy
control <- list(
  dt = 0.01,      # unit
  dt_cs0 = 10,    # how often should the CSs be compared to the unit (untreated strategy)
  dt_cs1 = 1,     # how often should the CSs be compared to the unit (treated strategy)
  nt0 = 100*12,   # max num of small intervals before time of entry (and max num of CSs) 
  t_rmst = 3*100, # units for RMST
  rounding = 2
) 
control$nt <-  control$nt0 + control$t_rmst 

# Now rates refer to weeks
b <- 1
lambda <- 0.2 * control$dt
lambda_c <- 1/4 * control$dt
lambda_pt <- 0.1 * control$dt #0.05
beta <- 1.5
beta_c <- 2
beta_pt <- -2
gamma <- 1


pars <- list(
  var_tot = 1,       # total variance of X(t)
  tau = 0,           # standard deviation of measurement error
  rho = 0.9999,      # proportion of total variance represented by random person effect variance
  theta = 0.1,       # degree of mean reversal (to zero) of the Wiener process
  mu = 1,            # mean MELD
  delta = 0,         # average linear drift of MELD
  delta_sd =0.5,     # standard deviation of delta (individual drift of MELD)
  a = lambda,        # baseline Weibull rate WL survival
  b = 1,             # Weibull shape WL survival
  alph = 0.00,       # effect of age on hazard of event
  bet = beta,        # effect of MELD on hazard of event
  a_treat = lambda_c,# baseline Weibull rate treatment chances (time dep)
  b_treat = 1,       # Weibull shape treatment chances (time dep)
  alph_treat = 0,    # effect of age on hazard of treatment (time dep)
  bet_treat = beta_c,# effect of MELD on hazard of treatment (time dep)
  a_pt = lambda_pt,  # baseline Weibull rate treatment chances
  b_pt = 1,          # Weibull shape treatment chances,
  alph_pt = 0,       # effect of age on hazard of PT
  bet_pt = beta_pt,  # effect of MELD on hazard of PT
  gamm_pt = gamma    # effect of liver quality on hazard of PT
)

pars$omega <- sqrt(pars$rho * pars$var_tot)
pars$sigma <- sqrt((pars$var_tot - pars$omega^2) * 2 * pars$theta)

meta <- list(
  n = 2500,  # (we have 21398 IDs)
  M = 200,
  seed = 123
)

set.seed(meta$seed)
seeds <- round(runif(meta$M, 1000, 100000))

coeff_summ <- NULL
summ <- NULL
summ_cs <- NULL

for (i in 1:meta$M) {
  set.seed(seeds[i]) 
  # Data generation
  gd <- gen_data(control, pars, meta, monitor=FALSE, progress=TRUE)
  dfr <- gd$obs_data
  dfr$meld1 <- ave(dfr$meld, dfr$id, FUN = function (a) last(a))
  fumax <- control$nt * control$dt
  dfr <- dfr[dfr$tstart < fumax,] 
  ind <- (dfr$treatment == 0 & dfr$tstop > fumax)
  dfr$event[ind] <- 0
  dfr$tstop[ind] <- fumax
  dfr$active <- (dfr$meld > -1) # selection criteria
  # Transform data in the cross-section format: untreated strategy
  cs_df <- make_cs(dfr, tstart, tstop, t_entry, event, treatment, id, active, control$dt, control$dt_cs0, control$nt0)
  cs_df <- cs_df[!(cs_df$tstart != 0 & cs_df$treatment == 1), ] #remove lines where patients are treated after CSk
  # Transform data in the cross-section format: untreated strategy
  cs_df1 <- make_cs(dfr, tstart, tstop, t_entry, event, treatment, id, active, control$dt, control$dt_cs1, control$nt0)
  cs_df1 <- cs_df1[!(cs_df1$tstart != 0 & cs_df1$treatment == 1), ] #remove lines where patients are treated after CSk
  cs_df1$logtwait <- log(cs_df1$twait + 1)
  # Remove lines we no longer need from dfr:
  dfr <- dfr[!eq(dfr$tstart, dfr$tstop, control$rounding)| eq(0, dfr$tstop, control$rounding), ]
  # Add weights for the untreated strategy
  pred_vars_denom <- c("meld")
  pred_vars_num <- c("meld")
  pred_vars_w <- unique(c(pred_vars_denom, pred_vars_num))
  vars_pt <- c("meld1", "qliver") # PT-survival variables
  wdenom <- weights_denom0(dfr, id, pred_vars_denom, tstart, tstop, treatment, inclusion = active, control = control) 
  wnum <- weights_num0(cs_df, id, cs, pred_vars_num, tstop, treatment, control = control)
  w_df <- weights0(wdenom, wnum, cs_df, pred_vars_w, status, treatment, twait, 
                  vars_pt, control = control)
  w_df <- w_df[w_df$treatment==0,]
  # Add weights for the treated strategy
  cs_df1_0 <- cs_df1 |> group_by(id, cs) |> mutate(treated = last(treatment)) |> filter(treatment == 0) |> ungroup() 
  cs_df1_0$ntreat <- 1 # initialize covariate that counts how many people received treatment at that CS
  indxt <- (cs_df1_0$treated == 1)
  cs_df1_0$ntreat[indxt] <- ave(cs_df1_0$ntreat[indxt], cs_df1_0$cs[indxt], FUN = cumsum)
  qlivers <- cs_df1_0[indxt, c("id","cs","ntreat",vars_pt, "treated")]
  grid <- cs_df1_0 |> select(id, cs, ntreat) |>  group_by(cs) |>  expand(id, ntreat) |> ungroup()
  cs_df1_0 <- left_join(grid, select(cs_df1_0, -all_of(c(vars_pt, "treated", "ntreat"))), by = c("id", "cs")) |>
    left_join(select(qlivers, -c(id,treated)), by = c("cs", "ntreat")) |>
    left_join(select(qlivers, id, cs, ntreat, treated), by = c("id", "cs", "ntreat"))
  cs_df1_0$treated[is.na(cs_df1_0$treated)] <- 0
  fmla_w1 <- formula(paste("treated ~", paste(c("logtwait",pred_vars_denom), collapse = "+")))
  w_denom1 <-glm(fmla_w1, data=cs_df1_0, family=binomial(link="cloglog"))
  w1_df <- cs_df1[cs_df1$treatment==1,]
  pred.denom <- predict(w_denom1, newdata=w1_df, type="response")
  pred.num <- 1 - exp(-exp(w_denom1$coefficients["(Intercept)"]))
  w1_df$weights <- pred.num/pred.denom
  
  w1_df$event <- w1_df$status
  
  # Dataset for naive model
  dfr0 <- dfr
  dfr0$time <- ave(dfr0$tstop, dfr0$id, dfr$treatment, FUN = function (a) last(a))
  dfr0$event <- ave(dfr0$event, dfr0$id, dfr$treatment, FUN = function (a) last(a))
  dfr0 <- dfr0[dfr0$tstart==0 | dfr0$treatment == 1, ]
  dfr0$tstop <- dfr0$time
  dfr0$tstop[dfr0$treatment==1] <- (dfr0$tstop-dfr0$tstart)[dfr0$treatment==1] # Time reset
  dfr0$tstart[dfr0$treatment == 1] <- 0
  # Model estimation naive
  cox_naive0 <- coxph(Surv(tstop, event) ~ meld, ties = "breslow", data = dfr0[dfr0$treatment==0,])
  cox_naive1 <- coxph(Surv(tstop, event) ~ meld1 + qliver, ties = "breslow", 
                      data = dfr0[dfr0$treatment==1,])
  # Model estimation cross-sections
  # With weights
  cox_fit_w0 <- coxph(Surv(tstart, tstop, event) ~ meld,
                     ties = "breslow", data = w_df, weights = weights)
  cox_fit_w1 <- coxph(Surv(tstart, tstop, event) ~ meld1 + qliver,
                      ties = "breslow", data = w1_df, weights = weights)
  # Without weights
  cox_fit0 <- coxph(Surv(tstart, tstop, event) ~ meld, ties = "breslow", data = w_df)
  cox_fit1 <- coxph(Surv(tstart, tstop, event) ~ meld1 + qliver, ties = "breslow", data = w1_df)
  
  # RMSTs
  pred_df <- cs_df[cs_df$tstart == 0, ] %>%
    select(id, cs, meld, twait) %>%
    left_join(gd$counterfactual_df, by = c("id", "cs")) %>%
    left_join(gd$dfr_qliver, by = "cs")
  ## Naive approach (from baseline):
    meld_c <- cox_naive0$coefficients["meld"]
  meld_pt_c <- cox_naive1$coefficients[grep("meld1", names(cox_naive1$coefficients))]
  qliver_c <- cox_naive1$coefficients["qliver"]
  bh0 <- basehaz(cox_naive0, centered = FALSE)
  cumhaz0 <- c(0, bh0$hazard[bh0$time < control$t_rmst*control$dt])
  dt0 <- diff(c(0, bh0$time[bh0$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  bh1 <- basehaz(cox_naive1, centered = FALSE)
  cumhaz1 <- c(0, bh1$hazard[bh1$time < control$t_rmst*control$dt])
  dt1 <- diff(c(0, bh1$time[bh1$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  pred_df$esrmst0_naive <-  c(exp(-exp(pred_df$meld * meld_c)%*%t(cumhaz0)) %*%dt0)
  #lp1 <- rowSums(t(t(meldns_new)*meld_pt_c)) + pred_df$qliver * qliver_c
  lp1 <- pred_df$meld*meld_pt_c + pred_df$qliver * qliver_c
  pred_df$esrmst1_naive <- c(exp(-exp(lp1)%*%t(cumhaz1)) %*% dt1)

  ## With cross-sections
  # Weighted
  meld_c <- cox_fit_w0$coefficients["meld"]
  meld_pt_c <- cox_fit_w1$coefficients[grep("meld1", names(cox_fit_w1$coefficients))]
  qliver_c <- cox_fit_w1$coefficients["qliver"] 
  # Time reset
  bh0 <- basehaz(cox_fit_w0, centered = FALSE)
  bh1 <- basehaz(cox_fit_w1, centered = FALSE)
  cumhaz0 <- c(0, bh0$hazard[bh0$time < control$t_rmst*control$dt])
  dt0 <- diff(c(0, bh0$time[bh0$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  cumhaz1 <- c(0, bh1$hazard[bh1$time < control$t_rmst*control$dt])
  dt1 <- diff(c(0, bh1$time[bh1$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  lp0_w <- pred_df$meld * meld_c  
  pred_df$esrmst0_w <-  c(exp(-exp(lp0_w)%*%t(cumhaz0)) %*% dt0)
  lp1_w <- pred_df$meld*meld_pt_c + pred_df$qliver * qliver_c
  pred_df$esrmst1_w <- c(exp(-exp(lp1_w)%*%t(cumhaz1)) %*%dt1)
  
  # Unweighted
  meld_c <- cox_fit0$coefficients["meld"]
  meld_pt_c <- cox_fit1$coefficients[grep("meld1", names(cox_fit1$coefficients))]
  qliver_c <- cox_fit1$coefficients["qliver"]
  # Time reset
  bh0 <- basehaz(cox_fit0, centered = FALSE)
  bh1 <- basehaz(cox_fit1, centered = FALSE)
  cumhaz0 <- c(0, bh0$hazard[bh0$time < control$t_rmst*control$dt])
  dt0 <- diff(c(0, bh0$time[bh0$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  cumhaz1 <- c(0, bh1$hazard[bh1$time < control$t_rmst*control$dt])
  dt1 <- diff(c(0, bh1$time[bh1$time < control$t_rmst*control$dt], control$t_rmst*control$dt))
  lp0 <- pred_df$meld * meld_c
  pred_df$esrmst0_cs <-  c(exp(-exp(lp0)%*%t(cumhaz0)) %*% dt0)
  lp1 <- pred_df$meld*meld_pt_c + pred_df$qliver * qliver_c
  pred_df$esrmst1_cs <- c(exp(-exp(lp1)%*%t(cumhaz1)) %*%dt1)
  
  # Summary
  pred_df <- round(pred_df, 2)
  pred_df <- pred_df %>%
    select(rmst0, rmst1, esrmst0_naive, esrmst1_naive, esrmst0_cs, esrmst1_cs, 
           esrmst0_w, esrmst1_w, cs)
  bias_df <- pred_df
  bias_df$nrun <- i
  
  summ_cs <- rbind(summ_cs, bias_df)
  gc()
  
}

saveRDS(summ_cs, file = "Simulation_livers/results.rds")

