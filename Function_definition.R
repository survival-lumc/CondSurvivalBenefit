survweibull <- function(t, a, b) exp(- a * t^(b))
cumhazweibull <- function(t, a, b) a * t^(b)

eq <- function(x, y, rounding = 10) {
  x <- round(x, rounding)
  y <- round(y, rounding)
  return (x == y)
}

step_lin <- function(x, beta, val = 0.6) {
  if (x <= -val) y <- -val * beta
  if (x >= val) y <- val * beta
  if (x > -val & x < val) y <- x * beta
  return(y)
} 

### This function generates covariates (MELD, age), as well as survival for one 
### patient on the waiting list. 
gen_OU <- function(control, pars, urand, monitor=FALSE) {
  # Generate changes in the Wiener process 
  burnin <- 100
  rn <- diff(c(0, rwiener(end = control$nt + burnin, frequency = 1)))
  # Generate the Ornstein-Uhlenbeck process Xstar
  Xstar <- rep(NA, control$nt)
  change <- pars$sigma * rn[1] ### RMK: add 2 yrs burn-in
  Xstar[1] <- change
  for (i in 2:(control$nt +burnin)){
    change <- (-pars$theta * Xstar[i-1]) * control$dt + pars$sigma * rn[i]
    Xstar[i] <- Xstar[i-1] + change
  }
  Xstar <- Xstar[((burnin+1):(burnin+control$nt))]
  # Generate longitudinal MELD score from time 0 to the end of follow-up time
  bi <- round(rnorm(1, mean = pars$mu, sd = pars$omega), 2) #Baseline MELD
  di <- rnorm(1, mean = 0, sd = pars$delta_sd)*control$dt #Individual difference from MELD average growth
  tt <- 1:control$nt  #Times 
  X <- bi + (pars$delta + di) * tt + Xstar #Longitudinal MELD
  age <- round(rnorm(1, mean = 54, sd = 10), 2) #Generate patient's age
  # Patient's survival, here assumed to follow a Weibull distribution)
  hr <- exp(pars$bet * X + pars$alph * age) 
  # We define Haz := integral(hazard) with Haz(beginning interval) = 0
  Haz <- pars$a * hr * (tt^(pars$b) -(tt - 1)^(pars$b))
  Cumhaz <- cumsum(Haz) #Cumulative hazard
  Surv <- exp(-Cumhaz) #Survival
  # Draw event time via inversion formula for the correct time interval
  n0 <- min(sum(Surv >= urand), control$nt-1)
  # t0 <= tev < t1
  # In this interval: Surv = exp(-(Cumhaz(t0) + a*t^b*hr - a*t0^b*hr))
  # This formula fails at the two extremes n0 = 0 and n0 = max follow-up time
  uev <- ((-log(urand)-Cumhaz[n0]+pars$a * n0^pars$b * hr[n0+1])/(pars$a*hr[n0+1]))^(1/pars$b)
  if (n0 == 0) uev <- ((-log(urand))/(pars$a*hr[n0+1]))^(1/pars$b)
  if (uev < 10^-control$rounding) uev <- 10^-control$rounding
  tev <- uev * control$dt
  Surv <- c(1,exp(-Cumhaz))
  # #Upper intergral
  du <- 1/round(control$dt*10^control$rounding)
  hr_rmst <- rep(hr, each = round(1/du))
  tt_rmst <- rep(0:(control$nt-1), each = round(1/du))
  tt <- tt_rmst + rep(0:(round(1/du)-1)*du, rep = control$nt)
  Haz_rmst <- pars$a * hr * ((tt+du)^(pars$b) -tt^(pars$b))
  Surv_rmst <- exp(-cumsum(Haz_rmst))
  Surv_int <- ave(Surv_rmst, tt_rmst, FUN = function(a) sum(a)*du*control$dt)
  Surv_int <- Surv_int[diff(c(-1,tt_rmst))!=0]
  
  return(list(X=X, age=age, tev=tev, Cumhaz=Cumhaz, Surv = Surv, Surv_int = Surv_int, Xstar = Xstar))
}

# Generate time of treatment
# RMK: time of treatment should depend on what we observe, not the truth. Fix later.
gen_treat <- function(OU_obj, control, pars, urand, monitor=FALSE) { #qliver
  
  X <- OU_obj$X
  age <- OU_obj$age
  tt <- 1:control$nt
  # Draw event time
  hr <- exp(pars$bet_treat * X + pars$alph_treat*age) # pars$gamm_treat*qliver*(max_X - X)) 
  Haz <- pars$a_treat * hr * (tt^(pars$b_treat) - (tt - 1)^(pars$b_treat))
  Cumhaz <- cumsum(Haz)
  Surv <- exp(-Cumhaz)
  u_tr <- sum(Surv >= urand)
  if (u_tr == control$nt) {
    u_tr <- ((-log(urand)-Cumhaz[u_tr-1]+pars$a * u_tr^pars$b * hr[u_tr])/(pars$a*hr[u_tr]))^(1/pars$b)
    u_tr <- round(u_tr)
  }
  u_tr <- u_tr + sample(0:1, 1, prob = c(0.5, 0.5)) # to avoid always rounding down
  t_treat <- u_tr * control$dt #Assumption: a patient can only receive treatment when available
  return(t_treat)
}

## Generate PT survival for a certain time of treatment
gen_pt <- function(OU_obj, control, qliver, t_treat, pars, urand) { 
  
  X <- OU_obj$X
  Xstar <- OU_obj$Xstar
  u_treat <- max(min(round(t_treat/control$dt), control$nt), 1)
  # We base PT survival on the "real MELD" at time of treatment, age and quality of liver
  X_treat <- X[u_treat] - Xstar[u_treat]
  # Draw event time based on inversion method:
  # surv_pt <- exp(- a * t^(b) * i_risk) ---> x = (- ln(u)/(a*i_risk))^(1/b)
  hr <- exp(step_lin(X_treat, pars$bet_pt) + pars$gamm_pt * qliver) 
  u_pt <- (-log(urand)/(pars$a_pt * hr))^(1/pars$b_pt)
  t_pt <- u_pt * control$dt
  if (t_pt < 10^-control$rounding) t_pt <- 10^-control$rounding
  # Compute RMST until time horizon (control$t_rmst)
  du <- 1/round(control$dt*10^control$rounding)
  t_rmst <- control$t_rmst*control$dt
  Surv <- exp(-hr * pars$a_pt * (seq(0, control$t_rmst-du, du)) ^(pars$b_pt))
  rmst <- sum(Surv)*du*control$dt
  l_pt <- list(t_pt = t_pt, rmst = rmst)
  return(l_pt)
}

# Generate data with observed marker values
gen_data <- function(control, pars, meta, monitor=FALSE, progress=TRUE) {
  n <- meta$n #Number of subjects
  idx_train <- 1:n
  dfr <- NULL #Here we will store the synthetic data
  rmst0 <- rep(NA_real_, n*control$nt0) #Here we will store the counterfactual RMSTs
  rmst1 <- rep(NA_real_, n*control$nt0) #at all CSs for all patients
  t1 <- rep(NA_real_, n*control$nt0) #Counterfactual PT survival time at all CSs for all patients
  # Generate liver quality, one per CS (at each CS one liver comes in)
  qliver <- round(rnorm(control$nt, mean = 0, sd = 1), 2)
  if (progress) {
    cat("Generating training data ...\n")
    pb <- progress_bar$new(total = n)
  }
  for (i in idx_train) {
    if (progress) pb$tick()
    # Generate WL survival
    urand_wl <- runif(1)
    dd <- gen_OU(control, pars, urand_wl, monitor= FALSE)
    tev <- dd$tev
    uev <- tev / control$dt
    meld <- dd$X
    meld_obs <- meld + rnorm(control$nt, mean=0, sd=pars$tau)
    ni <- ceiling(uev) #number of observations for patient i
    # Generate entry time on the WL
    u_entry <- round(runif(1, min = 1, max = control$nt0))
    t_entry <- u_entry * control$dt
    # Generate treatment chances
    urand_c <- runif(1)
    t_treat <- gen_treat(dd, control, pars, urand_c) # RMK: maybe make this dependent on qliver later
    u_treat <- round(t_treat/control$dt)
    last_ob <- max(min(ni, u_treat, control$nt),1)
    last_t <- min(tev, t_treat) 
    # Treatment is only observed if it is administered before time of death
    # t_treat <- if_else(t_treat >= dd$tev, NA_real_, t_treat)
    
    # We generate all PT counterfactuals at all times we would have seen if treatments was not there.
    # We also generate the observed PT survival time in case treatment is administered
    # (if treatment is not given, PT survival time is unknown, hence NA)
    urand_pt <-  runif(1)
    pt_t <- NA_real_
    for (j in u_entry:(min(u_entry+ni-1, control$nt0))) {
      t_tr <- (j - u_entry +1) * control$dt
      pt <- gen_pt(dd, control, qliver[j], t_tr, pars, urand_pt)
      rmst1[control$nt0*(i-1) + j] <- pt$rmst #RMST1 for patient i at CS j
      t1[control$nt0*(i-1) + j] <- pt$t_pt #PT survival time for patient i at CS j
      rmst0j <- sum(dd$Surv_int[((j-u_entry+1):(j-u_entry+control$t_rmst))])/dd$Surv[j-u_entry+1]
      rmst0[control$nt0*(i-1) + j] <- rmst0j  #RMST0 for patient i at CS j
      if(j == u_entry+u_treat) pt_t <- pt$t_pt
    }
    qliveri <- qliver[u_entry + u_treat]
    
    if (u_treat < uev ) {
      meld_obs <- c(meld_obs[1:last_ob], meld_obs[last_ob])
      tstart <- c((0:(last_ob-1))*control$dt, t_treat)
      tstop <- c(if (last_ob != 1) (1:(last_ob-1))*control$dt, t_treat, pt_t + t_treat)
      treatment <- c(rep(0, last_ob), 1)
    }
    else {
      meld_obs <- meld_obs[1:last_ob]
      tstart <- c((0:(last_ob-1))*control$dt)
      tstop <- c(if (last_ob != 1) (1:(last_ob-1))*control$dt, tev)
      treatment <- rep(0, last_ob)
    }
    ni1 <- length(meld_obs)
    status <- rep(uev <= u_treat, ni1)
    dfri <- matrix(c(rep(i, ni1), rep(t_entry, ni1), tstart, tstop, status, 
                     treatment, meld_obs, round(dd$age + tstart, 2), 
                     rep(qliveri, ni1), rep(pt_t, ni1)), 
                   ni1, 10)
    dfr <- rbind(dfr, dfri)
  }
  dfr <- as.data.frame(dfr)
  names(dfr) <- c("id", "t_entry", "tstart", "tstop", "event", "treatment", 
                  "meld", "age", "qliver", "t_pt")
  dfr <- dfr[!is.na(dfr$tstop),] # some treated patients do not have t_pt (treated too late)
  
  dfr <- dfr %>% 
    mutate_at(vars(t_entry, tstart, tstop, t_pt), ~ round(.x, control$rounding)) 
  dfr <- dfr[!eq(dfr$tstart, dfr$tstop, control$rounding) | dfr$tstop == 0,]
  dfr$event <- ave(dfr$event,dfr$id,FUN= function(a) c(if (length(a)!=1) rep(0,length(a)-1),unique(a)))
  dfr$event[dfr$treatment == 1]  <- 1
  
  # liver quality at each CS
  dfr_qliver = data.frame(cs = 1:control$nt, qliver = qliver)
  
  counterfactual_df <- data.frame(
    id = rep(idx_train, each = control$nt0),
    cs = rep(1:control$nt0, rep = n),
    rmst0 = rmst0,
    rmst1 = rmst1,
    t_pt = t1
  )
  counterfactual_df <- counterfactual_df[!is.na(rmst1),] 
  
  return(list(obs_data=dfr, counterfactual_df = counterfactual_df, dfr_qliver = dfr_qliver))
}

### ----------------------------------------------------------------------------

make_cs <- function(data, tstart, tstop, tentry, status, treatment, ids, inclusion, 
                    dt, dt_cs, nt0, progress = TRUE) {
  
  if (missing(inclusion)) data$inclusion <- TRUE
  
  data <- data %>% 
    rename(id := {{ids}}, tentry := {{tentry}}, tstart := {{tstart}}, tstop := {{tstop}}, 
           status:= {{status}}, treatment := {{treatment}}, inclusion := {{inclusion}}) %>%
    arrange(id, tstop) %>%
    filter(inclusion) # only keep rows where patient is active 
  
  data$treated <- ave(data$treatment, data$id, FUN = function(a) last(a))
  data$event <- ave(data$status, data$id, FUN = function(a) last(a))
  data$last_tstop <-  ave(data$tstop, data$id, data$treatment, FUN = function(a) last(a))
  data$t_pt <- ave((data$tstop - data$tstart) * data$treatment, data$id, FUN = function(a) last(a))
  
  uentry <- round(data$tentry/dt)
  ustart <- round(data$tstart/dt)
  ustop <-  round(data$tstop/dt)
  csi <- seq(dt_cs, nt0, dt_cs)
  n <- length(csi)
  
  df_cs <- data.frame(matrix(NA_real_, ncol = length(data), nrow = 0))
  names(df_cs) <- names(data)
  df_cs <- add_column(df_cs, cs = NA_real_, .after = "id")
  
  if (progress) {
    cat("Generating cross-sections ...\n")
    pb <- progress_bar$new(total = n)
  }
  
  for (i in 1:n) {
    if (progress) pb$tick()
    ind <- (uentry + ustart <= csi[i] & uentry + ustop > csi[i] & data$treatment == 0)  # patient in the CS ind
    ind1 <- (uentry + ustart == csi[i] & data$treatment == 1) # patients for whom I want to take the previous line setting tstart = tstop
    df_csi <- data[ind, ] 
    df_csi1 <- data[c(ind1[-1], FALSE),]
    df_csi1$tstart <- df_csi1$tstop
    df_csi <- add_column(rbind(df_csi, df_csi1), cs = csi[i], .after = "id")
    df_cs <- rbind(df_cs, df_csi)
  }
  
  df_cs$tstart <- 0
  uentry_cs <- round(df_cs$tentry/dt)
  utimes_i <- df_cs$last_tstop/dt # Time since entry to either treatment or death in weeks
  uwait <- (df_cs$cs - uentry_cs)
  df_cs$tstop <- (round(utimes_i - uwait))*dt # Time since CS to treatment or death in yrs
  df_cs <- add_column(df_cs, twait = round(uwait*dt, 3), .after = "cs")
  df_cs$status <- (1-df_cs$treated)*df_cs$event
  df_cs_pt <- df_cs[df_cs$treated == 1, ]
  df_cs_pt$tstart <- df_cs_pt$tstop # Time to PT-death since CS
  df_cs_pt$tstop <- df_cs_pt$t_pt + df_cs_pt$tstop # Time to PT-death since CS
  df_cs_pt$treatment <- 1
  df_cs_pt$status <- df_cs_pt$treated*df_cs_pt$event
  
  df_cs <- rbind(df_cs, df_cs_pt)
  df_cs <- select(df_cs, -c(tentry, inclusion, treated, event, last_tstop, t_pt))
  # Order the dataframe correctly
  df_cs <- df_cs %>% arrange(cs,id,treatment)
  
  return(df_cs)
}


weights_denom0 <- function(obs_data, id, pred_vars, tstart, tstop, treatment, 
                          strat, inclusion, control) {
  
  dt <- 10^(-control$rounding)
  
  if (!missing(strat)) cox_strata <- factor(eval(substitute(strat), obs_data))
  
  if (missing(inclusion)) {
    obs_data <- obs_data %>% 
      select(id = {{id}}, tstart = {{tstart}}, tstop = {{tstop}}, treatment = {{treatment}},
             all_of(pred_vars)) %>%
      mutate(inclusion = TRUE)
  }
  
  if (!missing(inclusion)) {
    obs_data <- obs_data %>% 
      select(id = {{id}}, tstart = {{tstart}}, tstop = {{tstop}}, treatment = {{treatment}},
             all_of(pred_vars), inclusion = {{inclusion}}) 
  }
  
  if (!missing(strat)) obs_data$strata <- cox_strata
  
  obs_data$tstart <- round(obs_data$tstart/dt)
  obs_data$tstop <- round(obs_data$tstop/dt)
  cox_df <- obs_data[obs_data$inclusion,]
  cox_df$treated <- ave(cox_df$treatment, cox_df$id, FUN = function(a) last(a))
  cox_df <- cox_df[cox_df$treatment == 0,]
  cox_df$status_treat <- ave(cox_df$treated, cox_df$id, 
                             FUN = function(a) c(rep(0, length(a)-1), last(a)))
  cox_df <- select(cox_df, -c(treated))
  cox_df$tstop[cox_df$tstop==0] <- 10^{-3}
  cox_df <- cox_df[cox_df$tstart != cox_df$tstop, ]
  
  # Fit the cox model
  if (missing(strat)) fmla <- formula(paste("Surv(tstart, tstop, status_treat) ~ ",
                                            paste(pred_vars, collapse = "+" )))
  if (!missing(strat)) fmla <- formula(paste("Surv(tstart, tstop, status_treat) ~ ",
                                             paste(pred_vars, collapse = "+" ),
                                             " + strata(strata)"))
  
  cox_fit <- coxph(fmla, data = cox_df)
  
  # Choosing relevant time points from the baseline cumulative hazard
  
  if (missing(strat)) {
    bh <-basehaz(cox_fit, centered = TRUE)
    bh$cumhaz <- round(bh$hazard, 2*control$rounding)
    bh$haz <- diff(c(0, bh$cumhaz))
    bh <- bh[bh$haz != 0,]
    bh_times <- bh$time
  }
  
  if (!missing(strat)) {
    bh <-basehaz(cox_fit, centered = TRUE) %>% arrange(strata, time)
    bh$cumhaz <- round(bh$hazard, 2*control$rounding)
    bh$haz <- ave(bh$cumhaz, bh$strata, FUN = function(a) diff(c(0,a)))
    bh <- bh[bh$haz != 0,]
    bh_times <- bh %>% select(strata, time)
  }
  
  # Merge the time points at bh_times to the time points where X_obs changes per strata
  obs_data <- obs_data[obs_data$treatment == 0,]
  obs_data$tstop[obs_data$tstop==0] <- 10^(-2*control$rounding)
  obs_data$status_treat <- 0
  if (missing(strat)) {
    w_df <- survSplit(obs_data, cut = bh_times, start = "tstart",
                      end = "tstop", event = "status_treat")
  }
  
  if (!missing(strat)) {
    w_df <- data.frame()
    for (i in 1:length(levels(obs_data$strata))) {
      w_dfi <- survSplit(obs_data[obs_data$strata == levels(obs_data$strata)[i], ],
                         cut = bh_times$time[bh_times$strata == levels(obs_data$strata)[i]],
                         start = "tstart", end = "tstop", event = "status_treat")
      w_df <- rbind(w_df, w_dfi)
    }
  }
  # Predicting survival
  w_df$haz <- predict(cox_fit, newdata = w_df, type = "expected")
  w_df$dweight <- ave(w_df$haz, w_df$id, FUN = function(a) exp(-cumsum(a)))
  w_df <- select(w_df, -c(status_treat))
  w_df$tstop <- round(w_df$tstop)*dt
  w_df$tstart <- w_df$tstart*dt
  
  return(list(dweight_df = w_df, cox_denom = cox_fit, cox_df = cox_df))
}

weights_num0 <- function(cs_data, id, cs, pred_vars, tstop, treatment, strat, control) {
  
  dt <- 10^(-control$rounding)
  
  if (!missing(strat)) cox_strata <- eval(substitute(strat), cs_data)
  
  # Data for the cox model
  cs_data <- cs_data %>% 
    select(id = {{id}}, cs = {{cs}}, tstop = {{tstop}}, treatment = {{treatment}},
           all_of(pred_vars), strata = {{cs}}) %>%
    mutate(strata = paste("cs", cs, sep =""))
  
  if (!missing(strat)) cs_data$strata <- paste(cs_data$strata, cox_strata, sep = "_")
  
  cs_data$tstop <- round(cs_data$tstop/dt)
  cs_data$status_treat <- ave(cs_data$treatment, cs_data$id, FUN = function(a) last(a))
  cs_data <- cs_data[cs_data$treatment == 0,]
  cs_data <- select(cs_data, -c(treatment))
  
  cox_df <- cs_data#[cs_data$tstop != 0,] #REMOVE!!
  
  # Fit the cox model
  fmla <- formula(paste("Surv(tstop, status_treat) ~ ", 
                        paste(pred_vars, collapse = "+"),
                        " + strata(strata) + cluster(id)"))
  cox_fit <- coxph(fmla, data = cox_df)
  
  # Choosing relevant time points from the baseline cumulative hazard
  bh <- basehaz(cox_fit, centered = TRUE) 
  if (is.null(bh$strata)) bh$strata <- unique(cs_data$strata)
  bh <- bh %>%  arrange(strata, time)
  bh$cumhaz <- round(bh$hazard, 2*control$rounding)
  bh$haz <- ave(bh$cumhaz, bh$strata, FUN = function(a) diff(c(0,a)))
  bh <- bh[bh$haz != 0 & bh$time != 0,]
  bh_times <- bh %>% select(strata, time) 
  bh_times$time <- round(bh_times$time, control$rounding)
  bh_times <- unique(bh_times)
  bh_times$tstart <- ave(bh_times$time, bh_times$strata, FUN = function(a) c(0, a[-length(a)]))
  
  # Merge the time points at bh_times to the time points where X_obs changes
  w_df <- left_join(cs_data, bh_times, by = "strata", relationship = "many-to-many") %>%  arrange(id, cs, time) 
  w_df$tstart[is.na(w_df$time)] <- 0
  w_df$time[is.na(w_df$time)] <- max(bh_times$time ) + 1
  w_df <- w_df[w_df$tstart <= w_df$tstop,]
  w_df <- w_df[w_df$tstart < w_df$tstop | w_df$tstop == 0,]
  
  indx <- (diff(c(w_df$id, 0))==0 & diff(c(w_df$cs, 0))==0)
  w_df$tstop[indx] <- w_df$time[indx]
  w_df$tstop <- round(w_df$tstop, control$rounding)
  
  # Predicting survival
  w_df$cumhaz <- predict(cox_fit, newdata = w_df, type = "expected")
  w_df$nweight = exp(-w_df$cumhaz)
  
  w_df <- w_df %>% select(-c(time)) %>% rename(treated = status_treat)
  w_df$tstop <- w_df$tstop*dt
  w_df$tstart <- w_df$tstart*dt
  
  return(list(nweight_df = w_df, cox_num = cox_fit, cox_df = cox_df))
}

weights0 <- function(weights_denom_obj, weights_num_obj, cs_data, pred_vars, event, 
                    treatment, twait, vars_at_treat, control, truncation, cens) {
  
  dt <- 10^(-control$rounding)
  wdenom <- weights_denom_obj
  wnum <- weights_num_obj
  nw_df <- wnum$nweight_df |> select(-c(cumhaz, nweight)) 
  nw_df$tstart <- round(nw_df$tstart/dt)
  nw_df$tstop <- round(nw_df$tstop/dt)
  dw_df <- wdenom$dweight_df |> select(-c(haz, dweight))
  dw_df$tstart <- round(dw_df$tstart/dt)
  dw_df$tstop <- round(dw_df$tstop/dt)
  cox_num <- wnum$cox_num
  cox_denom <- wdenom$cox_denom
  if(missing(vars_at_treat)) vars_at_treat <- NULL
  #dfr <- cs_df |> rename(event = status) |>mutate_at(vars( tstart, tstop, twait), ~ round(.x/dt))
  dfr <- cs_data |> # post-treatment info
    select(id, cs, event = {{event}}, treatment = {{treatment}}, tstart, tstop,
           twait = {{twait}}, all_of(vars_at_treat)) |>
   mutate_at(vars( tstart, tstop, twait), ~ round(.x/dt))
  dfr$treated <- ave(dfr$treatment, dfr$id, dfr$cs, FUN = function(a) last(a))
  dfr$event <- ave(dfr$event, dfr$id, dfr$cs, FUN = function(a) last(a))
  dfr$t_pt <- ave(dfr$tstop - dfr$tstart, dfr$id, dfr$cs, FUN = function(a) last(a)) 
  dfr <- select(dfr, - c(tstart, tstop))
  eps <- 10^{-2}
  
  # Df with all the patients IDs and the CSs they cross during their time on WL
  id_cs <- unique(dfr[dfr$treatment == 0, c("id", "cs", "twait")]) 
  id_cs <- left_join(unique(nw_df[, c("id", "cs")]), id_cs, by = c("id", "cs")) 
  
  # Adding CS to the df dweights (which is in observational time)
  dw_df <- left_join(dw_df, id_cs, by = "id", relationship = "many-to-many") |> arrange(id, cs, tstop)
  # Removing patients that stayed on the WL too shortly and were therefore not included in any CS
  dw_df <- dw_df[!is.na(dw_df$cs),]
  # Reset tstop for dweight, so that it has the same time scale as nweight
  dw_df$tstop <- dw_df$tstop - dw_df$twait
  dw_df$last_tstop <- ave(dw_df$tstop, dw_df$id, dw_df$cs, FUN = function(a) last(a))
  dw_df <- dw_df[dw_df$tstop >= 0,]
  dw_df <- dw_df[dw_df$tstop > 0 | dw_df$last_tstop==0,] #Include patients with txp at cs
  
  dw_df <- merge(dw_df, nw_df[, c("id", "cs","tstop")], by = c("id", "cs", "tstop"), all = TRUE) 
  dw_df <- dw_df |>   # fill NAs and add tstart:
    arrange(id, cs, tstop) |> 
    group_by(id, cs) |> 
    fill(-c(id, cs, tstop, tstart), .direction = "downup") |>
    fill(inclusion, .direction = "down") |>
    mutate(tstart = c(0, tstop[-n()])) |>
    ungroup() 
  
  dw_df <- dw_df[dw_df$tstart != dw_df$tstop | dw_df$tstop == 0,] 
  
  nw_df <- nw_df[diff(c(0,nw_df$id))!=0 | diff(c(0,nw_df$cs))!=0,] |>
    select(-c(tstart, tstop)) |>
    full_join(dw_df[, c("id","cs","tstart","tstop","inclusion")], by = c("id","cs")) |>
    relocate(tstart, tstop, inclusion, .after = "cs")
  
  # For prediction purposes, dweight needs to go back to follow-up time
  dw_df$tstop <- dw_df$tstop + dw_df$twait
  dw_df$tstart <- dw_df$tstart + dw_df$twait
  
  if (!is.null(dw_df$strata)) {
    keep <- (!is.na(dw_df$strata))
    dw_df <- dw_df[keep, ]
    nw_df <- nw_df[keep, ]
  }
  
  # Weights
  cox_df <- wdenom$cox_df
  indx <- (dw_df$tstart == dw_df$tstop)
  dw_df$tstart <- dw_df$tstart - eps
  dw_df$tstop <- dw_df$tstop - eps
  dw_df$tstop[indx] <- dw_df$tstop[indx] + 2*eps
  # Probability of staying treatment-free in a small interval of time (denom)
  dw_df$status_treat <- 0
  hazdenom <- predict(cox_denom, newdata = dw_df, type = "expected")
  # Probability of receiving treatment in a small interval of time (denom)
  pdenom0 <- exp(-hazdenom)
  # Probability of staying treatment-free in a small interval of time (num)
  cox_df <- wnum$cox_df
  nw_df$status_treat <- 0
  nw_df$status_treat <- ave(nw_df$treated, nw_df$id, nw_df$cs,
                            FUN = function(a) c(if (length(a)>1) rep(0, length(a)-1), last(a)))
  nw_df <- select(nw_df, -c(treated))
  cumhaznum <- predict(cox_num, newdata = nw_df, type = "expected")
  haznum <- ave(cumhaznum, nw_df$id, nw_df$cs, FUN = function(a) diff(c(0, a)))
  # Probability of receiving treatment in a small interval of time (num)
  pnum0 <- exp(-haznum)
  # # Denominator probabilities
  nw_df$pdenom <- if_else(nw_df$status_treat==0, pdenom0, 1 - pdenom0)
  nw_df$pnum <- if_else(nw_df$status_treat==0, pnum0, 1 - pnum0)
  rm(dw_df, cox_num, cox_denom, pdenom0, pnum0, hazdenom)
  nw_df$tstop <- round(nw_df$tstop)
  
  # Weights
  nw_df <- nw_df |> filter(inclusion) |> select(-c(status_treat, inclusion))
  nw_df$wdenom <- ave(nw_df$pdenom, nw_df$id, nw_df$cs, FUN = function(a) cumprod(a))
  nw_df$wnum <- ave(nw_df$pnum, nw_df$id, nw_df$cs, FUN = function(a) cumprod(a))
  nw_df$weights <- nw_df$wnum/nw_df$wdenom
  nw_df <- select(nw_df, -c(pdenom, pnum))
  
  # Add weights to nw_df
  nw_df <- left_join(nw_df, dfr[dfr$treatment==0, ], by = c("id","cs"), relationship = "many-to-many")
  nw_df2 <- nw_df[(diff(c(nw_df$id,0)) != 0 | diff(c(nw_df$cs,0)) != 0) & nw_df$treated == 1,]
  nw_df$treatment <- 0
  nw_df2$treatment <- 1
  nw_df2$tstart <- nw_df2$tstop
  nw_df2$tstop <- nw_df2$tstop + nw_df2$t_pt
  nw_df <- rbind(nw_df, nw_df2) |> select(-c(t_pt, treated))# Add one line to nw_df for treated patients
  nw_df$event <- ave(nw_df$event, nw_df$id, nw_df$cs, 
                     FUN = function(a) c(if (length(a)>1) rep(0, length(a)-1), last(a)))
  nw_df$wdenom <- ave(nw_df$wdenom, nw_df$id, nw_df$cs, FUN = function(a) c(1, a[-length(a)]))
  nw_df$wnum <- ave(nw_df$wnum, nw_df$id, nw_df$cs, FUN = function(a) c(1, a[-length(a)]))
  nw_df$weights <- ave(nw_df$weights, nw_df$id, nw_df$cs, FUN = function(a) c(1, a[-length(a)]))
  
  nw_df <- nw_df[nw_df$tstop > eps,] # Remove lines with tstop == 0
  nw_df$tstop <- nw_df$tstop*dt
  nw_df$tstart <-nw_df$tstart*dt
  nw_df$twait <-nw_df$twait*dt
  
  return(nw_df)
}

pred_rmst <- function(coxmodel, newdata, treat_var, control, progress = TRUE, treat_vars=NULL) { 
  
  cox_terms <- attr(coxmodel$terms, "variables")
  strat <- grep("strata", cox_terms)
  if (length(strat) == 0) strat <- "nostrata"
  cox_strata <- as.character(cox_terms[[strat]])
  cox_strata <- cox_strata[-which(cox_strata =="strata")]
  n <- nrow(newdata)
  rmst <- rep(NA_real_, n)
  newdata <- newdata %>% mutate_at(vars(eval(cox_strata)), ~ as.character(.))
  newdata$treat_var <- eval(substitute(treat_var), newdata)
  
  ### Baseline hazard
  basehazard <- data.table(basehaz(coxmodel, centered = FALSE))
  basehazard <- basehazard[time < control$t_rmst*control$dt]
  names(basehazard)[which(names(basehazard)=="strata")] <- cox_strata
  ### Risk 
  newdata$pred_risk <- predict(coxmodel, newdata = newdata, type = "risk", reference = "zero")
  newdata <- newdata %>% select("id", eval(cox_strata), "pred_risk", treat_var, all_of(treat_vars))
  newdata <- data.table(newdata)
  
  ### Survival per patient
  rmst <- rep(NA_real_, nrow(newdata))
  if (length(cox_strata) !=0) {
    strata_levels <- unique(newdata[[cox_strata]])
    num_levels <- length(strata_levels)
    for (j in 1:num_levels) {
      basehazard_j <- basehazard[basehazard[[cox_strata]] == strata_levels[j]]
      times_j <- diff(c(0, basehazard_j$time, control$t_rmst*control$dt))
      which_j <- which(newdata[[cox_strata]] == strata_levels[j])
      newdata_j <- newdata[newdata[[cox_strata]] == strata_levels[j]]
      if (nrow(newdata_j) != 0) {
        for (i in 1:nrow(newdata_j)) {
          surv_ij <- c(1, exp(-newdata_j[i]$pred_risk*basehazard_j$hazard))
          if (!missing(treat_vars)) {
            hr_fixed_ij <-  newdata_j[i]$pred_risk
            hr_timedep_ij <- newdata_j[i, ..treat_vars]*coxmodel$coefficients[treat_vars]
            hr_timedep_ij <- t(t(as.vector(hr_timedep_ij))%*%c(log(times_j[-1]))) 
            risk_ij <- hr_fixed_ij + rowSums(hr_timedep_ij)
            haz <- diff(c(0, basehazard_j$hazard))
            surv_ij <- c(1, exp(-cumsum(risk_ij*haz)))
          }
          rmst[which_j[i]] <- sum(surv_ij*times_j)
        }
      }
    }
  }
  else {
    times <- diff(c(0, basehazard$time, control$t_rmst*control$dt))
    for (i in 1:nrow(newdata)) {
      surv_i <- c(1, exp(-newdata[i]$pred_risk*basehazard$hazard))
      rmst[i] <- sum(surv_i*times)
    }
  }
  
  return(list(rmst))
}


