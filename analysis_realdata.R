library(tidyverse)
library(survival)
library(progress)
library(janitor)
library(latex2exp)
source(here::here("Simulation","Function_definition_new.R"))

# Function generation settings: no treatment strategy
control <- list(
  dt = 0.01,     # unit
  dt_cs = 10,    # how often should the CSs be compared to the unit
  nt0 = 100*12,  # max num of small intervals before time of entry (and max num of CSs) 
  t_rmst = 3*100, # units for RMST
  rounding = 2
) 
control$nt <-  control$nt0 + control$t_rmst 

# Function generation settings: treated strategy
control1 <- list(
  dt = 0.01,     # unit
  dt_cs = 1,    # how often should the CSs be compared to the unit
  nt0 = 100*12,  # max num of small intervals before time of entry (and max num of CSs) 
  t_rmst = 3*100, # units for RMST
  rounding = 2
) 
control1$nt <-  control1$nt0 + control1$t_rmst 


### 0a) Load data --------------------------------------------------------------
### ----------------------------------------------------------------------------
load(here::here("survival", "longdata.rda")) #load longitudinal data

### 0b) Add some variables -----------------------------------------------------
### ----------------------------------------------------------------------------
longdata$recipient_age <- round(longdata$recipient_age + longdata$tstart)
longdata$log_bili <- log(longdata$bilirubine_mg_dl +1)
longdata$log_creat <- log(longdata$creatinine_mg_dl +1)
longdata$log_inr <- log(longdata$inr +1)

longdata <- select(longdata, - c(bilirubine_mg_dl, creatinine_mg_dl, inr))

# Expand data ------------------------------------------------------------------
# Variables that we center and standardize
cont_vars <- c("recipient_age", "r_weight", "r_height", "meld_nat_match", "meld_lab",
               "log_bili", "log_inr", "log_creat", "graft_age", "graft_et_dri") # Maybe add cumulative time        
# Logical variable that I want to turn into dummies
lgl_vars <- c("patient_male", "dummy_biweekly_dialysis", "txp_rescue", "graft_dcd")  
fact_vars <- c("recipient_bloodgroup", "disease_group", "donor_death_cause_group",
               "txp_destination")                      # expand

covariates <- fact_vars
longdata2 <- model.frame(~ ., longdata, na.action=na.pass)
for (j in 1:length(covariates)) {
  vars <- longdata[, covariates[j]][[1]]
  lev <- levels(vars)
  dummies <- model.matrix(formula(paste("~", covariates[j], "-1")), 
                          na.action = "na.pass", data = longdata2)
  dummies <- data.frame(dummies)
  dummies <- dummies[,-1, drop = FALSE]
  names(dummies) <- paste(covariates[j], lev[-1], sep = "_")
  names(dummies) <- paste(covariates[j], lev[-1], sep = "")
  longdata <- cbind(select(longdata, -c(covariates[j])), dummies)
}
rm(dummies)

longdata$meld_int <- longdata$meld_lab
longdata <- longdata |>
  mutate_at(vars(all_of(lgl_vars)), ~ if_else(.x, 1, 0)) |>
  mutate_at(vars(all_of(cont_vars)), ~ (.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)) |>
  clean_names()
longdata$tentry <- round(longdata$tentry,2)
rm(longdata2)

# There are a few duplicate lines for some reason, keep only one line
# (these are patients who had multiple lines in the very first day, we take the
# last line as correct)
dupl <- longdata[, c("id", "tstart", "txp")] |> group_by(id,tstart,txp) |> mutate(dupl = n():1) |> pull(dupl)
longdata <- longdata[dupl == 1,]

### 1) Make cross-sections -----------------------------------------------------
### ----------------------------------------------------------------------------
# For the untreated strategy
cs_df <- make_cs(longdata, tstart, tstop, tentry, event, txp, id, active, control)
cs_df <- cs_df[!(cs_df$tstart != 0 & cs_df$treatment == 1), ] #remove lines where patients are not treated after CSk
cs_df$logtwait <- log(cs_df$twait + 1)
# For the treated strategy
cs_df1 <- make_cs(long, tstart, tstop, tentry, event, txp, id, active, control1)
cs_df1 <- cs_df1[!(cs_df1$tstart != 0 & cs_df1$treatment == 1), ] #remove lines where patients are not treated after CSk
cs_df1$logtwait <- log(cs_df1$twait + 1)

### 2) Add Weights -------------------------------------------------------------
### ----------------------------------------------------------------------------
pred_vars <- c("recipient_age", "patient_male", "r_weight", "r_height", "meld_nat_match", 
               "dummy_biweekly_dialysis", "recipient_bloodgroup_ab", 
               "recipient_bloodgroup_b", "recipient_bloodgroup_o")
pt_vars <- c("txp_rescue", "graft_age", "graft_dcd", "disease_group_metabolic", 
             "donor_death_cause_group_cva_stroke", "donor_death_cause_group_trauma", 
             "donor_death_cause_group_anoxia", "txp_destination_h", "txp_destination_l") 
wl_vars <- c("log_bili", "log_inr", "log_creat", "meld_lab")
inter_vars <- c("graft_age:meld_nat_match", "txp_rescue:meld_nat_match",
                "recipient_age:graft_age", "r_weight:graft_age", "r_height:graft_age")
# Weights for treated strategy
cs_df1_0 <- cs_df1 |> group_by(id, cs) |> mutate(treated = last(treatment)) |> filter(treatment == 0) |> ungroup() 
cs_df1_0$ntreat <- 1 # initialize covariate that counts how many people received treatment at that CS
indxt <- (cs_df1_0$treated == 1)
cs_df1_0$ntreat[indxt] <- ave(cs_df1_0$ntreat[indxt], cs_df1_0$cs[indxt], FUN = cumsum)
qlivers <- cs_df1_0[indxt, c("id","cs","ntreat",pt_vars, "treated")]
grid <- cs_df1_0 |> select(id, cs, ntreat) |>  group_by(cs) |>  expand(id, ntreat) |> ungroup()
cs_df1_0 <- left_join(grid, select(cs_df1_0, -all_of(c(pt_vars, "treated", "ntreat"))), by = c("id", "cs")) |>
  left_join(select(qlivers, -c(id,treated)), by = c("cs", "ntreat")) |>
  left_join(select(qlivers, id, cs, ntreat, treated), by = c("id", "cs", "ntreat"))
cs_df1_0$treated[is.na(cs_df1_0$treated)] <- 0
fmla_w1 <- formula(paste("treated ~", paste(c("logtwait",pred_vars, pt_vars, inter_vars), collapse = "+")))
w_denom1 <-glm(fmla_w1, data=cs_df1_0, family=binomial(link="cloglog"))
# Note that the model removes some observations containing NAS: these are the
# cross-sections without any transplantations, which are correctly removed
w1_df <- cs_df1[cs_df1$treatment==1,]
pred.denom <- 1-predict(w_denom1, newdata=w1_df, type="response")
#num <- 1 - exp(-exp(w_denom1$coefficients["(Intercept)"]))
w1_df$weights <- 1/pred.denom
rm(cs_df1_0, cs_df1, grid)

# Weights for the untreated strategy
longdata <- longdata[!eq(longdata$tstart, longdata$tstop, control$rounding)| 
                       eq(0, longdata$tstop, control$rounding), ]
w_denom0 <- weights_denom0(longdata, id, pred_vars, tstart, tstop, txp, 
                         rec_country, active, control) 
# add twait
w_num0 <- weights_num0(cs_df, id, cs, pred_vars, tstop, treatment, rec_country, control)
gc()
w_df <- weights0(w_denom0, w_num0, cs_df, pred_vars, status, treatment, twait, 
                  c(pt_vars, wl_vars), control)
w_df <- w_df[w_df$tstart < control$t_rmst, ] # remove lines after the end of follow-up
w_df <- w_df[w_df$treatment == 0,] # remove post-treatment lines
summary(w_df$weights)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.2231   1.0276   1.0767   1.1779   1.1808 944.5366 
upq <- quantile(w_df$weights, probs=0.9995)   
w_df <- w_df %>% # Truncate weights at 0.0001
  mutate(weights = if_else(weights > upq, upq, weights))
hist(w_df$weights, breaks = 1000, xlim = c(0.75,1.25))

w_df$log_twait <- log(w_df$twait+1)
#w_df$tstop <- if_else(w_df$tstop == 0 & w_df$treatment == 0, 10^-4, w_df$tstop)
save(w_df, w1_df, file = here::here("survival","rda_output","weights_df1.rda"))
save(w_num0, w_denom0, w_denom1, file = here::here("survival","rda_output","txp_models.rda"))
summ_num <- summary(w_num0$cox_num)
summ_denom <- summary(w_denom0$cox_denom)
summ_denom1 <- summary(w_denom1)
save(summ_num, summ_denom, summ_denom1, file = here::here("survival","rda_output","txp_summ.rda"))
rm(w_num, w_denom)
gc()

### 3) Model -------------------------------------------------------------------
### ----------------------------------------------------------------------------
cox_wl_vars <- c("recipient_age", "patient_male", "r_weight", "r_height", 
                "log_bili", "log_inr", "log_creat", "dummy_biweekly_dialysis", 
                "recipient_bloodgroup_ab", "recipient_bloodgroup_b", 
                "recipient_bloodgroup_o", "log_twait")
donor_vars <- c("graft_age", "graft_dcd",  "donor_death_cause_group_cva_stroke", 
                "donor_death_cause_group_trauma", "donor_death_cause_group_anoxia", 
                "txp_rescue", "txp_destination_h", "txp_destination_l")
cox_pt_vars <- c("recipient_age", "patient_male","meld_lab", 
                 "dummy_biweekly_dialysis",  donor_vars)

w_df$strat <- paste(sub(".*\\_", "", w_df$strata), w_df$treatment, sep="_")
fmla0 <- formula(paste("Surv(tstart, tstop, event) ~", paste(cox_wl_vars, collapse = "+"),
                      "+ strata(strat)"))
cox_fit0 <- coxph(fmla0,ties = "breslow", data = w_df, weights = weights)
bh0 <- basehaz(cox_fit0, centered = FALSE)
coeffs0 <- cox_fit0$coefficients
summ_cox0 <- summary(cox_fit0)
w1_df$strat <- w1_df$rec_country
fmla1 <- formula(paste("Surv(tstop, status) ~", paste(cox_pt_vars, collapse = "+"),
                       "+ strata(strat)"))
cox_fit1 <- coxph(fmla1, ties = "breslow", data = w1_df, weights = weights)
bh1 <- basehaz(cox_fit1, centered = FALSE)
coeffs1 <- cox_fit1$coefficients
summ_cox1 <- summary(cox_fit1)
save(bh0, bh1, coeffs0, coeffs1, file = here::here("survival","rda_output","cox_w.rda"))
save(summ_cox0, summ_cox1, file = here::here("survival","rda_output","cox_summ.rda"))

### 4) Prediction --------------------------------------------------------------
### ----------------------------------------------------------------------------
# First we create the prediction data. We predict at the same cross-sections 
# as the the development data. The recipient variables are therefore taken
# from the cross-section data. For the donor variables, we sample (with repetition) 
# from the donors in the registry. Each cross-section is matched with one donor.
pred_df <- cs_df[cs_df$tstart == 0, ] |>  
  mutate(log_twait = log(twait+1)) |>
  select(id, cs, all_of(unique(c(cox_wl_vars, "meld_int", "rec_country", "meld_lab"))))
set.seed(123) #To make the sampling reproducible
smpl <- sample(1: length(unique(cs_df$cs)), replace = TRUE)
donor_df <- longdata[longdata$txp==1,] %>%  
  select(all_of(donor_vars)) %>%
  .[smpl,] %>%
  mutate(cs = unique(cs_df$cs))
pred_df <- pred_df |> left_join(donor_df, by = "cs")
levs <- unique(bh1$strata)
pred_df$rmst0 <- NA_real_
pred_df$rmst1 <- NA_real_
for (i in 1:length(levs)) {
  lev <- levs[i] 
  ll0 <- (bh0$time < control$t_rmst*control$dt & bh0$strata==paste(lev, "0", sep = "_"))
  ll1 <- (bh1$time < control$t_rmst*control$dt & bh1$strata==lev)
  cumhaz0 <- c(0, bh0$hazard[ll0])
  dt0 <- diff(c(0, bh0$time[ll0], control$t_rmst*control$dt))
  cumhaz1 <- c(0, bh1$hazard[ll1])
  dt1 <- diff(c(0, bh1$time[ll1], control$t_rmst*control$dt))
  lp0 <- rowSums(t(t(pred_df[pred_df$rec_country == lev,cox_wl_vars]) * coeffs0))
  pred_df$rmst0[pred_df$rec_country == lev] <-  c(exp(-exp(lp0)%*%t(cumhaz0)) %*% dt0)
  lp1 <- rowSums(t(t(pred_df[pred_df$rec_country == lev,cox_pt_vars]) * coeffs1))
  pred_df$rmst1[pred_df$rec_country == lev] <- c(exp(-exp(lp1)%*%t(cumhaz1)) %*%dt1)
}
pred_df$Benefit <- pred_df$rmst1 - pred_df$rmst0

### 5) Picture -----------------------------------------------------------------
### ----------------------------------------------------------------------------
pred_df <- pred_df[!is.na(pred_df$meld_int),]
pred_df$MELD <- round(pred_df$meld_int)
pred_df$Country <- pred_df$rec_country

set.seed(123)  
p <- ggplot(pred_df, aes(x = MELD, y = Benefit)) +
  geom_jitter(width = 0.4, size = 0.85, color = "#5E81ACFF") +
  #scale_color_manual(values=c("#9EB0FFFF", "#5E81ACFF", "#8FB0BBFF"))+
  #guides(colour = guide_legend(override.aes = list(size=2)))+
  labs(
    title = "", 
    y = TeX("Benefit"),
    x = TeX("MELD"),
    caption = "") +
  ylim(c(-1,3)) +
  theme(plot.caption = element_text(hjust = 0))
theme_set(theme_bw(base_size = 14))
p


### 6) Top picks ---------------------------------------------------------------
### ----------------------------------------------------------------------------
top_ben <- ave(pred_df$Benefit, pred_df$cs, FUN = function(a) max(a)) 
topben_df <- pred_df[top_ben == pred_df$Benefit,]  |>
  group_by(cs) |>
  summarise_all(~first(.x)) |>
  ungroup()
top_meld <- ave(pred_df$MELD, pred_df$cs, FUN = function(a) max(a))
topmeld_df <- pred_df[top_meld == pred_df$MELD,]  |>
  group_by(cs) |>
  summarise_all(~first(.x)) |>
  ungroup()
5000*(round(mean(topben_df$Benefit), 2) - round(mean(topmeld_df$Benefit), 2))
