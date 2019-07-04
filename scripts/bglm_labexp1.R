library(tidyverse)
library(rstanarm)
library(loo)
source("functions/functions.R")

# Load DF
load(file="data/VoiceLineUp.Rda")

# Assign contrasts
d$TP <- factor(d$TP)
contrasts(d$TP) <- contr.sum(2)

d$exposure <- factor(d$exposure)
contrasts(d$exposure) <- contr.sum(2)

d$order <- ifelse(d$order == 3, -1, 
                  ifelse(d$order == 7, 1, 0))

d %>% mutate(acc = acc*100) %>% group_by(TP, exposure, order) %>%
  summarise(M = mean(acc),
            SD = sd(acc))

# Fit model
nchains = ncores = 3
iter = 2000

# Null model
# random intercept and slodes for items (called lineup)
m <- stan_glmer(formula = acc ~ 1
                   + (TP*order*exposure|lineup) 
                   , prior_intercept = student_t(df = 4, location = 0), # Priors
                   prior = student_t(df = 4, location = 0),
                   prior_covariance = decov(regularization = 2),
                   data = d,
                   chains = nchains,
                   iter = iter,
                   cores = ncores,
                   seed = 22,
                   adapt_delta = .99,
                   family = binomial("probit"));m
# Save model
save(m,
     file="stanout/BGLMNULL2000iter3chains.rda",
     compress="xz")


# Order model
m.order <- update(m, formula. = . ~ . + order)
# Save model
save(m.order,
     file="stanout/BGLMorder2000iter3chains.rda",
     compress="xz")


# TP model
m.tp <- update(m, formula. = . ~ . + TP)
# Save model
save(m.tp,
     file="stanout/BGLMTP2000iter3chains.rda",
     compress="xz")


# Exposure model
m.exposure <- update(m, formula. = . ~ . + exposure)
# Save model
save(m.exposure,
     file="stanout/BGLMexposure2000iter3chains.rda",
     compress="xz")


# TP and exposure model (main effects model)
m.tp.exposure <- update(m.tp, formula. = . ~ . + exposure)
# Save model
save(m.tp.exposure,
     file="stanout/BGLMTPexposure2000iter3chains.rda",
     compress="xz")


# TP and exposure model (interaction model)
m.tp.exposure.inter <- update(m.tp.exposure, formula. = . ~ . + TP : exposure)
# Save model
save(m.tp.exposure.inter,
     file="stanout/BGLMTPexposureinter2000iter3chains.rda",
     compress="xz")

#m.tp.exposure.order <- update(m.tp.exposure, formula. = . ~ . + order)
#m.tp.exposure.order.inter1 <- update(m.tp.exposure, formula. = . ~ . + TP : exposure)
#m.tp.exposure.order.inter2 <- update(m.tp.exposure, formula. = . ~ . + TP : order)
#m.tp.exposure.order.inter3 <- update(m.tp.exposure, formula. = . ~ . + exposure : order)
#m.tp.exposure.order.inter4 <- update(m.tp.exposure, formula. = . ~ . + exposure : order : TP)


# Model comparisons
# Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) 
#load(file="stanout/BGLMNULL2000iter3chains.rda")
#load(file="stanout/BGLMTP2000iter3chains.rda")
#load(file="stanout/BGLMexposure2000iter3chains.rda")
#load(file="stanout/BGLMTPexposure2000iter3chains.rda")
#load(file="stanout/BGLMTPexposureinter2000iter3chains.rda")
loo.m <- loo(m)
loo.m.order <- loo(m.order)
loo.m.tp <- loo(m.tp)
loo.m.exposure <- loo(m.exposure)
loo.m.tp.exposure <- loo(m.tp.exposure)
loo.m.tp.exposure.inter <- loo(m.tp.exposure.inter, k_threshold = 0.7)

#loo.m.tp.exposure.order.inter1 <- loo(m.tp.exposure.order.inter1)
#loo.m.tp.exposure.order.inter2 <- loo(m.tp.exposure.order.inter2)
#loo.m.tp.exposure.order.inter3 <- loo(m.tp.exposure.order.inter3)
#loo.m.tp.exposure.order.inter4 <- loo(m.tp.exposure.order.inter4)

# verify that the posterior is not too sensitive to any particular observation in the dataset.
par(mfrow = c(3,2))#mar = c(5,3.8,1,0) + 0.1, las = 3
plot(loo.m, label_points = TRUE)
plot(loo.m.order, label_points = TRUE)
plot(loo.m.tp, label_points = TRUE)
plot(loo.m.exposure, label_points = TRUE)
plot(loo.m.tp.exposure, label_points = TRUE)
plot(loo.m.tp.exposure.inter, label_points = TRUE)
#plot(loo.m.tp.exposure.order.inter1, label_points = TRUE)
#plot(loo.m.tp.exposure.order.inter2, label_points = TRUE)
#plot(loo.m.tp.exposure.order.inter3, label_points = TRUE)
#plot(loo.m.tp.exposure.order.inter4, label_points = TRUE)


# There are a couple of moderate outliers (whose statistics are greater than 0.5),
# should not have too much of an effect on the resulting model comparison (at least not model 2 and 3)
#compare_models(loo.m,loo.m.order)
compare_models(loo.m,loo.m.tp)
#compare_models(loo.m.order,loo.m.tp)
compare_models(loo.m,loo.m.exposure)
compare_models(loo.m.tp,loo.m.exposure)
compare_models(loo.m.tp,loo.m.tp.exposure)
compare_models(loo.m.tp.exposure,loo.m.tp.exposure.inter)

#compare_models(loo.m,loo.m.tp.exposure.order.inter1)
#compare_models(loo.m,loo.m.tp.exposure.order.inter2)
#compare_models(loo.m,loo.m.tp.exposure.order.inter3)
#compare_models(loo.m,loo.m.tp.exposure.order.inter4)


# Difference between expected log pointwise deviance between the two models: difference after taking into account that the second model estimates an additional parameter. The “LOO Information Criterion (LOOIC)”
compare_models(loo.m,loo.m.tp,loo.m.exposure,loo.m.tp.exposure,loo.m.tp.exposure.inter)


# WAIC
#(waic1 <- waic(m.tp.exposure))
#(waic2 <- waic(log_lik_2))
#print(compare(waic1, waic2), digits = 2)

#yrep <- posterior_predict(m)
#ppc_loo_pit_overlay(
#  y = d$acc,
#  yrep = yrep,
#  lw = weights(loo.m$psis_object)
#)



