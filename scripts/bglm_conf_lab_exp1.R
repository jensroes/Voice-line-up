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

d$acc <- factor(d$acc)
contrasts(d$acc) <- contr.sum(2)

d$exposure <- factor(d$exposure)
contrasts(d$exposure) <- contr.sum(2)

d %>% group_by(TP, exposure,  acc) %>%
  summarise(M = mean(conf, na.rm = T))

# Fit model
nchains = ncores = 3
iter = 2000

# Null model
# random intercept and slodes for items (called lineup)
m <- stan_glmer(formula = conf ~ 1
                   + (TP*order*exposure*acc|lineup) 
                   , prior_intercept = student_t(df = 4, location = 1), # Priors
                   prior = student_t(df = 4, location = 0),
                   prior_covariance = decov(regularization = 2),
                   data = d,
                   chains = nchains,
                   iter = iter,
                   cores = ncores,
                   seed = 22,
                   adapt_delta = .99,
                   family = neg_binomial_2)

# Save model
save(m,
     file="stanout/CONFBGLMNULL2000iter3chains.rda",
     compress="xz")

# acc model
m.acc <- update(m, formula. = . ~ . + acc)
# Save model
save(m.acc,
     file="stanout/CONFBGLMacc2000iter3chains.rda",
     compress="xz")


# TP model
m.tp <- update(m, formula. = . ~ . + TP)
# Save model
save(m.tp,
     file="stanout/CONFBGLMTP2000iter3chains.rda",
     compress="xz")


# Exposure model
m.exposure <- update(m, formula. = . ~ . + exposure)
# Save model
save(m.exposure,
     file="stanout/CONFBGLMexposure2000iter3chains.rda",
     compress="xz")



# All main effects
m.tp.acc.exposure <- update(m.tp, formula. = . ~ . + acc + exposure)
# Save model
save(m.tp.acc.exposure,
     file="stanout/CONFBGLMTPaccexposure2000iter3chains.rda",
     compress="xz")


# All 2-way interactions
m.2.way.inter <- update(m.tp.acc.exposure, formula. = . ~ . + acc : TP + acc : exposure + exposure : TP)
# Save model
save(m.2.way.inter,
     file="stanout/CONFBGLM2wayinters2000iter3chains.rda",
     compress="xz")

# Three way interaction
m.3.way.inter <- update(m.2.way.inter, formula. = . ~ . + acc : TP : exposure )
# Save model
save(m.3.way.inter,
     file="stanout/CONFBGLM3wayinters2000iter3chains.rda",
     compress="xz")


# Compare main effects models
loo.m <- loo(m)
loo.m.acc <- loo(m.acc)
loo.m.tp <- loo(m.tp)
loo.m.exposure <- loo(m.exposure)
loo.m.tp.acc.exposure <- loo(m.tp.acc.exposure)
loo.m.2.way.inter <- loo(m.2.way.inter)
loo.m.3.way.inter <- loo(m.3.way.inter)

compare(loo.m,loo.m.acc,loo.m.tp,loo.m.exposure, loo.m.tp.acc.exposure, loo.m.2.way.inter, loo.m.3.way.inter)


# verify that the posterior is not too sensitive to any particular observation in the dataset.
par(mfrow = c(3,3))#mar = c(5,3.8,1,0) + 0.1, las = 3
plot(loo.m, label_points = TRUE)
plot(loo.m.acc, label_points = TRUE)
plot(loo.m.tp, label_points = TRUE)
plot(loo.m.exposure, label_points = TRUE)
plot(loo.m.tp.acc.exposure, label_points = TRUE)
plot(loo.m.2.way.inter, label_points = TRUE)
plot(loo.m.3.way.inter, label_points = TRUE)



# Exposure model
#m.exposure <- update(m, formula. = . ~ . + exposure)
# Save model
#save(m.exposure,
#     file="stanout/CONFBGLMexposure2000iter3chains.rda",
#     compress="xz")




# TP and exposure model (main effects model)
#m.tp.exposure <- update(m.tp, formula. = . ~ . + exposure)
# Save model
#save(m.tp.exposure,
#     file="stanout/CONFBGLMTPexposure2000iter3chains.rda",
#     compress="xz")

# TP and exposure model (interaction model)
#m.tp.exposure.acc.inter1 <- update(m.tp.exposure.acc, formula. = . ~ . + TP : exposure)
#m.tp.exposure.acc.inter2 <- update(m.tp.exposure.acc.inter1, formula. = . ~ . + TP : acc)
#m.tp.exposure.acc.inter3 <- update(m.tp.exposure.acc.inter2, formula. = . ~ . + exposure : acc)
#m.tp.exposure.acc.inter4 <- update(m.tp.exposure.acc.inter3, formula. = . ~ . + exposure : acc : TP)



# Model comparisons
# Pareto smoothed importance-sampling leave-one-out cross-validation (PSIS-LOO) 
loo.m <- loo(m)
loo.m.order <- loo(m.order)
loo.m.tp <- loo(m.tp)
loo.m.acc <- loo(m.acc)
loo.m.tp.exposure <- loo(m.tp.exposure)
loo.m.tp.exposure.acc <- loo(m.tp.exposure.acc)
loo.m.tp.exposure.acc.inter1 <- loo(m.tp.exposure.acc.inter1)
loo.m.tp.exposure.acc.inter2 <- loo(m.tp.exposure.acc.inter2)
loo.m.tp.exposure.acc.inter3 <- loo(m.tp.exposure.acc.inter3)
loo.m.tp.exposure.acc.inter4 <- loo(m.tp.exposure.acc.inter4)

# verify that the posterior is not too sensitive to any particular observation in the dataset.
par(mfrow = c(4,4))#mar = c(5,3.8,1,0) + 0.1, las = 3
plot(loo.m, label_points = TRUE)
plot(loo.m.acc, label_points = TRUE)
plot(loo.m.tp, label_points = TRUE)
plot(loo.m.exposure, label_points = TRUE)
plot(loo.m.tp.exposure, label_points = TRUE)
plot(loo.m.tp.exposure.acc, label_points = TRUE)
plot(loo.m.tp.exposure.acc.inter1, label_points = TRUE)
plot(loo.m.tp.exposure.acc.inter2, label_points = TRUE)
plot(loo.m.tp.exposure.acc.inter3, label_points = TRUE)
plot(loo.m.tp.exposure.acc.inter4, label_points = TRUE)


# There are a couple of moderate outliers (whose statistics are greater than 0.5),
# should not have too much of an effect on the resulting model comparison (at least not model 2 and 3)
#compare_models(loo.m,loo.m.order)
compare_models(loo.m,loo.m.tp)
compare_models(loo.m,loo.m.acc)
compare_models(loo.m,loo.m.exposure)
compare_models(loo.m,loo.m.tp.exposure)
compare_models(loo.m,loo.m.tp.exposure.acc)
compare_models(loo.m,loo.m.tp.exposure.acc.inter1)
compare_models(loo.m,loo.m.tp.exposure.acc.inter2)
compare_models(loo.m,loo.m.tp.exposure.acc.inter3)
compare_models(loo.m,loo.m.tp.exposure.acc.inter4)

# Difference between expected log pointwise deviance between the two models: difference after taking into account that the second model estimates an additional parameter. The “LOO Information Criterion (LOOIC)”
compare_models(loo.m,loo.m.tp,
               loo.m.exposure,
               loo.m.acc, 
               loo.m.tp.exposure,
               loo.m.tp.exposure.acc,
               loo.m.tp.exposure.acc.inter1,
               loo.m.tp.exposure.acc.inter2,
               loo.m.tp.exposure.acc.inter3,
               loo.m.tp.exposure.acc.inter4)


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



