# msmod

[![Build Status](https://travis-ci.org/wkmor1/msmod.svg?branch=master)](https://travis-ci.org/wkmor1/msmod)

##Fitting multispecies models

First install the software package JAGS for your platform from [this link](http://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/) 

To install run the following
```R
for (i in c('arm', 'devtools', 'dplyr', 'ggplot2', 'lme4', 'MASS', 'mclust', 'R2jags', 'rstan'))
if(!require(i, char = TRUE)) install.packages(i, char = TRUE)

devtools::install_github('wkmor1/msmod')
```
Fit a multispecies trait model
```R
library(msmod)
head(eucs)
```
```
##   plot logit_rock   ln_mrvbf ln_prec_yr ln_cv_temp ln_rad_d21 sand  loam
## 1    1  0.3846743 -0.7874579   6.665684   4.955827   8.719972 TRUE FALSE
## 2    2  1.0330150 -0.7874579   6.665684   4.955827   8.719972 TRUE FALSE
## 3    3  0.8001193 -0.7874579   6.659294   4.955827   8.087321 TRUE FALSE
## 4    4 -0.3846743 -0.7874579   6.656727   4.955827   7.400395 TRUE FALSE
## 5    5 -1.6034499 -0.7874579   6.648985   4.955827   7.209399 TRUE FALSE
## 6    6  1.2950457 -0.7874579   6.734592   4.955827   8.520882 TRUE FALSE
##   present ln_seed_wt ln_sla ln_ht species
## 1   FALSE       0.89  35.57    35     ALA
## 2   FALSE       0.89  35.57    35     ALA
## 3   FALSE       0.89  35.57    35     ALA
## 4   FALSE       0.89  35.57    35     ALA
## 5   FALSE       0.89  35.57    35     ALA
## 6   FALSE       0.89  35.57    35     ALA
```
```R
msm_glmer <- msm(y = "present", sites = "plot", x = "logit_rock",
  species = "species", traits = "ln_sla", data = eucs, type = "mstm",
  method = "glmer")

summary(msm_glmer)
```
```
## Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
##  Family: binomial  ( logit )
## Formula: present ~ logit_rock + ln_sla + logit_rock:ln_sla + ln_sla:ln_sla +  
##     (1 + logit_rock + ln_sla | species)
##    Data: data
## Control: control
## 
##      AIC      BIC   logLik deviance df.resid 
##   4660.7   4731.9  -2320.3   4640.7     9150 
## 
## Scaled residuals: 
##     Min      1Q  Median      3Q     Max 
## -0.8191 -0.3628 -0.2014 -0.0725 16.9390 
## 
## Random effects:
##  Groups  Name        Variance Std.Dev. Corr       
##  species (Intercept)  0.7140  0.8450              
##          logit_rock   0.6507  0.8066    0.15      
##          ln_sla      11.8727  3.4457   -0.98  0.07
## Number of obs: 9160, groups:  species, 20
## 
## Fixed effects:
##                   Estimate Std. Error z value Pr(>|z|)    
## (Intercept)        -2.5200     0.2485 -10.142  < 2e-16 ***
## logit_rock         -0.8069     0.2563  -3.148  0.00165 ** 
## ln_sla             -1.4362     0.8390  -1.712  0.08694 .  
## logit_rock:ln_sla  -4.7811     0.6000  -7.969  1.6e-15 ***
## ---
## Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
## 
## Correlation of Fixed Effects:
##             (Intr) lgt_rc ln_sla
## logit_rock   0.146              
## ln_sla      -0.769  0.147       
## lgt_rck:ln_  0.122  0.292  0.123
```
```R
coef_plot(msm_glmer)
```
![Coefficient plot](inst/coef_plot.png?raw=true)
```R
te_plot(msm_glmer)
```
![Trait-environment plot](inst/te_plot.png?raw=true)
