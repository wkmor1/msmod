% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msm.r
\name{msm}
\alias{msm}
\title{Fit a Multispecies Model}
\usage{
msm(y, sites, x, species, traits, data, site_re = FALSE, type = c("mstm",
  "jsdm"), method = c("glmer", "jags", "stan"), ...)
}
\arguments{
\item{y}{Character. Column name of response variable, species occurrences.}

\item{sites}{Character. Column name of site variable.}

\item{x}{Character. Column names of site level predictor variables.}

\item{species}{Character. Column name of species variable.}

\item{traits}{Character. Column names of species trait variables. Ignore}

\item{data}{A data.frame containing the variables for the model.}

\item{site_re}{Logical. Should a site level random effect be included.}

\item{type}{The type of model to fit either 'mstm', multispecies trait model,
or 'jsdm', joint species distribution model.}

\item{method}{The method used to fit the model.}

\item{...}{Further arguments to pass on to model fitting functions.}
}
\description{
Fit a Multispecies Model
}
\examples{
msm_glmer <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla', eucs)
msm_glmer_probit <- msm('present', 'plot', 'logit_rock', 'species', 'ln_sla',
  eucs, family = binomial(link = 'probit'))
msm_jags <- msm('present', 'plot', 'logit_rock', 'species', data = eucs, type = 'jsdm',
  method = 'jags')
}

