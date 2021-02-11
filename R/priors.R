# This file holds all semi hardcoded variables
#
# These can be modified with epp.set. Objects created in the package namespace
# are locked, so you can't change their values after loading. While the
# environment object is locked, its contents are still modifiable.
eppenv <- R6::R6Class("eppenv", class=F, cloneable=F, portable=F, lock_objects=F,
	public = list(
		 # r(t) prior
		 priors.rlog_pr_mean = c(log(0.35), log(0.09), log(0.2), 1980),
		 priors.rlog_pr_sd = c(0.5, 0.3, 0.5, 10),

		 # Seeding, see transf_iota 
		 priors.logiota.unif.prior = log(c(1e-13, 0.0025)),
		 priors.r0logiotaratio.unif.prior = c(-25, -5),

		 ## r-spline prior parameters
		 ## tau2.prior.rate <- 0.5  # initial sampling distribution for tau2 parameter
		 priors.tau2_init_shape = 3,
		 priors.tau2_init_rate = 4,
		 # Inverse gamma parameter for tau^2 prior for spline
		 priors.tau2_prior_shape = 0.001,   
		 priors.tau2_prior_rate = 0.001,
		 #1/duration for r steady state prior, in the code use epp::muSS not this value 
		 priors.muSS = 1/11.5,               

		 priors.rw_prior_shape = 300,
		 priors.rw_prior_rate = 1.0,
		 priors.rw_prior_sd = 0.06,

		 ## r-trend prior parameters
		 priors.t0.unif.prior = c(1970, 1990),
		 ## t1.unif.prior <- c(10, 30)
		 ## logr0.unif.prior <- c(1/11.5, 10)
		 priors.t1.pr.mean = 20.0,
		 priors.t1.pr.sd = 4.5,
		 priors.logr0.pr.mean = 0.42,
		 priors.logr0.pr.sd = 0.23,
		 ## rtrend.beta.pr.mean <- 0.0
		 ## rtrend.beta.pr.sd <- 0.2
		 priors.rtrend.beta.pr.mean = c(0.46, 0.17, -0.68, -0.038),
		 priors.rtrend.beta.pr.sd = c(0.12, 0.07, 0.24, 0.009),

		 ################################################
		 ####  Prior for incidence rate ratio model  ####
		 ################################################

		 ## RW2 model
		 NPARAM_RW2 = 13,

		 sexincrr.pr.mean = log(1.38),
		 sexincrr.pr.sd = 0.2,

		 mf_transm_rr.pr.mean = log(1.9),
		 mf_transm_rr.pr.sd = 0.3,  # change default to 0.3,

		 ## ageincrr.pr.mean <- c(-1.40707274, -0.23518703, 0.69314718, 0.78845736, -0.39975544, -0.70620810, -0.84054571, -0.02101324, -0.16382449, -0.37914407, -0.59639985, -0.82038300)
		 ## ageincrr.pr.sd <- 0.5

		 ## Informative priors based on estimates for 11 countries with 3+ surveys
		 ageincrr.pr.mean = c(-1.4, -0.28, 0.3, 0.3, -0.3, -0.6, -0.2, 0.05, -0.4, -0.45, -0.6, -0.7),
		 ageincrr.pr.sd = c(0.5, 0.4, 0.23, 0.3, 0.3, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.2),

		 NPARAM_LININCRR = 6,
		 ## incrr_trend_mean <- c(0, 0, 0, 0, 0, 0)
		 ## incrr_trend_sd <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

		 ## Informative priors based on estimates for 11 countries with 3+ surveys
		 incrr_trend_mean = c(0.0, 0.035, -0.02, -0.09, -0.016, -0.06),
		 incrr_trend_sd = c(0.07, 0.07, 0.1, 0.1, 0.08, 0.08),

		 ## Incidence rate ratios for age 50 plus, relative to 15-49
		 incrr_50plus_logdiff = cbind(male = log(4.493510) - log(c(0.358980, 0.282400, 0.259240, 0.264920, 0.254790, 0.164140, 0.000000)), 
																	female = log(2.440260) - log(c(0.336720, 0.239470, 0.167890, 0.146590, 0.171350, 0.000000, 0.000000))),

		 ## Beers coefficient matrix
		 beers_Amat = create_beers(17)[16:81, 4:17],
		 sigma_agepen = 0.1,

		 # kincrr
		 kincrr.mean = 1,
		 kincrr.sd = .5, #

		 initialize = function(fp) {
		 }
	)
)

.epp.env <- eppenv$new()

#' Set eppasm options, such as priors
#'
#' @export
epp.set <- function(name, value) {
	if (exists(".epp.env") && is.environment(.epp.env))
		assign(name, value, envir=.epp.env)
	else
		stop('likelihood calculate needs the priors set in .epp.env')
}

#' Get eppasm options, such as priors
#'
#' Trivial but for completeness sake
#' @export
epp.get <- function(name) {
	if (exists(".epp.env") && is.environment(.epp.env))
		if (exists(name, where = .epp.env))
			return(.epp.env[[name]])
	else
		stop('variable not exist or check spelling')
}

#' List eppasm options, such as priors
#'
#' Trivial but for completeness sake
#' @export
epp.ls <- function() {
	if (exists(".epp.env") && is.environment(.epp.env))
		return(sort(names(.epp.env)))
}


# Spline
.epp.env$gamS <- mgcv::gam(rr ~ -1 + s(age, bs = 'ts', pc=60),
													 data=data.frame(rr=1, age=15:80), fit=FALSE)
epp.set("spline.Q", .epp.env$gamS$smooth[[1]]$S[[1]])
epp.set("spline.X", .epp.env$gamS$X)
.epp.env$spline.X[60:80-14, ] <- 1e-14
epp.set("spline.n", nrow(.epp.env$spline.Q))
epp.set("spline.penalty", .1)
epp.set("spline.mean", 0)
epp.set("spline.stdv", 1)

.epp.env.original <- .epp.env
#' Reset eppasm options, such as priors
#'
#' @export
epp.reset <- function() {
	for (name in names(.epp.env.original))
		.epp.env[[name]] <<- .epp.env.original[[name]]
}

