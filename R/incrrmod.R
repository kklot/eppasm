## Beers coefficients to distribute infections from 5-year age groups to single-year of age
create_beers <- function(n5yr){

  ## Beer's coefficients for disaggregating 5 year age groups into
  ## single-year age groups (from John Stover)
  Afirst <- rbind(c(0.3333, -0.1636, -0.0210,  0.0796, -0.0283),
                  c(0.2595, -0.0780,  0.0130,  0.0100, -0.0045),
                  c(0.1924,  0.0064,  0.0184, -0.0256,  0.0084),
                  c(0.1329,  0.0844,  0.0054, -0.0356,  0.0129),
                  c(0.0819,  0.1508, -0.0158, -0.0284,  0.0115))
  Asecond <- rbind(c( 0.0404,  0.2000, -0.0344, -0.0128,  0.0068),
                   c( 0.0093,  0.2268, -0.0402,  0.0028,  0.0013),
                   c(-0.0108,  0.2272, -0.0248,  0.0112, -0.0028),
                   c(-0.0198,  0.1992,  0.0172,  0.0072, -0.0038),
                   c(-0.0191,  0.1468,  0.0822, -0.0084, -0.0015))
  Amid <- rbind(c(-0.0117,  0.0804,  0.1570, -0.0284,  0.0027),
                c(-0.0020,  0.0160,  0.2200, -0.0400,  0.0060),
                c( 0.0050, -0.0280,  0.2460, -0.0280,  0.0050),
                c( 0.0060, -0.0400,  0.2200,  0.0160, -0.0020),
                c( 0.0027, -0.0284,  0.1570,  0.0804, -0.0117))
  Apenult <- rbind(c(-0.0015, -0.0084,  0.0822,  0.1468, -0.0191),
                   c(-0.0038,  0.0072,  0.0172,  0.1992, -0.0198),
                   c(-0.0028,  0.0112, -0.0248,  0.2272, -0.0108),
                   c( 0.0013,  0.0028, -0.0402,  0.2268,  0.0093),
                   c( 0.0068, -0.0128, -0.0344,  0.2000,  0.0404))
  Aultim <- rbind(c( 0.0115, -0.0284, -0.0158,  0.1508,  0.0819),
                  c( 0.0129, -0.0356,  0.0054,  0.0844,  0.1329),
                  c( 0.0084, -0.0256,  0.0184,  0.0064,  0.1924),
                  c(-0.0045,  0.0100,  0.0130, -0.0780,  0.2595),
                  c(-0.0283,  0.0796, -0.0210, -0.1636,  0.3333))

  A <- do.call(rbind,
               c(list(cbind(Afirst, matrix(0, 5, n5yr-5)),
                      cbind(Asecond, matrix(0, 5, n5yr-5))),
                 lapply(0:(n5yr-6), function(i) cbind(matrix(0, 5, i), Amid, matrix(0, 5, (n5yr-5)-i))),
                 list(cbind(matrix(0, 5, n5yr-6), Apenult, matrix(0, 5, 1)),
                      cbind(matrix(0, 5, n5yr-6), Aultim, matrix(0, 5, 1)),
                      c(rep(0, n5yr-1), 1))))
  return(round(A, 4))
}

getnparam_incrr <- function(fp){
  value <- switch(as.character(fp$fitincrr),
                  "FALSE"  = 0,
                  "TRUE"   = .epp.env$NPARAM_RW2,
									sex_only = 1,
                  linincrr = .epp.env$NPARAM_RW2+.epp.env$NPARAM_LININCRR,
                  lognorm  = 7,
                  kincrr   = 1+2*length(16:60),
                  spline   = 1+.epp.env$spline.n*2,
                  relbehav = NPAR_RELBEHAV, NA)
  if(is.na(value))
    stop(paste0("fitincrr model '", fp$fitincrr, "' is not recognized"))
  value
}

transf_incrr <- function(theta_incrr, param, fp){
  incrr_nparam <- getnparam_incrr(fp)
  if(fp$incidmod == "eppspectrum"){
    param$incrr_sex <- fp$incrr_sex
    param$incrr_sex[] <- exp(theta_incrr[1])
  } else if(fp$incidmod == "transm") {
    param$mf_transm_rr <- rep(exp(theta_incrr[1]), fp$ss$PROJ_YEARS)
  }
  if(fp$fitincrr %in% c(TRUE ,"linincrr")){
    ## param$sigma_agepen <- exp(theta_incrr[incrr_nparam])
    param$sigma_agepen <- .epp.env$sigma_agepen
    param$logincrr_age <- array(0, c(14, 2))
    param$logincrr_age[c(1:2, 4:7), ] <- theta_incrr[2:13]
		param$logincrr_age[8:14, ] <- sweep(-.epp.env$incrr_50plus_logdiff, 2,
																				param$logincrr_age[7, ], "+")
    ## Smooth 5-year age group IRRs to 1-year IRRs
    incrr_age <- .epp.env$beers_Amat %*% exp(param$logincrr_age)
    incrr_age[incrr_age < 0] <- 0
    param$incrr_age <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))
    years <- with(fp$ss, proj_start+1:PROJ_YEARS-1)
    if(fp$fitincrr == "linincrr"){
      par <- theta_incrr[.epp.env$NPARAM_RW2+1:.epp.env$NPARAM_LININCRR]
      param$logincrr_trend <- par
      sexadjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[1], 0, par[2]), years, rule=2)$y
      if(fp$incidmod == "eppspectrum")
        param$incrr_sex <- param$incrr_sex * exp(sexadjust)
      else if(fp$incidmod == "transm")
        param$mf_transm_rr <- param$mf_transm_rr * exp(sexadjust)

      ## adjustment to age IRRs among 15-24
      m15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[3], 0, par[4]), years, rule=2)$y
      f15to24_adjust <- approx(c(2002, 2007, 2012), c(-5, 0, 5)*c(par[5], 0, par[6]), years, rule=2)$y
      param$incrr_age[1:10,,] <- sweep(param$incrr_age[1:10,,,drop=FALSE], 2:3, exp(rbind(m15to24_adjust, f15to24_adjust)), "*")      
    }
  } else if(fp$fitincrr=="lognorm"){
    param$logincrr_age <- cbind(calc_lognorm_logagerr(theta_incrr[2:4]),
                                calc_lognorm_logagerr(theta_incrr[5:7]))
    ## Smooth 5-year age group IRRs to 1-year IRRs
    incrr_age <- .epp.env$beers_Amat %*% exp(param$logincrr_age)
    incrr_age[incrr_age < 0] <- 0
    param$incrr_age <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))
  } else if(fp$fitincrr == "relbehav"){
    stop("relbehav is not implemented currently")
  } else if (fp$fitincrr=="kincrr") {
    param$sigma_agepen       <- .epp.env$sigma_agepen
    incrr_age                <- array(1, c(66, 2))
    incrr_age[a2i(61:80), ]  <- 0
    incrr_age[a2i(16:60), ]  <- theta_incrr[2:(1+2*length(16:60))]
    incrr_age[incrr_age < 0] <- 0
    param$incrr_age          <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))
  } else if (fp$fitincrr=="spline") {
    incrr_age         <- exp(spline_predict(theta_incrr[-1]))
    param$incrr_age   <- array(incrr_age, c(dim(incrr_age), fp$ss$PROJ_YEARS))
		param$spline.coef <- theta_incrr[-1]
	}
  return(param)
}

spline_predict <- function(theta, X = .epp.env$spline.X, npar = .epp.env$spline.n) {
		cbind(X %*% theta[1:npar], X %*% theta[(npar+1):(npar*2)])
}

lprior_incrr <- function(theta_incrr, fp){

  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  lpr <- 0
	#   Male-female ratio
  if(fp$incidmod == "eppspectrum")
    lpr <- lpr + dnorm(theta_incrr[1], 
											 .epp.env$sexincrr.pr.mean,
											 .epp.env$sexincrr.pr.sd, log=TRUE)
  else if(fp$incidmod == "transm")
    lpr <- lpr + dnorm(theta_incrr[1], 
											 .epp.env$mf_transm_rr.pr.mean,
											 .epp.env$mf_transm_rr.pr.sd, log=TRUE)
	#   Age ratio
	if(fp$fitincrr %in% c(TRUE, "linincrr")){
		lpr <- lpr + sum(dnorm(theta_incrr[2:13], .epp.env$ageincrr.pr.mean, .epp.env$ageincrr.pr.sd, log=TRUE))
		if(fp$fitincrr == "linincrr")
			lpr <- lpr + sum(dnorm(theta_incrr[.epp.env$NPARAM_RW2+1:.epp.env$NPARAM_LININCRR], 
														 .epp.env$incrr_trend_mean,
														 .epp.env$incrr_trend_sd, log=TRUE))
	} else if(fp$fitincrr=="lognorm"){
		lpr <- lpr +
			sum(dnorm(theta_incrr[c(2,5)], lognorm.a0.pr.mean, lognorm.a0.pr.sd, log=TRUE)) +
			sum(dnorm(theta_incrr[c(3,6)], lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd, log=TRUE)) +
			sum(dnorm(theta_incrr[c(4,7)], lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd, log=TRUE))
	} else if(fp$fitincrr=="relbehav") {
		lpr <- lpr + sum(dnorm(theta_incrr[2:NPAR_RELBEHAV], 0, relbehav_adjust_sd, log=TRUE));
	} else if(fp$fitincrr=="kincrr") {
		lpr <- lpr + sum(dnorm(theta_incrr[2:(1+length(16:60)*2)], .epp.env$kincrr.mean, .epp.env$kincrr.sd, log=TRUE))
	} else if(fp$fitincrr=="spline") {
		lpr <- lpr + sum(dnorm(theta_incrr[-1], .epp.env$spline.mean, .epp.env$spline.stdv, log=TRUE))
	}

  return(lpr)
}

sample_incrr <- function(n, fp){
  if(exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  incrr_nparam <- getnparam_incrr(fp)
  mat <- matrix(NA, n, incrr_nparam)
  
  if(fp$incidmod == "eppspectrum")
    mat[,1] <- rnorm(n, .epp.env$sexincrr.pr.mean, .epp.env$sexincrr.pr.sd)
  else if(fp$incidmod == "transm")
    mat[,1] <- rnorm(n, .epp.env$mf_transm_rr.pr.mean, .epp.env$mf_transm_rr.pr.sd)
  
  if(fp$fitincrr %in% c(TRUE, "linincrr")){
    mat[,2:13] <- t(matrix(rnorm(n*12, .epp.env$ageincrr.pr.mean, .epp.env$ageincrr.pr.sd), nrow=12))
    ## mat[,14] <- rnorm(n, -1, 0.7)  # log variance of ageincrr difference penalty
    if(fp$fitincrr == "linincrr")
      mat[,.epp.env$NPARAM_RW2+1:.epp.env$NPARAM_LININCRR] <- t(matrix(rnorm(n*.epp.env$NPARAM_LININCRR, .epp.env$incrr_trend_mean, .epp.env$incrr_trend_sd), nrow=.epp.env$NPARAM_LININCRR))
  } else if(fp$fitincrr=="lognorm"){
    mat[,c(2,5)] <- t(matrix(rnorm(n*2, lognorm.a0.pr.mean, lognorm.a0.pr.sd), nrow=2))
    mat[,c(3,6)] <- t(matrix(rnorm(n*2, lognorm.meanlog.pr.mean, lognorm.meanlog.pr.sd), nrow=2))
    mat[,c(4,7)] <- t(matrix(rnorm(n*2, lognorm.logsdlog.pr.mean, lognorm.logsdlog.pr.sd), nrow=2))
  } else if(fp$fitincrr=="relbehav"){
    incrr_nparam <- NPAR_RELBEHAV
    mat[,2:NPAR_RELBEHAV] <- rnorm(n*(NPAR_RELBEHAV-1), 0, relbehav_adjust_sd)
  } else if(fp$fitincrr=="kincrr") {
		mat[,2:incrr_nparam] <- rnorm(n*(incrr_nparam-1), .epp.env$kincrr.mean, .epp.env$kincrr.sd)
  } else if(fp$fitincrr=="spline") {
		mat[,2:incrr_nparam] <- rnorm(n*(incrr_nparam-1), .epp.env$spline.mean, .epp.env$spline.stdv)
  }

  return(mat)
}

ldsamp_incrr <- lprior_incrr
