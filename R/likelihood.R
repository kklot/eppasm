###########################
####  EPP model prior  ####
###########################

ldinvgamma <- function(x, alpha, beta){
  log.density <- alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - (beta/x)
  return(log.density)
}

bayes_lmvt <- function(x, shape, rate){
  mvtnorm::dmvt(x, sigma=diag(length(x)) / (shape / rate), df=2*shape, log=TRUE)
}

bayes_rmvt <- function(n, d, shape, rate){
  mvtnorm::rmvt(n, sigma=diag(d) / (shape / rate), df=2*shape)
}

## Binomial distribution log-density permitting non-integer counts
ldbinom <- function(x, size, prob){
  lgamma(size+1) - lgamma(x+1) - lgamma(size-x+1) + x*log(prob) + (size-x)*log(1-prob)
}

###################################
####  Age/sex incidence model  ####
###################################

fnCreateParam <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }
  
  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1
    if (fp$eppmod == "rspline"){
      u <- theta[1:fp$numKnots]
      if (fp$rtpenord == 2){
        beta <- numeric(fp$numKnots)
        beta[1] <- u[1]
        beta[2] <- u[1]+u[2]
        for(i in 3:fp$numKnots)
          beta[i] <- -beta[i-2] + 2*beta[i-1] + u[i]
      } else # first order penalty
        beta <- cumsum(u)
    } else if (fp$eppmod %in% "logrw")
      beta <- theta[1:fp$numKnots]

    param <- list(beta = beta,
                  rvec = as.vector(fp$rvec.spldes %*% beta))
    
    if (fp$eppmod %in% "logrw")
      param$rvec <- exp(param$rvec)

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      param$iota <- exp(param$rvec[fp$proj.steps == fp$tsEpidemicStart] * theta[fp$numKnots+1])
    else
      param$iota <- transf_iota(theta[fp$numKnots+1], fp)
  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    par <- theta[1:4]
    par[3] <- exp(theta[3])
    param <- list()
    param$rvec <- exp(rlogistic(fp$proj.steps, par))
    param$iota <- transf_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend
    epp_nparam <- 7
    param <- list(tsEpidemicStart = fp$proj.steps[which.min(abs(fp$proj.steps - (round(theta[1]-0.5)+0.5)))], # t0
                  rtrend = list(tStabilize = round(theta[1]-0.5)+0.5+round(theta[2]),  # t0 + t1
                                r0 = exp(theta[3]),              # r0
                                beta = theta[4:7]))
  } else {
    epp_nparam <- fp$rt$n_param+1
    param <- list()
    param$rvec <- create_rvec(theta[1:fp$rt$n_param], fp$rt)
    param$iota <- transf_iota(theta[fp$rt$n_param+1], fp)
  }

  nparam <- epp_nparam
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    ancmod.id <- nparam + 1:fp$ancmod$nparam
    nparam <- nparam + fp$ancmod$nparam
    param <- update_par(param, list = create_ancmod_param(theta[ancmod.id], fp$ancmod))
    param$frr_cd4 <- fp$frr_cd4 * exp(param$log_frr_adjust)
    param$frr_art <- fp$frr_art * exp(param$log_frr_adjust)
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE){
    fitincrr.id <- nparam + 1:getnparam_incrr(fp)
    nparam <- nparam + getnparam_incrr(fp)
    param  <- transf_incrr(theta[fitincrr.id], param, fp)
  }

  if (fp$ss$MODEL == 2) {
    if (exists("fitpcr", fp) && fp$fitpcr==TRUE) {
      fitpcr.id <- nparam + 1:8
      param$est_pcr <- sexual_rate(theta[fitpcr.id])
    } else
      param$est_pcr <- fp$est_pcr
    
    param$leading_ev <- NGM(fp$basepop, fp$db_rate[,,1], fp$mixmat,
                            param$est_pcr, 1-fp$Sx[,,1],
                            param$mf_transm_rr[[1]], param$rvec[1])

  }

  return(param)
}

#' Next generation matrix for the sexual mixing model
#'
#' this is simplified as the absolute size is not important to the relative size
#'
#' @param base_pop initial population
#' @param debut_rate sexual debut rate, need to adjust to included in F?
#' @param mixing_matrice 3-dimensions array size 66x66x2
#' @param partner_rate partner acquisition rate
#' @param death_rate mortality rate, currently using 1-survival
#' @param male_female_rate relative male to female transmission compared to other way around (not important)
#' @param rt_0 initial growth rate of the r_t model (not important)
#' @export 
NGM <- function(base_pop, debut_rate, mixing_matrice, partner_rate, death_rate, male_female_rate, rt_0) {
	base_pop[1:16, ] = base_pop[1:16, ] * debut_rate
  nc_m <- sweepx(mixing_matrice[,,1], 1, partner_rate[, 1])
  nc_f <- sweepx(mixing_matrice[,,2], 1, partner_rate[, 2])
  nc_m_total <- sweepx(nc_m, 1, base_pop[,1])
  nc_f_total <- sweepx(nc_f, 1, base_pop[,2])
  ratio_mf <- nc_m_total / t(nc_f_total)
  nc_m_adj <- nc_m * ratio_mf^(-0.5)
  nc_f_adj <- nc_f * t(ratio_mf)^0.5
  M1 = Matrix(1, 1, 66)
	Ffm = nc_m_adj * base_pop[,1] %*% M1 * (Matrix::t(M1) %*% (1/base_pop[,2]))
	Fmf = nc_f_adj * base_pop[,2] %*% M1 * (Matrix::t(M1) %*% (1/base_pop[,1]))
  Ffm[is.na(Ffm) | !is.finite(Ffm)] = min(Ffm, na.rm=TRUE)
  Fmf[is.na(Fmf) | !is.finite(Fmf)] = min(Fmf, na.rm=TRUE)
	FF = (kronecker(Matrix(c(0,1,rep(0, 2)), nrow=2), Fmf) +
				kronecker(Matrix(c(rep(0, 2),1,rep(0, 1)), nrow=2), Ffm) )
	fci = function(mui, q=0, adiff=1) (mui + q) / (exp((mui+q)*adiff) - 1)
	gm = fci(death_rate[,1])
	gf = fci(death_rate[,2])
	Vm = Diagonal(66, gm + death_rate[,1])
	Vf = Diagonal(66, gf + death_rate[,2])
	Vm[cbind(2:66, 1:65)] = - gm[-1]
	Vf[cbind(2:66, 1:65)] = - gf[-1]
	V = (kronecker(Matrix(c(1, rep(0, 3)), nrow=2), Vm) +
			 kronecker(Matrix(c(rep(0, 3), 1), nrow=2), Vf))
	NGM = FF %*% Matrix::solve(V)
	eig = eigen(NGM)
	eid = which.max(Re(eig$values))
  #   eig$vectors[,eid] %>% Re %>% { ./sum(.) } %>% matrix(66,2) %>% matplot
	eig$vectors[,eid] %>% Re %>% { ./sum(.) } %>% matrix(66,2)
}




########################################################
####  Age specific prevalence likelihood functions  ####
########################################################


#' Prepare age-specific HH survey prevalence likelihood data
#' 
#' @param hhsage age-specific HH survey prevalence likelihood data
#' @param fp fix parameters
prepare_hhsageprev_likdat <- function(hhsage, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsage$W.hhs <- qnorm(hhsage$prev)
  hhsage$v.hhs <- 2*pi*exp(hhsage$W.hhs^2)*hhsage$se^2
  hhsage$sd.W.hhs <- sqrt(hhsage$v.hhs)

  if(exists("deff_approx", hhsage))
    hhsage$n_eff <- hhsage$n/hhsage$deff_approx
  else if(exists("deff_approx", hhsage))
    hhsage$n_eff <- hhsage$n/hhsage$deff
  else
    hhsage$n_eff <- hhsage$prev * (1 - hhsage$prev) / hhsage$se ^ 2
  hhsage$x_eff <- hhsage$n_eff * hhsage$prev

  if(is.null(hhsage$sex))
    hhsage$sex <- rep("both", nrow(hhsage))

  if(is.null(hhsage$agegr))
    hhsage$agegr <- "15-49"

  startage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\1", hhsage$agegr))
  endage <- as.integer(sub("([0-9]*)-([0-9]*)", "\\2", hhsage$agegr))
  
  hhsage$sidx <- match(hhsage$sex, c("both", "male", "female")) - 1L
  hhsage$aidx <- startage - fp$ss$AGE_START+1L
  hhsage$yidx <- as.integer(hhsage$year - (anchor.year - 1))
  hhsage$agspan <- endage - startage + 1L

  return(subset(hhsage, aidx > 0))
}

#' Log likelihood for age-specific household survey prevalence
#' 
#' @param mod model simulation output
#' @param dat Output data from prepare_likdat
#' @param pointwise Point-wise likelihood
ll_hhsage <- function(mod, dat, pointwise = FALSE){
  qM.age <- suppressWarnings(qnorm(ageprev(mod,
                                           aidx = dat$aidx,
                                           sidx = dat$sidx, 
                                           yidx = dat$yidx, 
                                           agspan = dat$agspan)))
  if (any(is.na(qM.age)))
    val <- rep(-Inf, nrow(dat))
  else
    val <- dnorm(dat$W.hhs, qM.age, dat$sd.W.hhs, log=TRUE)

  if (pointwise)
    return(val)
  sum(val)
}

#' Log likelihood for age-specific household survey prevalence using binomial approximation
#' 
#' @param mod model simulation output
#' @param dat Output data from prepare_likdat
#' @param pointwise Point-wise likelihood
ll_hhsage_binom <- function(mod, dat, pointwise = FALSE){

  prevM.age <- suppressWarnings(ageprev(mod, dat$aidx, dat$sidx, dat$yidx, dat$agspan))

  if (any(is.na(prevM.age)) || any(prevM.age >= 1))
    val <- rep(-Inf, nrow(dat))
  else
    val <- ldbinom(dat$x_eff, dat$n_eff, prevM.age)
  val[is.na(val)] <- -Inf

  if (pointwise)
    return(val)

  sum(val)
}



#########################################
####  Incidence likelihood function  ####
#########################################

#' Prepare household survey incidence likelihood data
#' 
#' @param hhsincid household survey incidence likelihood data
#' @param fp fix parameters
prepare_hhsincid_likdat <- function(hhsincid, fp){
  anchor.year <- floor(min(fp$proj.steps))

  hhsincid$idx <- hhsincid$year - (anchor.year - 1)
  hhsincid$log_incid <- log(hhsincid$incid)
  hhsincid$log_incid.se <- hhsincid$se/hhsincid$incid

  return(hhsincid)
}

#' Log-likelhood for direct incidence estimate from household survey
#'
#' Calculate log-likelihood for nationally representative incidence
#' estimates from a household survey. Currently implements likelihood
#' for a log-transformed direct incidence estimate and standard error.
#' Needs to be updated to handle incidence assay outputs.
#'
#' @param mod model output, object of class `spec`.
#' @param hhsincid.dat prepared houshold survey incidence estimates (see perp
ll_hhsincid <- function(mod, hhsincid.dat){
  logincid <- log(incid(mod$data, fp))
  ll.incid <- sum(dnorm(hhsincid.dat$log_incid, logincid[hhsincid.dat$idx], hhsincid.dat$log_incid.se, TRUE))
  return(ll.incid)
}


###############################
####  Likelihood function  ####
###############################

prepare_likdat <- function(eppd, fp){

  likdat <- list()
  
  likdat$hhs.dat <- prepare_hhsageprev_likdat(eppd$hhs, fp)

  if (exists("ancsitedat", where=eppd) && nrow(eppd$ancsitedat)) {

    ancsitedat <- eppd$ancsitedat
    
    if (exists("ancmod", fp) && !fp$ancmod$has_ancrtsite)
      ancsitedat <- subset(ancsitedat, type == "ancss")
    
    likdat$ancsite.dat <- prepare_ancsite_likdat(ancsitedat, fp)
  }
 
  if (exists("ancrtcens", where=eppd)){
    if (exists("ancmod", fp) && !fp$ancmod$has_ancrtcens)
      eppd$ancrtcens <- NULL
    else
      likdat$ancrtcens.dat <- prepare_ancrtcens_likdat(eppd$ancrtcens, fp)
  }

  if (exists("hhsincid", where=eppd))
    likdat$hhsincid.dat <- prepare_hhsincid_likdat(eppd$hhsincid, fp)

  return(likdat)
}

lprior <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1

    nk <- fp$numKnots

    if (fp$eppmod == "logrw")
      lpr <- bayes_lmvt(theta[2:fp$numKnots], .epp.env$priors.rw_prior_shape, .epp.env$priors.rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[(1+fp$rtpenord):nk], .epp.env$priors.tau2_prior_shape, .epp.env$priors.tau2_prior_rate)

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif (theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + lprior_iota(theta[nk+1], fp)
  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend    
    epp_nparam <- 7
    
    lpr <- dunif (round(theta[1]), .epp.env$priors.t0.unif.prior[1], .epp.env$priors.t0.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), .epp.env$priors.t1.pr.mean, .epp.env$priors.t1.pr.sd, log=TRUE) +
      dnorm(theta[3], .epp.env$priors.logr0.pr.mean, .epp.env$priors.logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], .epp.env$priors.rtrend.beta.pr.mean, .epp.env$priors.rtrend.beta.pr.sd, log=TRUE))
  } else if (fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd, log=TRUE)) +
      sum(dnorm(theta[4+1:fp$rt$n_rw], 0, .epp.env$priors.rw_prior_sd, log=TRUE))
    lpr <- lpr + lprior_iota(theta[fp$rt$n_param+1], fp)
  }

  nparam <- epp_nparam
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    ancmod.id <- nparam + 1:fp$ancmod$nparam
    nparam <- nparam + fp$ancmod$nparam
    lpr <- lpr + lprior_ancmod(theta[ancmod.id], fp$ancmod, fp$prior_args)
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE){
    fitincrr.id  <- nparam + 1:getnparam_incrr(fp)
    nparam <- nparam + getnparam_incrr(fp)
    lpr <- lpr + lprior_incrr(theta[fitincrr.id], fp)
  }

  if(exists("fitpcr", where=fp)){
    fitpcr.id <- nparam + 1:8
    nparam <- nparam + 8
    lpr <- lpr + sexual_rate_prior(theta[fitpcr.id], fp)
  }

  return(lpr)
}

lgt <- function(x, initx, maxx, midx, nx) initx * (1 - maxx) + initx * maxx / (1 + (x/midx)^nx)
lgt_ <- function(x, p) p[1] * (1 - p[2]) + p[1] * p[2] / (1 + (x/p[3])^p[4])

# lgt2p <- function(x, K, M) 1 + (K-1)/(1+exp(x - M))
lgt2p <- function(x, p) 1 + (p[1]-1)/(1+exp(x - p[2]))

lgt_prior <- function(x) {
    dgamma(x[1], 1, scale=2, log=TRUE) + # increase hist(rgamma(3000, 1, scale=2))
    dgamma(x[2], 40, scale=0.6, log=TRUE) # # hist(rgamma(3000, 40, scale=.5))
}

lgt_sample <- function() {
    c(
      rgamma(1, 1, scale=2) , # increase hist(rgamma(3000, 1, scale=2))
      rgamma(1, 40, scale=0.6) # # hist(rgamma(3000, 40, scale=.5))
    )
}

ll_spline <- function(theta, 
											Q       = .epp.env$spline.Q,
											penalty = .epp.env$spline.penalty,
											N       = .epp.env$spline.n) {
	2 * 0.5 * N  * log(penalty) - 
		0.5 * penalty * theta[1:N] %*% (Q %*% theta[1:N]) - 
		0.5 * penalty * theta[1:N+N] %*% (Q %*% theta[1:N+N])
}

#' Full log-likelhood
#' 
#' @importFrom stats aggregate approx cov cov.wt density dexp dlnorm dnorm dunif ecdf mahalanobis median model.matrix na.omit optim optimHess pnorm qnorm quantile relevel rexp rgamma rnorm runif sd setNames update var
ll_all = function(theta, fp, likdat) {

  .epp.env$theta.last <<- theta

  nparam <- length(theta)
  
  fp <- update(fp, keep.attr=FALSE, list=fnCreateParam(theta, fp))

  if (fp$eppmod == "rspline")
    if (any(is.na(fp$rvec)) || min(fp$rvec) < 0 || max(fp$rvec) > 20) 
      return(-Inf)

  mod <- simmod(fp)

  .anc = .ancrt = .hhs = .incid = .rprior = .incpen = 0

  if (exists("fitincrr", where=fp)) {
    if (fp$fitincrr==TRUE)
      .incpen <- sum(dnorm(diff(fp$logincrr_age[1:7, ], differences=2), sd=fp$sigma_agepen, log=TRUE))
    if (fp$fitincrr=="kincrr")
			.incpen <- sum(dnorm(diff(fp$incrr_age[,,1], differences=2), sd=fp$sigma_agepen, log=TRUE))
		if (fp$fitincrr=="spline")
			.incpen <- ll_spline(fp$spline.coef)
  }

  ## ANC likelihood
  if (exists("ancsite.dat", likdat))
    .anc <- ll_ancsite(mod, fp, coef=c(fp$ancbias, fp$ancrtsite.beta),
                       vinfl=fp$v.infl, likdat$ancsite.dat)

  if (exists("ancrtcens.dat", likdat))
    .ancrt <- ll_ancrtcens(mod, likdat$ancrtcens.dat, fp)

  ## Household survey likelihood
  if (exists("hhs.dat", where=likdat)) {
    if (exists("ageprev", fp) && fp$ageprev=="binom")
      .hhs <- ll_hhsage_binom(mod$data, likdat$hhs.dat)
    else ## use probit likelihood
      .hhs <- ll_hhsage(mod$data, likdat$hhs.dat) # probit-transformed model
  }

  if (!is.null(likdat$hhsincid.dat))
    .incid <- ll_hhsincid(mod, likdat$hhsincid.dat)

  if (exists("equil.rprior", where=fp) && fp$equil.rprior) {
    if (fp$eppmod != "rspline")
      stop("error in ll(): equil.rprior is only for use with r-spline model")

    lastdata.idx <- max(likdat$ancsite.dat$df$yidx,
                        likdat$hhs.dat$yidx,
                        likdat$ancrtcens.dat$yidx,
                        likdat$hhsincid.dat$idx)
    
    qM.all <- suppressWarnings(qnorm(mod$prev15to49))

    if (any(is.na(qM.all[lastdata.idx - 9:0]))) {
      .rprior <- -Inf
    } 
    else {
      rvec.ann <- fp$rvec[fp$proj.steps %% 1 == 0.5]
      equil.rprior.mean <- epp:::muSS/(1-pnorm(qM.all[lastdata.idx]))
      if (!is.null(fp$prior_args$equil.rprior.sd))
        equil.rprior.sd <- fp$prior_args$equil.rprior.sd
      else
        equil.rprior.sd <- sqrt(mean((epp:::muSS / (1-pnorm(qM.all[lastdata.idx - 9:0])) - rvec.ann[lastdata.idx - 9:0])^2))  # empirical sd based on 10 previous years
      
      .rprior <- sum(dnorm(rvec.ann[(lastdata.idx+1L):length(qM.all)], equil.rprior.mean, equil.rprior.sd, log=TRUE))  # prior starts year after last data
    }
  }

  c(anc = .anc, ancrt = .ancrt, hhs = .hhs, incid = .incid, rprior = .rprior, incpen = .incpen) 
}


##########################
####  IMIS functions  ####
##########################

sample.prior <- function(n, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  ## Calculate number of parameters
  if (fp$eppmod %in% c("rspline", "logrw"))
    epp_nparam <- fp$numKnots+1L
  else if (fp$eppmod == "rlogistic")
    epp_nparam <- 5
  else if (fp$eppmod == "rtrend")
    epp_nparam <- 7
  else if (fp$eppmod == "rhybrid")
    epp_nparam <- fp$rt$n_param+1

  nparam <- epp_nparam
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    ancmod.id <- nparam + 1:fp$ancmod$nparam
    nparam <- nparam + fp$ancmod$nparam
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE) {
    fitincrr.id <- nparam + 1:getnparam_incrr(fp)
    nparam <- nparam + getnparam_incrr(fp)
  }
  
  if(exists("fitpcr", where=fp)) {
    fitpcr.id <- nparam + 1:8
    nparam <- nparam + 8
  }

  ## Create matrix for storing samples
  mat <- matrix(NA, n, nparam)

  if (fp$eppmod %in% c("rspline", "logrw")){
    if (fp$eppmod == "rspline")
      mat[,1] <- rnorm(n, 1.5, 1)                                                   # u[1]
    else # logrw
      mat[,1] <- rnorm(n, 0.2, 1)                                                   # u[1]
    if (fp$eppmod == "logrw"){
      mat[,2:fp$rt$n_rw] <- bayes_rmvt(n, fp$rt$n_rw-1, .epp.env$priors.rw_prior_shape, .epp.env$priors.rw_prior_rate)  # u[2:numKnots]
    } else {
      mat[,2:fp$numKnots] <- bayes_rmvt(n, fp$numKnots-1,.epp.env$priors.tau2_init_shape, .epp.env$priors.tau2_init_rate)  # u[2:numKnots]
    }

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      mat[,fp$numKnots+1] <-  runif (n, r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2])  # ratio r0 / log(iota)
    else
      mat[,fp$numKnots+1] <- sample_iota(n, fp)
  } else if (fp$eppmod == "rlogistic"){
    mat[,1:4] <- t(matrix(rnorm(4*n, .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd), 4))
    mat[,5] <- sample_iota(n, fp)
  } else if (fp$eppmod == "rtrend"){ # r-trend
    mat[,1] <- runif(n, .epp.env$priors.t0.unif.prior[1], .epp.env$priors.t0.unif.prior[2])           # t0
    mat[,2] <- rnorm(n, .epp.env$priors.t1.pr.mean, .epp.env$priors.t1.pr.sd)
    mat[,3] <- rnorm(n, .epp.env$priors.logr0.pr.mean, logr0.pr.sd)  # r0
    mat[,4:7] <- t(matrix(rnorm(4*n, .epp.env$priors.rtrend.beta.pr.mean, .epp.env$priors.rtrend.beta.pr.sd), 4, n))  # beta
  } else if (fp$eppmod == "rhybrid") {
    mat[,1:4] <- t(matrix(rnorm(4*n, .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd), 4))
    mat[,4+1:fp$rt$n_rw] <- rnorm(n*fp$rt$n_rw, 0, rw_prior_sd)  # u[2:numKnots]
    mat[,fp$rt$n_param+1] <- sample_iota(n, fp)
  }

  ## sample ANC model parameters
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0)
    mat[ , ancmod.id] <- sample_prior_ancmod(n, fp$ancmod, fp$prior_args)

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE)
    mat[, fitincrr.id] <- sample_incrr(n, fp)
  
  if (exists("fitpcr", fp))
    mat[, fitpcr.id] <- sexual_rate_sample(n, fp)
  
  return(mat)
}

ldsamp <- function(theta, fp){

  if (exists("prior_args", where = fp)){
    for(i in seq_along(fp$prior_args))
      assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  if (fp$eppmod %in% c("rspline", "logrw")){
    epp_nparam <- fp$numKnots+1
    nk <- fp$numKnots
    if (fp$eppmod == "rspline")  # u[1]
      lpr <- dnorm(theta[1], 1.5, 1, log=TRUE)
    else # logrw
      lpr <- dnorm(theta[1], 0.2, 1, log=TRUE)

    if (fp$eppmod == "logrw")
      bayes_lmvt(theta[2:fp$rt$n_rw], .epp.env$priors.rw_prior_shape, .epp.env$priors.rw_prior_rate)
    else
      lpr <- bayes_lmvt(theta[2:nk], .epp.env$priors.tau2_prior_shape, .epp.env$priors.tau2_prior_rate)

    if (exists("r0logiotaratio", fp) && fp$r0logiotaratio)
      lpr <- lpr + dunif (theta[nk+1], r0logiotaratio.unif.prior[1], r0logiotaratio.unif.prior[2], log=TRUE)
    else
      lpr <- lpr + ldsamp_iota(theta[nk+1], fp)
  } else if (fp$eppmod == "rlogistic") {
    epp_nparam <- 5
    lpr <- sum(dnorm(theta[1:4], .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[5], fp)
  } else if (fp$eppmod == "rtrend"){ # rtrend
    epp_nparam <- 7
    lpr <- dunif (round(theta[1]), .epp.env$priors.t0.unif.prior[1], .epp.env$priors.t0.unif.prior[2], log=TRUE) +
      dnorm(round(theta[2]), .epp.env$priors.t1.pr.mean, .epp.env$priors.t1.pr.sd, log=TRUE) +
      dnorm(theta[3], logr0.pr.mean, logr0.pr.sd, log=TRUE) +
      sum(dnorm(theta[4:7], .epp.env$priors.rtrend.beta.pr.mean, .epp.env$priors.rtrend.beta.pr.sd, log=TRUE))
  } else if (fp$eppmod == "rhybrid"){
    epp_nparam <- fp$rt$n_param+1
    lpr <- sum(dnorm(theta[1:4], .epp.env$priors.rlog_pr_mean, .epp.env$priors.rlog_pr_sd, log=TRUE)) +
      sum(dnorm(theta[4+1:fp$rt$n_rw], 0, rw_prior_sd, log=TRUE))
    lpr <- lpr + ldsamp_iota(theta[fp$rt$n_param+1], fp)
  }

  nparam <- epp_nparam
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
    ancmod.id <- nparam + 1:fp$ancmod$nparam
    nparam <- nparam + fp$ancmod$nparam
    theta_anc <- theta[ancmod.id]
    lpr <- lpr + lprior_ancmod(theta_anc, fp$ancmod, fp$prior_args)
  }

  if(exists("fitincrr", where=fp) && fp$fitincrr==TRUE){
    fitincrr.id <- nparam + 1:getnparam_incrr(fp)
    nparam <- nparam + getnparam_incrr(fp)
    lpr <- lpr + lprior_incrr(theta[cols], fp)
  }
  
  if(exists("fitpcr", where=fp)){
    fitpcr.id <- nparam + 1:8
    nparam <- nparam + 8
    lpr <- lpr + sexual_rate_prior(theta[fitpcr.id], fp)
  }


  return(lpr)
}

prior <- function(theta, fp, log=FALSE){
  if (is.vector(theta))
    lval <- lprior(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (lprior(theta[i,], fp))))
  if (log)
    return(lval)
  else
    return(exp(lval))
}

likelihood <- function(theta, fp, likdat, log=FALSE, doParallel=FALSE) {
  if (is.vector(theta)) {
    lval <- sum(ll_all(theta, fp, likdat))
  } else {
    theta_id <- seq_len(nrow(theta))
    ll_fn    <- function(i) sum(ll_all(theta[i,], fp, likdat))
    if (!.Platform$OS.type=='unix' | !doParallel) {
      lval <- unlist(lapply(theta_id, ll_fn))
    } else {
      n_cores <- max(parallel::detectCores()-2, 1)
      if (is.na(n_cores)) n_cores <- 1
      # cat('calculate ll on', n_cores, 'cores)\n')
      lval <- unlist(parallel::mclapply(theta_id, ll_fn, mc.cores = n_cores))
    }
  }
  if (any(!is.finite(lval))) 
    lval[!is.finite(lval)] <- -1e6
  if (log)
    return(lval)
  else
    return(exp(lval))
}

dsamp <- function(theta, fp, log=FALSE){
  if (is.vector(theta))
    lval <- ldsamp(theta, fp)
  else
    lval <- unlist(lapply(seq_len(nrow(theta)), function(i) (ldsamp(theta[i,], fp))))
  if (any(!is.finite(lval))) 
    lval[!is.finite(lval)] <- -1e6
  if (log)
    return(lval)
  else
    return(exp(lval))
}
