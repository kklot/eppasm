#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom magrittr %$%
#' @export
magrittr::`%$%`

#' @importFrom magrittr %<>%
#' @export
magrittr::`%<>%`

#' Recursively convert list embeded to numeric
#' 
#' Recursively convert list embeded to numeric
#' 
#' @param a_list a list
#' @return a_list with all element converted to numeric
#' @export
recursively_double <- function(a_list) {
    rapply(a_list, 
        function(x) {
            atx = attributes(x)
            x = as.numeric(x)
            attributes(x) = atx
            x
        }, 
        c("integer", 'logical', 'matrix', 'array'), 
        how='replace'
    )
}

#' Model's attributes to list
#' 
#' We do this in R because naming output in C++ is not as flexible as assigning
#' attribbutes when different models return different set of outputs. 
#' See Rinlinedfuns.h's mkNamed
#' 
#' @param x x
#' @return a list with name
.asList <- function(mod) {
    att = attributes(mod)
    keep = which(names(att) %in% c('dim', 'class'))
    att = att[-keep]
    attributes(mod)[-keep] = NULL
    att[['data']] = mod
    class(att) = class(mod)
    att
}

sweepx <- function(...) sweep(..., FUN='*') # missing FUN too many times!

# Natural age to index
a2i <- function(x, min=15, max=80) which(min:max %in% x)

# Non number to value
na2num <- function(x, y) {x[is.na(x)] <- y; return(x)}

# hazard function for log-logistic distribution (parameterize as in INLA)
hz_llogis <- function (x, alpha = 8, lambda = 1/18) {
    num <- alpha * lambda
    den <- (lambda * x)^(1-alpha) + lambda * x
    num/den
}

mu_llogis <- function(alpha, lambda) { # convert to other
    (1 / lambda * pi * 1 / alpha) /(sin(pi / alpha))
}

mu_llogis2 <- function(alpha, lambda) { # convert to other
    pi / (alpha*lambda* sin(pi / alpha))
}

# density function for log-logistic distribution (parameterize as in INLA)
d_llogis <- function (x, alpha = 8, lambda = 1/18) {
    num <- alpha
    den <- x^(1+alpha)*lambda^alpha + x^(1-alpha)*lambda^(-alpha) + 2*x
    num/den
}
d_llogis2 <- function (x, alpha = 8, lambda = 1/18) {
    num <- alpha * lambda * (x * lambda)^(alpha-1)
    den <- ( (lambda*x)^alpha + 1 )^2
    num/den
}

# Cumulative function for log-logistic distribution (parameterize as in INLA)
cdf_llogis <- function (x, alpha = 8, lambda = 1/18) {
    1 / ( 1 + (lambda * x)^(-alpha) )
}

logit <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))
ldinvlogit <- function(x){v <- invlogit(x); log(v) + log(1-v)}

#' Replace R naming of fp for C++ 
#' 
#' dots are replaced with underscore, char is replaced with integer, 
#' set default value for mixing model and dummy variables
#' 
#' @param fp Fix parameters
prepare_fp_for_Cpp <- function(fp) {
    # for NGM calculation
    fp$cd4_prog_xp = abind(
        fp$cd4_prog[,,1] %>% apply(1, rep, times=fp$ss$h.ag.span),
        fp$cd4_prog[,,2] %>% apply(1, rep, times=fp$ss$h.ag.span), along = 3) 
    fp$cd4_initdist_xp = abind(
        fp$cd4_initdist[,,1] %>% apply(1, rep, times=fp$ss$h.ag.span),
        fp$cd4_initdist[,,2] %>% apply(1, rep, times=fp$ss$h.ag.span), along = 3)
    names(fp$ss) <- gsub('\\.', '_', names(fp$ss))
    names(fp) <- gsub('\\.', '_', names(fp))
    if (exists("rt", where=fp))
        names(fp$rt) <- gsub('\\.', '_', names(fp$rt))
    recursively_double(fp)
}
# Converting prior assumption to parameter boundary for DE
prior_to_DE_bounds <- function(fp) {

  if (exists("prior_args", where = fp)){
	for(i in seq_along(fp$prior_args))
	  assign(names(fp$prior_args)[i], fp$prior_args[[i]])
  }

  up <- c()
  lo <- c()

  if (fp$eppmod == "rhybrid") {
	up <- rlog_pr_mean + 2.58*rlog_pr_sd
	lo <- rlog_pr_mean - 2.58*rlog_pr_sd
	up <- c(up, rep(+2.58*rw_prior_sd, fp$rt$n_rw))
	lo <- c(lo, rep(-2.58*rw_prior_sd, fp$rt$n_rw))
	if (exists("logitiota", fp) && fp$logitiota) {
		up <- c(up, logit(.9999))
		lo <- c(lo, logit(.0001))
	} else {
		up <- c(up, logiota.unif.prior[2])
		lo <- c(lo, logiota.unif.prior[1])
	}
  }
  ## sample ANC model parameters
  if (exists("ancmod", fp) && fp$ancmod$nparam > 0) {
	  if(fp$ancmod$fit_ancbias) {
		up <- c(up, ancbias.pr.mean + 2.58 * ancbias.pr.sd)
		lo <- c(lo, ancbias.pr.mean - 2.58 * ancbias.pr.sd)
	  }
	  if(fp$ancmod$fit_vinfl){
		bs <- range(log(rexp(1000, ancrtcens.vinfl.pr.rate)))
		lo <- c(lo, bs[1])
		up <- c(up, bs[2])
	  }
	  if(fp$ancmod$fit_logfrr){
		lo <- c(lo, log_frr_adjust.pr.mean - 2.58 * log_frr_adjust.pr.sd)
		up <- c(up, log_frr_adjust.pr.mean + 2.58 * log_frr_adjust.pr.sd)
	  }
	  if(fp$ancmod$fit_ancrtcens_vinfl){
		bs <- range(log(rexp(1000, ancrtcens.vinfl.pr.rate)))
		lo <- c(lo, bs[1])
		up <- c(up, bs[2])
	  }
	  if(fp$ancmod$fit_ancrtsite_beta){
	  	lo <- c(lo, ancrtsite.beta.pr.mean - 2.58 * ancrtsite.beta.pr.sd)
		up <- c(up, ancrtsite.beta.pr.mean + 2.58 * ancrtsite.beta.pr.sd)
	  } 

      if (fp$ss$MIX) { # balancing parameter ~ beta(2,2)
            lo <- c(lo, 0)
            up <- c(up, 1)
      }

    return(cbind(lo, up))
    }
}

# -----------------------------------------------------------------------------
# Extracting bits from eppasm.R
# -----------------------------------------------------------------------------

# pop x by age groups
# When x is a vector
#' @importFrom fastmatch ctapply
sumByAG <- function(x, ag.idx, fertile=FALSE, p.fert.idx=NULL) {
  if (!fertile) fastmatch::ctapply(x, ag.idx, sum) 
    else fastmatch::ctapply(x, ag.idx[p.fert.idx], sum)
}

# pop x by age groups
# When x is a data.frame/matrix
sumByAGs <- function(k, ag.idx, fertile=FALSE, p.fert.idx=NULL) 
  apply(k, 2, sumByAG, ag.idx=ag.idx, fertile=fertile, p.fert.idx=p.fert.idx)

# Scale mortality cd4
scale_cd4_mort <- function(hivpop, artpop) {
  year  <- hivpop$year
  if (hivpop$p$scale_cd4_mort && (year >= hivpop$p$tARTstart) ) {
    num   <- hivpop$get(year)
    den   <- artpop$get(year)
    if (hivpop$MODEL == 2) {
      num = num + hivpop$data_db[,,,year]
      den = den + artpop$data_db[,,,,year]
    }
    den   <- colSums(den)
    cd4mx <- num / (num + den)
    cd4mx[!is.finite(cd4mx)] <- 1.0
    cd4_mort_ts <- cd4mx * hivpop$p$cd4_mort
    return(cd4_mort_ts)
  } 
  else
    return(hivpop$p$cd4_mort)
}

# hiv deaths at ts
calc.agdist <- function(x, ag.idx, h.ag.span) {
  d <- x / rep(sumByAG(x, ag.idx), h.ag.span) # percent of each age/age group
  d[is.na(d)] <- 0 
  d
}

# End bits from eppasm.R
# Continue in progression.R
# -----------------------------------------------------------------------------