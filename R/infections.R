#' Mixing function as in Tim's model
fn_mixmat_f <- function(j, i, kappa = 2, rho = 0.25, r = 0.15) {
	if (j >= i) {
		num = kappa * rho^kappa * (j - i + r)^(kappa - 1)
		den = (1 + (rho * (j - i + r))^kappa)^2
		return(num/den)
	} 
	else {
		return(0)
	}
}
# Annualized number of new infections
# infect_spec <- function(pop, hivpop, artpop, ii, fp)
# K moved to popClass method
infectFns <- c(
infect_mix = function(hivpop, artpop, ii) {
    ts <- (year-2)/DT + ii
    update_active_pop_to(year)
		
    # calculate and distribute relative infection by CD4 stages
		N1 <- colSums(artpop$data[,,,,year],,2) + colSums(hivpop$data[,,,year])
    N2 <- hivpop$stage0[,,year] + (N1 - hivpop$stage0[,,year]) * p$rel_vl[2]

		N1 <- sapply(1:2, function(x) rep(N1[, x], h.ag.span)) 
		N2 <- sapply(1:2, function(x) rep(N2[, x], h.ag.span))
		
    PP <- data_active[,,hivp.idx]/N1
    PP[which(is.na(PP) | !is.finite(PP), TRUE)] <- 0
		hiv_cd4_adj <- N2 * PP
    
    data_active <<- sweepx(data_active, 1:2, p$est_senesence)
    # sweep over sexual mixing matrices, this results in the number of partnerships
    nc_m <- sweepx(p$mixmat[,,m.idx], 1, p$est_pcr[, 1])
    nc_f <- sweepx(p$mixmat[,,f.idx], 1, p$est_pcr[, 2])
		# get the total numner of partnerships formed by HIV negative population
    nc_m_total <- sweepx(nc_m, 1, rowSums(data_active[, m.idx, ]))
    nc_f_total <- sweepx(nc_f, 1, rowSums(data_active[, f.idx, ]))
		# the balancing ratio
    ratio_mf <- nc_m_total / t(nc_f_total)
    # adjusted number of partnerships
    nc_m_adj <- nc_m * ratio_mf^(-1 + p$balancing)
    nc_f_adj <- nc_f * t(ratio_mf)^p$balancing
    # adjusted number of partnerships
    nc_m_adj <- sweepx(nc_m_adj, 1, data_active[,m.idx,hivn.idx]/rowSums(data_active[,m.idx,]))
    nc_f_adj <- sweepx(nc_f_adj, 1, data_active[,f.idx,hivn.idx]/rowSums(data_active[,f.idx,]))
		# at this point we have the number of partnerships in each age combinations
    art_cov <- matrix(0, pAG, NG)
    if (year >= p$tARTstart) {
      art_ <- colSums(artpop$data[,,,,year] + artpop$data_db[,,,,year],,2)
      hiv_ <- colSums(hivpop$data[,,,year] + hivpop$data_db[,,,year],,1)
      art_cov <- art_/(art_+hiv_)      # when none infected can be NaN
			if(any(is.nan(art_cov)))
				art_cov[which(is.nan( art_cov ), TRUE)] <- 0
      art_cov <- sapply(1:2, function(x) rep(art_cov[, x], h.ag.span))
    }
		transm_prev <- hiv_cd4_adj * (1 - art_cov * p$relinfectART) / rowSums(data_active,,2) # prevalence adjusted for art
    # other one is "transm"
    sex_factor = ifelse(p$incidmod == "eppspectrum", p$incrr_sex[year], p$mf_transm_rr[year])
    if (p$proj.steps[ts] == p$tsEpidemicStart)
      transm_prev <- p$leading_ev * p$iota
		# r(t) x HIV prevalence
    inc_r <- rvec[ts] * sweepx(transm_prev, 2, c(sex_factor, 1))
		# x contact rate adjusted
    inc_m <- sweepx(nc_m_adj, 2, inc_r[, f.idx])
    inc_f <- sweepx(nc_f_adj, 2, inc_r[, m.idx])
		# x reduction using condom in men
    inc_m <- sweepx(inc_m, 1, 1-p$est_condom[, m.idx, year])
    inc_f <- sweepx(inc_f, 2, 1-p$est_condom[, m.idx, year])
		# number of new infections
    infections.ts <- cbind(
      rowSums(sweepx(inc_m, 1, data_active[, m.idx, hivn.idx])), 
      rowSums(sweepx(inc_f, 1, data_active[, f.idx, hivn.idx]))
    )
    
    # incrate15to49_ts[,,ts] <<- inc_r
    prev15to49_ts[ts] <<- prevlast <<- sum(data[,,hivp.idx,year])/sum(data[,,,year])
    infections.ts
},

infect_spec = function(hivpop, artpop, time_step) {
    ts    <- (year-2) / DT + time_step
    dt_ii <- 1 - DT * (time_step - 1) # transition of population    
    update_active_pop_to(year)
  
    # counting all negative including virgin
    hivn_both       <- n_NEG(dt_ii)
    hivp_inactive   <- 0
    if (MODEL==2) # add safe positive 
      hivp_inactive <- VIRGIN$n_HIV(dt_ii)
    hivp_active     <- n_HIV(dt_ii) - hivp_inactive
    all_pop         <- N(dt_ii)
  
    art.ii  <- 0
    if (year >= p$tARTstart)
      art.ii  <- artpop_adj(hivpop, artpop, dt_ii)

    transm_prev <- (hivp_active - art.ii*(1 - p$relinfectART)) / all_pop
    w <- p$iota * (p$proj.steps[ts] == p$tsEpidemicStart)
    inc_rate <- rvec[ts] * transm_prev + w

    sus_age_sex <- data_active[p.age15to49.idx,,hivn.idx]
    adj_sex <- sum(sus_age_sex) /
      ( sum(data_active[p.age15to49.idx,m.idx,hivn.idx]) +
        sum(data_active[p.age15to49.idx,f.idx,hivn.idx]) * 
        p$incrr_sex[year] )
    sexinc15to49.ts <- inc_rate * c(1, p$incrr_sex[year]) * adj_sex
    # New infections distributed by age: ratio age_i/ 25-29 age
    adj_age <- sexinc15to49.ts * colSums(sus_age_sex) /
      colSums(sus_age_sex * p$incrr_age[p.age15to49.idx,,year])
    agesex.inc <- sweep(p$incrr_age[,,year], 2, adj_age, "*")

    ## Adjust age-specific incidence among men for circumcision coverage
    agesex.inc[, m.idx] <- agesex.inc[, m.idx] * (1 - p$circ_incid_rr * p$circ_prop[,year])
    infections.ts <- agesex.inc * data_active[,,hivn.idx]
    # saving
    incrate15to49_ts[ts] <<- inc_rate
    prev15to49_ts[ts] <<- prevlast <<- (hivp_active + hivp_inactive) / all_pop
    infections.ts
},

epp_disease_model_direct = function(hivpop, artpop) {
  if (p$incidpopage) 
    age_id <- p.age15plus.idx # incidence for 15+ population
  else
    age_id <- p.age15to49.idx # incidence for 15 -49 population
  update_active_pop_to(year-1)
  num <- c(1, p$incrr_sex[year]) * sum(data_active[age_id,, hivn.idx])
  den <- sum(data_active[age_id, m.idx, hivn.idx]) + 
         sum(data_active[age_id, f.idx, hivn.idx]) * 
         p$incrr_sex[year]
  sexinc <- p$incidinput[year] * num / den
  ageinc <- colSums( data_active[age_id,,hivn.idx] * 
                     p$incrr_age[age_id,,year] ) /
            colSums(data_active[age_id,,hivn.idx])
  agesex.inc <- sweep(p$incrr_age[,,year], 2, sexinc / ageinc, "*")
  infections[,,year]     <<- agesex.inc * data[,,hivn.idx,year-1]
  data[,,hivn.idx,year]  <<- data[,,hivn.idx,year] - infections[,,year]
  data[,,hivp.idx,year]  <<- data[,,hivp.idx,year] + infections[,,year]
  infect_agrp            <- sumByAGs(infections[,,year], ag.idx)
  hivpop$data[,,,year] <- hivpop$data[,,,year] + 
                          sweep(p$cd4_initdist, 2:3, infect_agrp, "*")
  incid15to49[year]      <<- sum(infections[p.age15to49.idx,,year])
})

setMembers(popEPP, "public", names(infectFns), infectFns)
infectFns <- NULL # optional

calc_infections_simpletransm <- function(fp, pop, hivpop, artpop, i, ii, r_ts){

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  ts <- (i-2)/DT + ii
  
  ## Calculate prevalence of unsuppressed viral load among sexually active population

  contacts_ii <- sweep((pop[c(2:pAG, pAG), , , i] * (1-DT*(ii-1)) + pop[ , , , i] * DT*(ii-1)),
                       1:2, fp$relbehav_age, "*")

  
  hivpop_w_ha <- colSums(sweep(hivpop[ , , , i], 1, fp$relsexact_cd4cat, "*"), , 1)
  hivpop_ha <- colSums(hivpop[ , , , i],,1)
  artpop_ha <- colSums(artpop[ , , , , i],,2)

  hivcontacts_ha <- (hivpop_w_ha + artpop_ha) / (hivpop_ha + artpop_ha)
  hivcontacts_ha[is.na(hivcontacts_ha)] <- 0

  hivtransm_ha <- (hivpop_w_ha + fp$relinfectART * artpop_ha) / (hivpop_ha + artpop_ha)
  hivtransm_ha[is.na(hivtransm_ha)] <- 0

  hivn_ii <- contacts_ii[ , , hivn.idx]
  hivcontacts_ii <- contacts_ii[ , , hivp.idx] * hivcontacts_ha[fp$ss$ag.idx, ]
  hivtransm_ii <- contacts_ii[ , , hivp.idx] * hivtransm_ha[fp$ss$ag.idx, ]

  hivtransm_prev <- colSums(hivtransm_ii) / (colSums(hivn_ii) + colSums(hivcontacts_ii))

  ## r_sex[1:2] is the transmission rate by (Men, Women)
  r_sex <- c(sqrt(fp$mf_transm_rr[i]), 1/sqrt(fp$mf_transm_rr[i])) * r_ts

  sexinc15to49.ts <- (r_sex * hivtransm_prev)[2:1] + fp$mf_transm_rr[i]^c(-0.25, 0.25) * fp$iota * (fp$proj.steps[ts] == fp$tsEpidemicStart)
  agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc15to49.ts/(colSums(pop[p.age15to49.idx,,hivn.idx,i] * fp$incrr_age[p.age15to49.idx,,i])/colSums(pop[p.age15to49.idx,,hivn.idx,i])), "*")

  ## Adjust age-specific incidence among men for circumcision coverage
  agesex.inc[ , m.idx] <- agesex.inc[ , m.idx] * (1 - fp$circ_incid_rr * fp$circ_prop[ , i])
  
  infections.ts <- agesex.inc * pop[,,hivn.idx,i]

  attr(infections.ts, "incrate15to49.ts") <- 0 # sum(infections.ts[p.age15to49.idx,]) / sum(hivn.ii)
  ## attr(infections.ts, "prevcurr") <- sum(hivp_noart.ii+art.ii) / sum(hivn.ii+hivp_noart.ii+art.ii)

  attr(infections.ts, "prevcurr") <- 0 # sum(hivp.ii) / sum(hivp.ii + hivn.ii)

  return(infections.ts)
}
