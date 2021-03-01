#include "Classes.hpp"

void popC::infect_spec (const hivC& hivpop, const artC& artpop, int time_step,
                        Views& v, const Parameters& p, const StateSpace& s) {
  int ts = (s.year-1)/ s.DT + time_step,
      p_lo = s.p_age15to49_[0] - 1, h_lo = s.h_age15to49_[0] - 1;
  double dt_ii = 1 - s.DT * time_step, // transition of population in 1 year
         n_neg_mf = 0, n_pos_mf = 0, n_pos_inactive = 0, n_pos_inactive_lo = 0,
         n_pos_lo = 0, n_pos_up = 0, n_neg_lo = 0, n_neg_up = 0,
         n_hiv_lo = 0, n_art_lo = 0, n_hiv_up = 0, n_art_up = 0, art_ii = 0;

  update_active_pop_to(s.year, v, s); // substract virgin when needed

  for (int sex = 0; sex < s.NG; sex++) {
    for (int age = p_lo; age < s.pAG_1549; age++) {
      n_neg_mf += v.now_pop[s.N][sex][age]; //  this includes debut neg
      if (s.MODEL == 2 && age < s.pDB)
        n_pos_inactive += data_db[s.year][s.P][sex][age]; // "safe" positive
      n_pos_mf += data_active[s.P][sex][age]; // transmissible positive
    }
    n_neg_lo += v.now_pop[s.N][sex][p_lo];
      if (s.MODEL == 2)
        n_pos_inactive_lo += data_db[s.year][s.P][sex][p_lo]; // "safe" positive
    n_neg_up += v.now_pop[s.N][sex][s.pAG_1549];
    n_pos_lo += data_active[s.P][sex][p_lo];
    n_pos_up += data_active[s.P][sex][s.pAG_1549];
    if (s.year >= s.tARTstart-1) {
      for (int agr = h_lo; agr < s.hAG_1549; agr++)
        for (int cd4 = 0; cd4 < s.hDS; cd4++)
          for (int dur = 0; dur < s.hTS; dur++) {
            art_ii   += v.now_art[sex][agr][cd4][dur];
            n_art_lo += v.now_art[sex][h_lo][cd4][dur];
            n_art_up += v.now_art[sex][s.hAG_1549][cd4][dur];
          }
    }
    for (int cd4 = 0; cd4 < s.hDS; cd4++) {
      n_hiv_lo += v.now_hiv[sex][h_lo][cd4];
      n_hiv_up += v.now_hiv[sex][s.hAG_1549][cd4];
    }
  }
  
  double hivp_inactive = n_pos_inactive - n_pos_inactive_lo * dt_ii;
  double hivn_both = n_neg_mf - n_neg_lo * dt_ii + n_neg_up * dt_ii;
  double hivp_active = n_pos_mf - n_pos_lo * dt_ii + n_pos_up * dt_ii;
  
  if (s.year >= s.tARTstart-1) {
    if (n_hiv_lo + n_art_lo > 0) {
      for (int sex = 0; sex < s.NG; sex++) {
        double art_trans = 0, hiv_trans = 0;
        for (int cd4 = 0; cd4 < s.hDS; cd4++) {
          for (int dur = 0; dur < s.hTS; dur++)
            art_trans += v.now_art[sex][h_lo][cd4][dur];
          hiv_trans += v.now_hiv[sex][h_lo][cd4];
        }
        art_ii -= ( data_active[s.P][sex][p_lo] * 
                    art_trans / (hiv_trans + art_trans) ) * dt_ii;
      }
    }
    if (n_hiv_up + n_art_up > 0) {
      for (int sex = 0; sex < s.NG; sex++) {
        double art_trans = 0, hiv_trans = 0;
        for (int cd4 = 0; cd4 < s.hDS; cd4++) {
          for (int dur = 0; dur < s.hTS; dur++)
            art_trans += v.now_art[sex][s.hAG_1549][cd4][dur];
          hiv_trans += v.now_hiv[sex][s.hAG_1549][cd4];
        }
        art_ii += ( data_active[s.P][sex][s.pAG_1549] * 
                    art_trans / (hiv_trans + art_trans) ) * dt_ii;
      }
    }
  }
  
  // Prob of contacting a sexual active, H+ in total pop
  double
  transm_prev = (hivp_active - art_ii * (1 - p.ic.relinfectART)) / 
                (hivn_both + hivp_active + hivp_inactive);

  double w = (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) ? p.ic.iota : 0.0;
  double inc_rate = rvec[ts] * transm_prev + w;

  // Incidence: male = negative / female negative * sexratio + male negative; 
  //          female = male * sexratio
  double n_neg_m = 0, n_neg_f = 0;
  for (int age = p_lo; age < s.pAG_1549; age++) {
    n_neg_m += data_active[s.N][s.M][age];
    n_neg_f += data_active[s.N][s.F][age];
  }
  double adj_sex = (n_neg_m + n_neg_f) / (n_neg_m + n_neg_f * p.ic.incrr_sex[s.year]);
  double sex_inc[2] = {inc_rate * adj_sex, inc_rate * adj_sex * p.ic.incrr_sex[s.year]};
  // New infections distributed by age: ratio age_i/ 25-29 age
  for (int sex = 0; sex < s.NG; sex++) {
    double n_neg = 0, n_neg_rr = 0, adj_age;
    for (int age = p_lo; age < s.pAG_1549; age++) {
      n_neg += data_active[s.N][sex][age]; 
      n_neg_rr += p.ic.incrr_age[s.year][sex][age] * data_active[s.N][sex][age];
    }
    adj_age = sex_inc[sex] / ( n_neg_rr / n_neg );
    for (int age = 0; age < s.pAG; age++) {
      if (sex == s.M) // age-specific incidence among circumcised men
        adj_age *= (1 - p.ic.circ_incid_rr * p.ic.circ_prop[s.year][age]);
      infections_[sex][age] = p.ic.incrr_age[s.year][sex][age] * adj_age * 
        data_active[s.N][sex][age];
    }
  }

  // saving
  incrate15to49_ts[ts] = inc_rate;
  prev_last = (hivp_active + hivp_inactive) / (hivn_both + hivp_active + hivp_inactive);
  prev15to49_ts[ts] = prev_last;
}

void popC::infect_mix (hivC& hivpop, artC& artpop, int ii, Views& v, const Parameters& p, const StateSpace& s) {
  update_active_pop_to(s.year, v, s);
	
	boost2D 
		N1(extents[s.NG][s.hAG]),
		N2(extents[s.NG][s.hAG]);

  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++) {
				double tmp = v.now_hiv[sex][agr][cd4];
        for (int dur = 0; dur < s.hTS; dur++)
					tmp += v.now_art[sex][agr][cd4][dur];
				N1[sex][agr] += tmp;
				N2[sex][agr] += tmp * p.ic.rel_vl[cd4];
			}

	boost2D 
		N11(extents[s.NG][s.pAG]),
		N21(extents[s.NG][s.pAG]),
		hiv_cd4_adj(extents[s.NG][s.pAG]);

	for (int sex = 0; sex < s.NG; sex++) {
		int cusu = 0;
		for (int agr = 0; agr < s.hAG; agr++) {
			for (int id = 0; id < s.h_ag_span[agr]; id++) {
				int age = cusu + id;
				N11[sex][age] = (N1[sex][agr] == 0) ? 1 : N1[sex][agr];
				N21[sex][age] = N2[sex][agr];
				hiv_cd4_adj[sex][age] = N21[sex][age] * data_active[s.P][sex][age]/N11[sex][age];
			}
			cusu += s.h_ag_span[agr];
		}
	}

  for (int ds = 0; ds < s.pDS; ds++)
    for (int sex = 0; sex < s.NG; sex++)
      for (int age = 0; age < s.pAG; age++)
        data_active[ds][sex][age] *= p.ic.est_senesence[sex][age];

  // balancing number of sex acts
  boost2D
    nc_m(extents[s.pAG][s.pAG]), nc_m_total(extents[s.pAG][s.pAG]),
    nc_f(extents[s.pAG][s.pAG]), nc_f_total(extents[s.pAG][s.pAG]),
    nc_m_adj(extents[s.pAG][s.pAG]),
    nc_f_adj(extents[s.pAG][s.pAG]);
  /* number of partnerships */
  for (int r = 0; r < s.pAG; ++r) { // over my ages
    for (int c = 0; c < s.pAG; ++c) { // over partner ages
      nc_m[c][r] = p.ic.mixmat[s.M][c][r] * (p.ic.incrr_age[s.year][s.M][r]);
      nc_f[c][r] = p.ic.mixmat[s.F][c][r] * (p.ic.incrr_age[s.year][s.F][r]);
    }
  }
	/* number of pns formed by total pop */
  for (int r = 0; r < s.pAG; ++r) { // over my ages
    for (int c = 0; c < s.pAG; ++c) { // over partner ages
      nc_m_total[c][r] = nc_m[c][r] * (data_active[s.N][s.M][r] + data_active[s.P][s.M][r]);
      nc_f_total[c][r] = nc_f[c][r] * (data_active[s.N][s.F][r] + data_active[s.P][s.F][r]);
    }
  }
	/* balancing ratio */
  for (int r = 0; r < s.pAG; ++r) 
    for (int c = 0; c < s.pAG; ++c) {
      double ratio_mf = nc_m_total[c][r] / nc_f_total[r][c];
      nc_m_adj[c][r] = nc_m[c][r] / pow(ratio_mf, 0.5);
      nc_f_adj[r][c] = nc_f[r][c] * pow(ratio_mf, 0.5);
    }
	/* on negative only  */
  for (int r = 0; r < s.pAG; ++r) 
    for (int c = 0; c < s.pAG; ++c) {
      nc_m_adj[c][r] *= (data_active[s.N][s.M][r] / (data_active[s.N][s.M][r] + data_active[s.P][s.M][r]));
      nc_f_adj[c][r] *= (data_active[s.N][s.F][r] / (data_active[s.N][s.F][r] + data_active[s.P][s.F][r]));
    }

  boost2D art_cov(extents[s.NG][s.pAG]);
  if (s.year >= s.tARTstart-1)
    art_cov = age_sex_cov(hivpop, artpop, v, p, s);

  int ts = (s.year-1)/s.DT + ii;

  boost2D transm_prev(extents[s.NG][s.pAG]);
  for (int sex = 0; sex < s.NG; sex++) {
    for (int age = 0; age < s.pAG; age++) {
      double all_pop = data_active[s.N][sex][age] + data_active[s.P][sex][age];
      transm_prev[sex][age] = (hiv_cd4_adj[sex][age] * (1 - art_cov[sex][age] * p.ic.relinfectART)) / all_pop;
    }
  }
  // only two modes: "eppspectrum"=0 or "transm"=1
  double sex_factor = (p.ic.incidmod == 0) ? p.ic.incrr_sex[s.year] : p.ic.mf_transm_rr[s.year];
  //+intervention effects and time epidemic start
  if (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) {
    for (int age = 0; age < s.pAG; ++age) {
      transm_prev[s.M][age] += p.ic.iota;
      transm_prev[s.F][age] += p.ic.iota * pow(sex_factor, 0.5);
    }    
  } 
	/* the incrr is acting as transm on the prevalence (column) not susceptible (row) */
	for (int age = 0; age < s.pAG; ++age) {
		transm_prev[s.M][age] *= p.ic.incrr_age[s.year][s.M][age];
		transm_prev[s.F][age] *= p.ic.incrr_age[s.year][s.F][age];
	}    

  multiply_with_inplace(transm_prev, rvec[ts]);
  for (int age = 0; age < s.pAG; ++age)
    transm_prev[s.M][age] *= sex_factor;

  boost2D inc_m(extents[s.pAG][s.pAG]), inc_f(extents[s.pAG][s.pAG]);

  // adjusted to IRRa
  for (int r = 0; r < s.pAG; ++r)
    for (int c = 0; c < s.pAG; ++c) {
      inc_m[c][r] = data_active[s.N][s.M][r] * nc_m_adj[c][r] * transm_prev[s.F][c] * 
        (1 - p.ic.est_condom[s.year][s.M][r]);
      inc_f[c][r] = data_active[s.N][s.F][r] * nc_f_adj[c][r] * transm_prev[s.M][c] * 
        (1 - p.ic.est_condom[s.year][s.M][c]); // male driver condom effect
    }

  boost1D inc_mv = rowSums(inc_m), inc_fv = rowSums(inc_f);

  for (int age = 0; age < s.pAG; ++age) {
    infections_[s.M][age] = inc_mv[age];
    infections_[s.F][age] = inc_fv[age];
  }

  if (p.ic.proj_steps[ts] == p.ic.tsEpidemicStart) {
    double ob = sumArray(infections_);
    double ex = p.ic.iota * sumArray(v.now_pop);
    multiply_with_inplace(infections_, ex/ob);
  }

  // prev15to49_ts_m should use this one! now just store as below
  boost2D n_pos = v.now_pop[ indices[s.P][_all][_all] ];
  prev15to49_ts[ts] = sumArray(n_pos) / sumArray(v.now_pop);
  prev_last = prev15to49_ts[ts];
}

void popC::epp_disease_model_direct(hivC& hivpop, artC& artpop, Views& v,
                                    const Parameters& p, const StateSpace& s) {
  int a_l, a_r;
  if (p.ic.incidpopage) { // incidence for 15+ population
    a_l = s.p_age15plus_[0] - 1;
    a_r = s.pAG_15plus;
  } else { // incidence for 15 -49 population
    a_l = s.p_age15to49_[0] - 1;
    a_r = s.pAG_1549;
  }
  update_active_last_year(v, s);
  double n_m = 0, n_f = 0;
  for (int age = a_l; age < a_r; ++age) {
    n_m += active_last_year_[s.N][s.M][age];
    n_f += active_last_year_[s.N][s.F][age];
  }
  dvec sex_inc(s.NG);
  sex_inc[s.M] = (n_m + n_f) * p.ic.incidinput[s.year] / 
                     (n_m + n_f  * p.ic.incrr_sex[s.year]);
  sex_inc[s.F] = (n_m + n_f) * p.ic.incidinput[s.year] * p.ic.incrr_sex[s.year] /
                     (n_m + n_f  * p.ic.incrr_sex[s.year]);
  dvec ageinc(s.NG);
  for (int sex = 0; sex < s.NG; sex++) {
    double neg_sa = 0, inc_sa = 0;
    for (int age = a_l; age < a_r; age++) {
      neg_sa += active_last_year_[s.N][sex][age];
      inc_sa += active_last_year_[s.N][sex][age] * p.ic.incrr_age[s.year][sex][age];
    }
    ageinc[sex] = inc_sa / neg_sa;
  }
  double new_infect = 0;
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = 0; age < s.pAG; age++) {
      new_infect =
        p.ic.incrr_age[s.year][sex][age] * ( sex_inc[sex] / ageinc[sex]) * 
        active_last_year_[s.N][sex][age];
      infections[s.year][sex][age]        = new_infect;
      v.now_pop[s.N][sex][age] -= new_infect;
      v.now_pop[s.P][sex][age] += new_infect;
    }
  boost2D infect_agrp = 
    sumByAG(infections[ indices[s.year][_all][_all]], s.ag_, s.hAG);
  for (int sex = 0; sex < s.NG; sex++)
    for (int agr = 0; agr < s.hAG; agr++)
      for (int cd4 = 0; cd4 < s.hDS; cd4++)
        v.now_hiv[sex][agr][cd4] +=
          p.nh.cd4_initdist[sex][agr][cd4] * infect_agrp[sex][agr];
  for (int sex = 0; sex < s.NG; sex++)
    for (int age = s.p_age15to49_[0] - 1; age < s.pAG_1549; age++)
      incid15to49[s.year] += infections[s.year][sex][age];
}
