#include "fpClass.hpp"
#include "hivClass.hpp"
#include "artClass.hpp"
#pragma once

// Pop class
template<typename T>
class popC {
public: // Pop fields
  StateSpace* const s;
  Parameters* const p;
  int N;
  Eigen::array<long, 3> dims;
  epp::array4D<T> data;
  epp::array4D<T> data_db; // debut only population
  epp::map_of_cube<T> data_curr;
  epp::map_of_cube<T> data_prev;
  epp::cube<T>    data_active;
  epp::cube<T>    active_last_year_;
  epp::vector<T>  birth_age;
  epp::vector<T>  birth_agrp;
  epp::vector<T>  prev15to49;  
  epp::vector<T>  incid15to49;
  epp::vector<T>  entrantprev;
  epp::vector<T>  pregprevlag;
  epp::vector<T>  pregprev;
  epp::vector<T>  incrate15to49_ts;  
  epp::vector<T>  prev15to49_ts;
  epp::vector<T>  rvec;
  epp::matrix<T>  birthslag;
  epp::matrix<T>  hivp_entrants_out;
  epp::matrix<T>  hiv_sx_prob;
  epp::matrix<T>  hiv_mr_prob;
  epp::matrix<T>  adj_prob;
  epp::cube<T>    infections;
  epp::cube<T>    hivdeaths;
  epp::cube<T>    natdeaths;
  epp::cube<T>    popadjust;
  epp::vector<T>  artcov; // initially no one on treatment
  T prev_last     = .0; // last time step prevalence
  epp::matrix<T>  infections_;
  epp::cube<T>    art_elig_;
  epp::cube<T>    art_init_;
  epp::matrix<T>  hiv_aging_prob_;
  epp::vector<T>  entrant_art_;
  epp::vector<T>  elig_art_;
  epp::vector<T>  n_art_ii_;
  epp::vector<T>  n_art_init_;
public: // Pop inits
  popC(StateSpace* const ss, Parameters* const pp) :
    s                  (ss), // save a copy of pointer to access
    p                  (pp), // save a copy of pointer to access
    N                  (s->pAG * s->NG * s->pDS),
    dims               ({s->pAG, s->NG, s->pDS}),
    data               (s->pAG, s->NG, s->pDS, s->PROJ_YEARS),
    data_db            (s->pDB, s->NG, s->pDS, s->PROJ_YEARS),
    data_curr          (data.data(), s->pAG, s->NG, s->pDS),
    data_prev          (data.data(), s->pAG, s->NG, s->pDS),
    data_active        (s->pAG, s->NG, s->pDS), // 1 year only
    active_last_year_  (s->pAG, s->NG, s->pDS), // 1 year only
    birth_age          (s->pAG_FERT),
    birth_agrp         (s->hAG_FERT),
    prev15to49         (s->PROJ_YEARS),
    incid15to49        (s->PROJ_YEARS),
    entrantprev        (s->PROJ_YEARS),
    pregprevlag        (s->PROJ_YEARS),
    pregprev           (s->PROJ_YEARS),
    incrate15to49_ts   (s->n_steps),
    prev15to49_ts      (s->n_steps),
    rvec               (s->n_steps),
    birthslag          (s->NG,  s->PROJ_YEARS),
    hivp_entrants_out  (s->NG,  s->PROJ_YEARS),
    hiv_sx_prob        (s->hAG, s->NG),
    hiv_mr_prob        (s->hAG, s->NG),
    adj_prob           (s->hAG, s->NG),
    infections         (s->pAG, s->NG, s->PROJ_YEARS),
    hivdeaths          (s->pAG, s->NG, s->PROJ_YEARS),
    natdeaths          (s->pAG, s->NG, s->PROJ_YEARS),
    popadjust          (s->pAG, s->NG, s->PROJ_YEARS),
    artcov             (s->NG), // 1 year only
    infections_        (s->pAG, s->NG), // avoid repeated allocations
    art_elig_          (s->hDS, s->hAG, s->NG),
    art_init_          (s->hDS, s->hAG, s->NG),
    hiv_aging_prob_    (s->hAG, s->NG),
    entrant_art_       (s->NG * s->pDS),
    elig_art_          (s->hDS),
    n_art_ii_          (s->NG),
    n_art_init_        (s->NG)
  {
    data.setZero();
    data_db.setZero();
    prev15to49.setZero();
    incid15to49.setZero();
    entrantprev.setZero();
    pregprevlag.setZero();
    pregprev.setZero();
    incrate15to49_ts.setZero();
    prev15to49_ts.setZero();
    rvec.setZero();
    if (s->MODEL == 2 && s->pDB == 1)
      Rf_error("Debut model state-space not exist, see update_fp_debut()");

    // Initialize base pop
		epp::map_of_cube<T> now(data.data(), s->pAG, s->NG, s->pDS);
    now.chip(s->N, 2) = p->dm.basepop;
    birthslag = p->dm.birthslag;
    if (p->ic.eppmod != 1) // existing rvec can be shorted than model dimension!
      for (int i = 0; i < p->ic.rvec.size(); ++i)
        rvec(i) = p->ic.rvec(i);

    if (s->MODEL==2)
      for (int sex = 0; sex < s->NG; sex++)
        for (int age = 0; age < s->pDB; age++)
          data_db(age, sex, s->N, 0) = p->dm.basepop(age, sex) * (1 - p->ic.db_rate(age, sex, 0));
  }

  void update_views(int t) { // #TutorialMapPlacementNew
    new (&data_curr) epp::map_of_cube<T>(&data(0) + N * t, dims);
    new (&data_prev) epp::map_of_cube<T>(&data(0) + N * (t - 1), dims);
  }

  void update_active_pop_to (int when) {
    data_active = data_curr;
    if (s->MODEL == 2) {
      epp::index<3> origin = {0,0,0}, xtents = {s->pDB,s->NG,s->pDS};
      data_active.slice(origin, xtents) -= data_db.chip(when, 3);
    }
  }

  void update_active_last_year() {
    active_last_year_ = data_prev;
    if (s->MODEL == 2) {
      epp::index<3> origin = {0,0,0}, xtents = {s->pDB,s->NG,s->pDS};
      active_last_year_.slice(origin, xtents) -= data_db.chip(s->year-1, 3);
    }
  }

  void aging() { // open ended
    epp::index<3>
      origin = {0,0,0}, origin_up = {1,0,0}, xtents = {s->pAG-1, s->NG, s->pDS};
    data_curr.slice(origin_up, xtents)  = data_prev.slice(origin, xtents);
    data_curr.chip(s->pAG-1, 0)         += data_prev.chip(s->pAG-1, 0);
    if (s->MODEL==2) { // sexual debut will move the last age
      epp::index<3> xtents_db = {s->pDB-1, s->NG, s->pDS};
      data_db.chip(s->year, 3).slice(origin_up, xtents_db) = 
        data_db.chip(s->year-1, 3).slice(origin, xtents_db);
    }
  }

  void add_entrants() 
  { // Add lagged births into youngest age group
    epp::vector<T> healthy(s->NG), positiv(s->NG);
    if (p->dm.flag_popadjust) {
      if (s->MODEL==0) {
        healthy = p->dm.entrantpop.chip(s->year-1, 1);
      } else {
        healthy = p->dm.entrantpop.chip(s->year-1, 1) * (1-p->ph.entrantprev.chip(s->year, 1));
        positiv = p->dm.entrantpop.chip(s->year-1, 1) *    p->ph.entrantprev.chip(s->year, 1);
      }
    } else {
      if (s->MODEL==0)
        healthy = birthslag.chip(s->year-1, 1) * p->dm.cumsurv.chip(s->year-1, 1) / 
                  p->ph.paedsurv_lag[s->year-1] + p->dm.cumnetmigr.chip(s->year-1, 1);
      else {
        healthy = birthslag.chip(s->year-1, 1) * p->dm.cumsurv.chip(s->year-1, 1) * 
          (1 - p->ph.entrantprev.chip(s->year, 1) / p->ph.paedsurv_lag[s->year-1]) +
          p->dm.cumnetmigr.chip(s->year-1, 1) * (1 - pregprevlag(s->year-1) * p->ph.netmig_hivprob);
        positiv = (birthslag.chip(s->year-1, 1) * p->dm.cumsurv.chip(s->year-1, 1) + 
          p->dm.cumnetmigr.chip(s->year-1, 1)) * p->ph.entrantprev.chip(s->year, 1);
      }
    }
    // save and update pop
    if (s->MODEL != 0) {
      data_curr.chip(s->N, 2).chip(0, 0) = healthy;
      data_curr.chip(s->P, 2).chip(0, 0) = positiv;
      hivp_entrants_out.chip(s->year, 1) = positiv;
      for (int sex = 0; sex < s->NG; ++sex) { // 1:2 ART+, 3:4 ART-
        entrant_art_(sex)   = positiv(sex) *      p->ph.entrantartcov(sex, s->year);
        entrant_art_(sex+2) = positiv(sex) * (1 - p->ph.entrantartcov(sex, s->year));
      }
      entrantprev(s->year) = sum_vector(positiv) / (sum_vector(positiv) + sum_vector(healthy));
    }
    if (s->MODEL==2) { // add to virgin to record
      data_db.chip(s->year, 3).chip(s->N, 2).chip(0, 0) = healthy;
      data_db.chip(s->year, 3).chip(s->P, 2).chip(0, 0) = positiv;
    }
  }

  void sexual_debut () {
    epp::matrix<T> _db = p->ic.db_rate.chip(s->year, 2);
    data_db.chip(s->year, 3) *= (1 - expand_right_2to3(_db, s->pDS));
  }

  void update_hiv_aging_prob () {
    epp::matrix<T> hiv_by_agrp_ = sumByAG2D<T>(data_prev.chip(s->P, 2), s->ag_, s->hAG);
    for (int sex = 0; sex < s->NG; sex++)
      for (int agr = 0; agr < s->hAG; agr++)
        hiv_aging_prob_(agr, sex) = (hiv_by_agrp_(agr, sex) == 0) ? 0 :
          data_prev(s->aglast_(agr) - 1, sex, s->P) / hiv_by_agrp_(agr, sex);
  }

  void deaths() {
    epp::matrix<T> 
      hiv_by_agrp_ = sumByAG2D<T>(data_curr.chip(s->P, 2), s->ag_, s->hAG),
      Sx           = p->dm.Sx.chip(s->year, 2),
      Sx_db        = Sx.slice(epp::index<2>({0,0}), epp::index<2>({s->pDB,s->NG}));
    epp::cube<T> num_death_ = data_curr * (1 - expand_right_2to3(Sx, s->pDS));
    data_curr *= expand_right_2to3(Sx, s->pDS);
    data_db.chip(s->year, 3) *= expand_right_2to3(Sx_db, s->pDS);
    // calculate survival prob for hivpop and artpop
    auto death_by_agrp_ = sumByAG2D<T>(num_death_.chip(s->P, 2), s->ag_, s->hAG);
    hiv_sx_prob    = 1 - death_by_agrp_ / hiv_by_agrp_;
    replace_na_with(hiv_sx_prob);
    //  save natural death outputs
    natdeaths.chip(s->year, 2) = num_death_.sum(epp::by_rack); // over DS
  }

  void migration () {
    epp::matrix<T> 
      netmigsurv       = p->dm.netmigr.chip(s->year, 2) * (1 + p->dm.Sx.chip(s->year, 2)) / 2.,
      migrate_prob_    = 1 + netmigsurv / data_curr.sum(epp::by_rack),
      num_migrate_     = migrate_prob_ * data_curr.chip(s->P, 2), 
      migrant_by_agrp_ = sumByAG2D<T>(num_migrate_, s->ag_, s->hAG), 
      hiv_by_agrp_     = sumByAG2D<T>(data_curr.chip(s->P, 2), s->ag_, s->hAG);
    hiv_mr_prob      = migrant_by_agrp_ / hiv_by_agrp_;
    replace_na_with(hiv_mr_prob);
    data_curr *= expand_right_2to3(migrate_prob_, s->pDS);
    if (s->MODEL == 2) {
      epp::index<2> db_start = {0,0}, db_xtend = {s->pDB, s->NG};
      epp::matrix<T> mg_db = migrate_prob_.slice(db_start, db_xtend);
      data_db.chip(s->year, 3) *= expand_right_2to3(mg_db, s->pDS);
    }
  }

  void update_fertile () // only on active pop
  {
    update_active_pop_to(s->year);
    update_active_last_year();

    epp::index<1> 
      sta_fert = {int(s->p_fert_(0) - 1)}, 
      len_fert = {int(s->pAG_FERT - (s->p_fert_(0) - 1))};

    auto f_now  = data_active.chip(s->F, 1).sum(epp::by_cols);
    auto f_pre  = active_last_year_.chip(s->F, 1).sum(epp::by_cols);
    auto _birth = (f_now + f_pre) / 2.;
    birth_age   = _birth.slice(sta_fert, len_fert) * p->dm.asfr.chip(s->year, 1);

    if (s->MODEL==2) { // adjust ASFR
      auto f_now2  = data_curr.chip(s->F, 1).sum(epp::by_cols);
      auto f_pre2  = data_prev.chip(s->F, 1).sum(epp::by_cols);
      auto N_mid   = (f_now2 + f_pre2) / 2.; // 66
      birth_age   *= p->dm.asfr.chip(s->year, 1) / (birth_age / N_mid.slice(sta_fert, len_fert));
    }
    epp::vector<T> sub_id = s->ag_.slice(sta_fert, len_fert);
    birth_agrp = sumByAG1D<T>(birth_age, sub_id, s->hAG_FERT);
    if ( (s->year + s->AGE_START) <= (s->PROJ_YEARS - 1) )
      birthslag.chip(s->year + s->AGE_START - 1, 1) = p->dm.srb.chip(s->year, 1) * sumArray(birth_agrp);
  }

  void adjust_pop () {
    popadjust.chip(s->year, 2) = p->dm.targetpop.chip(s->year, 2) / data_curr.sum(epp::by_rack);
    if (s->MODEL!=0) {  // calculate asjust prob for hiv and art pops
      epp::matrix<T> 
        num_adjust_         = popadjust.chip(s->year, 2) * data_curr.chip(s->P, 2),
        num_adjust_by_agrp_ = sumByAG2D<T>(num_adjust_, s->ag_, s->hAG), 
        hiv_by_agrp_        = sumByAG2D<T>(data_curr.chip(s->P, 2), s->ag_, s->hAG);
      adj_prob            = num_adjust_by_agrp_ / hiv_by_agrp_;
      replace_na_with(adj_prob);
    }
    // only then adjust myself
    epp::cube<T> adjxp = expand_right_2to3<T>(popadjust.chip(s->year, 2), s->pDS);
    data_curr *= adjxp;
    if (s->MODEL == 2)
      for (int ds = 0; ds < s->pDS; ds++)  
        for (int sex = 0; sex < s->NG; sex++)
          for (int age = 0; age < s->pDB; age++)
            data_db(age, sex, ds, s->year) *= popadjust(age, sex, s->year);
  }

  void cal_prev_pregant (const hivC<T>& hivpop, const artC<T>& artpop) 
  {
    update_active_pop_to(s->year); 
    update_active_last_year();
    epp::vector<T> n_mean(s->pAG_FERT); // 1 X 35
    for (int age = 0; age < s->pAG_FERT; ++age)
      n_mean(age) = (active_last_year_(age, s->F, s->N) + data_active(age, s->F, s->N)) / 2.;
    
    epp::index<1> sta = {int(s->p_fert_(0) - 1)}, len = {int(s->pAG_FERT - (s->p_fert_(0) - 1))};
    epp::vector<T> hivn = sumByAG1D<T>(n_mean, s->ag_.slice(sta, len), s->hAG_FERT);

    epp::index<2> O2 = {0,0}, X2 = {s->hDS, s->hAG_FERT};
    epp::matrix<T>
      mid_hiv = (hivpop.data_prev.chip(s->F, 2) + hivpop.data_curr.chip(s->F, 2)).slice(O2, X2) / 2.,
      frp = mid_hiv * p->nh.frr_cd4.chip(s->year, 2);
  
    epp::index<3> O3 = {0,0,0}, X3 = {s->hTS, s->hDS, s->hAG_FERT};
    epp::cube<T> 
      mid_art = (artpop.data_prev.chip(s->F, 3) + artpop.data_curr.chip(s->F, 3)).slice(O3, X3) / 2.,
      fra = mid_art * p->nh.frr_art.chip(s->year, 3);

    epp::vector<T> 
      frap = birth_agrp * (1. - hivn / (hivn + frp.sum(epp::by_rows) + fra.sum(epp::by_bay)));
    
    pregprev(s->year) = sum_vector(frap) / sum_vector(birth_age);
    if (s->year + s->AGE_START <= s->PROJ_YEARS - 1)
      pregprevlag(s->year + s->AGE_START - 1) = pregprev(s->year);
  }

  void save_prev_n_inc () 
  {
    if (s->year + s->AGE_START > s->PROJ_YEARS - 1) // repeat of cal_prev_pregant
      update_active_last_year();
    T n_positive = 0, everyone_now = 0, s_previous = 0;
    for (int sex = 0; sex < s->NG; sex++)
      for (int age = s->p_age15to49_(0) - 1; age < s->pAG_1549; age++) {
        for (int ds = 0; ds < s->pDS; ds++)
          everyone_now += data_curr(age, sex, ds);
        n_positive += data_curr(age, sex, s->P); // + included virgin positive
        s_previous += data_prev(age, sex, s->N); // susceptible including virgin
      }
    prev15to49(s->year) = n_positive / everyone_now;
    incid15to49(s->year) /= s_previous;
    // prev(s->year) = accu(data_all.slice(s.hivp_idx)) / accu(data_all);
    // incid(s->year) = incid15to49(s->year) / accu(data(s->year-1).slice(s.hivn_idx)); // toBfixed
  }

  T calc_rtrend_rt (int ts, T time_step) 
  {
    T rveclast = rvec(ts-1);
    Rf_error("K not write for r_trend, no fp template to write");
    // T dtii =  1 - s->DT * (time_step - 1);
    // int a_l = s->p_age15to49_(0) -1, a_r = a_l + s->p_age15to49_.num_elements();
    // boost2D A = data[ indices[s->year][s->N][in(0, s->NG)][in(a_l, a_r)] ];
    // boost1D B = data[ indices[s->year][s->N][in(0, s->NG)][a_l] ];
    // boost1D C = data[ indices[s->year][s->N][in(0, s->NG)][a_r] ];
    // T hivn_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
    // A = data[ indices[s->year][s->P][in(0, s->NG)][in(a_l, a_r)] ];
    // B = data[ indices[s->year][s->P][in(0, s->NG)][a_l] ];
    // C = data[ indices[s->year][s->P][in(0, s->NG)][a_r] ];
    // T hivp_ii = sumArray(A) - sumArray(B) * dtii + sumArray(C) * dtii;
    // T prevcurr = hivp_ii / (hivn_ii + hivp_ii);
    // T t_ii     = p->proj_steps[ts];
    // if (t_ii > p->tsEpidemicStart) {
      // par = p->rtrend;
      // if (t_ii < par.tStabilize)
      //   gamma_t = 0
      // else 
      //   gamma_t = (prevcurr-prevlast) * (t_ii - par$tStabilize) / (s->DT * prevlast)
      // logr.diff =  par$beta[2] * (par$beta[1] - rveclast) +
      //              par$beta[3] * prevlast + 
      //              par$beta[4] * gamma_t
      // return(exp(log(rveclast) + logr.diff))
    // }
    // else
      // return(p.rtrend$r0)
    return rveclast; // remove this
  }

  void update_rvec (T time_step) {
    int ts = (s->year-1) / s->DT + time_step;
    // if (p.eppmod %in% c("rtrend", "rtrend_rw")) <<- need to be fixed in FP
    if (p->ic.eppmod == 1) // rtrend see prepare_fp_for_Cpp
      rvec(ts) = calc_rtrend_rt(ts, time_step);
  }

  void update_infection () {
    epp::matrix<T> infection_ts = infections_ * infections_.constant(s->DT);
    data_curr.chip(s->N, 2) -= infection_ts;
    data_curr.chip(s->P, 2) += infection_ts;
    infections.chip(s->year, 2) += infection_ts;
    epp::index<2> origin = {int(s->p_age15to49_(0) - 1), 0}, 
                  xtents = {int(s->pAG_1549 - (s->p_age15to49_(0) - 1)), s->NG};
    epp::scalar<T> n1549 = infection_ts.slice(origin, xtents).sum(epp::by_bay);
    incid15to49(s->year) += n1549(0);
  }

  void remove_hiv_death (const hivC<T>& hivpop, const artC<T>& artpop) 
  {
    epp::matrix<T> 
      hivD_as = hivpop.death_.sum(epp::by_rows), 
      artD_as = artpop.death_.sum(epp::by_bay);
    if (s->MODEL==2) {
      hivD_as += hivpop.death_db_.sum(epp::by_rows);
      artD_as += artpop.death_db_.sum(epp::by_bay);
    }

    epp::matrix<T> 
      death_by_agrp_ = (hivD_as + artD_as) * s->DT,
      hiv_by_agrp_   = sumByAG2D<T>(data_curr.chip(s->P, 2), s->ag_, s->hAG),
      hiv_mx_agr = death_by_agrp_ / hiv_by_agrp_,
      hiv_mx_age = expand_age_group<T>(hiv_mx_agr, s->h_ag_span);
    replace_na_with(hiv_mx_age);

    hivdeaths.chip(s->year, 2) += data_curr.chip(s->P, 2) * hiv_mx_age;
    data_curr.chip(s->P, 2)    *= (1 - hiv_mx_age);

    if (s->MODEL==2)
      for (int sex = 0; sex < s->NG; sex++)
        for (int age = 0; age < s->pDB; age++) 
          data_db(age, sex, s->P, s->year) *= (1 - hiv_mx_age(age, sex));
  }

  void update_preg (const hivC<T>& hivpop, const artC<T>& artpop) 
  {
    update_active_pop_to(s->year);

    int p_lo = s->p_fert_(0) - 1, p_len = s->pAG_FERT - p_lo,
        h_lo = s->h_fert_(0) - 1, h_len = s->hAG_FERT - h_lo; 

    epp::index<1> start = {p_lo}, length = {p_len};
    epp::vector<T>
      sub_id     = s->ag_.slice(start, length),
      female_neg = data_active.chip(s->N, 2).chip(s->F, 1).slice(start, length),
      hivn       = sumByAG1D<T>(female_neg, sub_id, s->hAG_FERT);

    epp::cube<T> _art = artpop.data_curr.chip(s->F, 3).slice(epp::index<3>({0,0,0}), epp::index<3>({s->hTS, s->hDS, s->hAG_FERT})) * p->nh.frr_art.chip(s->year, 3);
    epp::matrix<T> _hiv = hivpop.data_curr.chip(s->F, 2).slice(epp::index<2>({0,0}), epp::index<2>({s->hDS,s->hAG_FERT})) * p->nh.frr_cd4.chip(s->year, 2);
    epp::vector<T> all_art = _art.sum(epp::by_bay) + _hiv.sum(epp::by_rows);

    for (int agr = h_lo; agr < s->hAG_FERT; agr++)
      for (int cd4 = 0; cd4 < (p->ad.artcd4elig_idx[s->year] - 1); cd4++)
        art_elig_(cd4, agr, s->F) +=
          hivpop.data_curr(cd4, agr, s->F) * p->nh.frr_cd4(cd4, agr, s->year) * 
            (birth_agrp(agr) / (hivn(agr) + all_art(agr)));
  }

  void art_initiate (const epp::vector<T>& art_curr, int time_step) 
  {
    int year_w = (s->DT * (time_step + 1) < 0.5) ? 0 : 1, 
        year_l = s->year - (2 - year_w), year_r = s->year - (1 - year_w);
    T trans = s->DT * (time_step + 1) + 0.5 - year_w;
    for (int sex = 0; sex < s->NG; ++sex) 
    {
      if( (p->ad.art15plus_isperc(sex, year_l) == 0) & 
          (p->ad.art15plus_isperc(sex, year_r) == 0) ) 
      { // both number
        n_art_ii_(sex) = p->ad.art15plus_num(sex, year_l) * (1 - trans) +
                         p->ad.art15plus_num(sex, year_r) * trans;
        if (s->MIX) {
          epp::matrix<T> art_elig_sex = art_elig_.chip(sex, 2);
          T cov = n_art_ii_(sex) / (sumArray(art_elig_sex) + art_curr(sex));
          artcov(sex) = (cov <= 1) ? cov : 1;
        }
      } 
      else if ( (p->ad.art15plus_isperc(sex, year_l) == 1) &
                (p->ad.art15plus_isperc(sex, year_r) == 1) ) 
      { // both percentage
        T cov = p->ad.art15plus_num(sex, year_l) * (1 - trans) + 
                     p->ad.art15plus_num(sex, year_r) * trans;
        if (s->MIX)
          artcov(sex) = (cov <= 1) ? cov : 1;
        epp::matrix<T> art_elig_sex = art_elig_.chip(sex, 2);
        n_art_ii_(sex) = cov * (sumArray(art_elig_sex) + art_curr(sex));
      }
      else if ( (p->ad.art15plus_isperc(sex, year_l) == 0) & 
                (p->ad.art15plus_isperc(sex, year_r) == 1) ) 
      { // transition number to percentage
        epp::matrix<T> art_elig_sex = art_elig_.chip(sex, 2);
        T 
          actual_cov = art_curr(sex) / (sumArray(art_elig_sex) + art_curr(sex)),
          diff_cov   = p->ad.art15plus_num(sex, year_r) - actual_cov,
          cov        = actual_cov + diff_cov * s->DT / (0.5 + year_w - s->DT * time_step);
        if (s->MIX)
          artcov(sex) = (cov < 1) ? cov : 1;
        n_art_ii_(sex) = cov * ( sumArray(art_elig_sex) + art_curr(sex) );
      }
    }
  }


  void update_eligible_for_art() { 
    // this one does not depend on model state
    for (int i = 0; i < s->hDS; ++i) {
      T A = (i >= p->ad.artcd4elig_idx[s->year] - 1) ? 1 : 0;
      T B = (i >= 2) ? p->ad.who34percelig : 0;
      elig_art_(i) = 1 - (1 - A) * (1 - B) * (1 - p->ad.specpop_percelig[s->year]);
    }
  }

  // calculate, distribute eligible for ART, update grad, gradART
  // -----------------------------------------------------------------------------
  void epp_art_init (hivC<T>& hivpop, artC<T>& artpop, int time_step) 
  {
    artpop.grad_progress();
    artpop.art_dropout(hivpop); // pass hivpop to receive the drop out
    update_eligible_for_art();
    epp::cube<T> elig_art_xp = elig_art_.broadcast(epp::index<1>({s->hAG * s->NG})).reshape(hivpop.dims);
    art_elig_ = hivpop.data_curr * elig_art_xp;
    if ( (p->ad.pw_artelig[s->year] == 1) && (p->ad.artcd4elig_idx[s->year] > 1) )
      update_preg(hivpop, artpop); // add pregnant?
    if (s->MODEL == 2) // add sexual inactive but eligible for treatment
      art_elig_ += hivpop.data_db.chip(s->year, 3) * elig_art_xp;
    // calculate number to initiate ART and distribute
    artpop.update_current_on_art();
    art_initiate(artpop.art_by_sex_, time_step); // update n_art_ii_
    auto n_afford = n_art_ii_ - artpop.art_by_sex_;
    n_art_init_ = (n_afford > 0.0).select(n_afford, n_afford.constant(0));
    art_distribute(n_art_init_); // update art_init_ as well
    replace_na_with(art_init_);
    if (s->MODEL == 1) {
      epp::cube<T> hiv_grad_dt = hivpop.data_curr + hivpop.grad * s->DT;
      art_init_ = (art_init_ < hiv_grad_dt).select(art_init_, hiv_grad_dt);
    }
    if (s->MODEL == 2) // split the number proportionally for active and idle pop
      hivpop.distribute_artinit(art_init_, artpop);
    hivpop.grad -= art_init_ / s->DT;
    artpop.grad_init(art_init_);
  }

  // calculate ART initiation distribution
  void art_distribute (const Eigen::Tensor<T, 1>& art_need) 
  {
    if (!p->ad.med_cd4init_input[s->year]) {
      if (p->ad.art_alloc_method == 4L) { // by lowest CD4
        // Calculate proportion to be initiated in each CD4 category
        epp::vector<T> init_pr(s->NG);
        for (int cd4 = s->hDS - 1; cd4 > 0; --cd4) { //6->0
          epp::vector<T> elig_hm(s->NG); elig_hm.setZero();
          for (int sex = 0; sex < s->NG; sex++)
            for (int age = 0; age < s->hAG; age++)
              elig_hm(sex) += art_elig_(cd4, age, sex);
          if ( (elig_hm(s->M) == 0) && (elig_hm(s->F) == 0) )
            init_pr = elig_hm;
          else {
            for (int sex = 0; sex < s->NG; ++sex) {
              T x = art_need(sex) / elig_hm(sex);
              init_pr(sex) = ( (x > 1) | std::isnan(x) | std::isinf(x)) ? 1 : x;
            }
          }
          for (int sex = 0; sex < s->NG; ++sex)
            for (int agr = 0; agr < s->hAG; ++agr)
              art_init_(cd4, agr, sex) = art_elig_(cd4, agr, sex) * init_pr(sex);
        }
      } 
      else { // Spectrum Manual p168--p169, 
        int A = s->h_age15plus_(0) - 1;
        dvec artX(s->NG), artY(s->NG);
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = A; agr < s->hAG_15plus; agr++)
            for (int cd4 = 0; cd4 < s->hDS; cd4++) {
              artX[sex] += art_elig_(cd4, agr, sex) * p->nh.cd4_mort(cd4, agr, sex);
              artY[sex] += art_elig_(cd4, agr, sex);
            }
        T xx;
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = A; agr < s->hAG_15plus; agr++)
            for (int cd4 = 0; cd4 < s->hDS; cd4++) {
              xx = (p->nh.cd4_mort(cd4, agr, sex) / artX[sex] *
                    p->ad.art_alloc_mxweight +
                    ((1 - p->ad.art_alloc_mxweight) / artY[sex]) ) *
                    art_elig_(cd4, agr, sex) * art_need(sex);
              art_init_(cd4, agr, sex) =
                (xx > art_elig_(cd4, agr, sex)) ? art_elig_(cd4, agr, sex) : xx;
            }
      }
    }
    else {
      int CD4_LO[] = {500,  350, 250, 200, 100, 50,  0 };
      int CD4_UP[] = {1000, 500, 350, 250, 200, 100, 50};

      int j = p->ad.med_cd4init_cat[s->year] - 1; // R to C++

      T pr_below = (p->ad.median_cd4init[s->year] - CD4_LO[j]) / (CD4_UP[j] - CD4_LO[j]);

      dvec elig_below(s->NG);
      for (int sex = 0; sex < s->NG; sex++)
        for (int agr = 0; agr < s->hAG; agr++)
          elig_below[sex] += art_elig_(j, agr, sex) * pr_below;
      if (j < (s->hDS - 1))
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = 0; agr < s->hAG; agr++)
            for (int cd4 = j+1; cd4 < s->hDS; cd4++)
              elig_below[sex] += art_elig_(cd4, agr, sex);

      dvec elig_above(s->NG);
      for (int sex = 0; sex < s->NG; sex++) 
        for (int agr = 0; agr < s->hAG; agr++)
          elig_above[sex] += art_elig_(j, agr, sex) * (1.0 - pr_below);
      if (j > 0)
        for (int sex = 0; sex < s->NG; sex++) 
          for (int agr = 0; agr < s->hAG; agr++)
            for (int cd4 = 0; cd4 <= (j - 1); cd4++)
              elig_above[sex] += art_elig_(cd4, agr, sex);

      dvec initpr_below(s->NG), initpr_above(s->NG), initpr_medcat(s->NG);
      T x, y;
      for (int sex = 0; sex < s->NG; ++sex) {
        x = art_need(sex) * 0.5 / elig_below[sex];
        y = art_need(sex) * 0.5 / elig_above[sex];
        initpr_below[sex] = (x > 1 || std::isnan(x) || std::isinf(x)) ? 1 : x;
        initpr_above[sex] = (y > 1 || std::isnan(y) || std::isinf(y)) ? 1 : y;
        initpr_medcat[sex] = initpr_below[sex] *      pr_below + 
                             initpr_above[sex] * (1 - pr_below);
      }
      if (j < (s->hDS - 1)) 
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = 0; agr < s->hAG; agr++)
            for (int cd4 = j + 1; cd4 < s->hDS; cd4++)
              art_init_(cd4, agr, sex) = art_elig_(cd4, agr, sex) * initpr_below[sex];
      
      for (int sex = 0; sex < s->NG; sex++)
        for (int agr = 0; agr < s->hAG; agr++)
          art_init_(j, agr, sex) = art_elig_(j, agr, sex) * initpr_medcat[sex];
      
      if (j == 0)
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = 0; agr < s->hAG; agr++)
              art_init_(0, agr, sex) = art_elig_(0, agr, sex) * initpr_above[sex];
      else if (j > 0)
        for (int sex = 0; sex < s->NG; sex++)
          for (int agr = 0; agr < s->hAG; agr++)
            for (int cd4 = 0; cd4 <= (j - 1); cd4++)
              art_init_(cd4, agr, sex) = art_elig_(cd4, agr, sex) * initpr_above[sex];
    }
  }

  epp::matrix<T> age_sex_cov (hivC<T>& hivpop, artC<T>& artpop) 
  {
    epp::matrix<T> 
      onn = artpop.n_by_agr(),
      off = hivpop.n_by_agr(),
      out(s->pAG, s->NG);
    for (int sex = 0; sex < s->NG; ++sex)
      for (int age = 0; age < s->pAG; ++age) {
        T cov = onn(s->ag_(age)-1, sex) / (onn(s->ag_(age)-1, sex) + off(s->ag_(age)-1, sex));
        out(age, sex) = std::isnan(cov) ? 0 : cov;
      }
    return out;
  }

  void infect_spec (const hivC<T>& hivpop, const artC<T>& artpop, int time_step) 
  { // TODO: rewrite in Tensor
    int ts = (s->year-1)/ s->DT + time_step,
        p_lo = s->p_age15to49_(0) - 1, h_lo = s->h_age15to49_(0) - 1;
    double dt_ii = 1 - s->DT * time_step, // transition of population in 1 year
           n_neg_mf = 0, n_pos_mf = 0, n_pos_inactive = 0, n_pos_inactive_lo = 0,
           n_pos_lo = 0, n_pos_up = 0, n_neg_lo = 0, n_neg_up = 0,
           n_hiv_lo = 0, n_art_lo = 0, n_hiv_up = 0, n_art_up = 0, art_ii = 0;

    update_active_pop_to(s->year); // substract virgin when needed

    for (int sex = 0; sex < s->NG; sex++) {
      for (int age = p_lo; age < s->pAG_1549; age++) {
        n_neg_mf += data_curr(age, sex, s->N); //  this includes debut neg
        if (s->MODEL == 2 && age < s->pDB)
          n_pos_inactive += data_db(age, sex, s->P, s->year); // "safe" positive
        n_pos_mf += data_active(age, sex, s->P); // transmissible positive
      }
      n_neg_lo += data_curr(p_lo, sex, s->N);
        if (s->MODEL == 2)
          n_pos_inactive_lo += data_db(p_lo, sex, s->P, s->year); // "safe" positive
      n_neg_up += data_curr(s->pAG_1549, sex, s->N);
      n_pos_lo += data_active(p_lo, sex, s->P);
      n_pos_up += data_active(s->pAG_1549, sex, s->P);
      if (s->year >= s->tARTstart-1) {
        for (int agr = h_lo; agr < s->hAG_1549; agr++)
          for (int cd4 = 0; cd4 < s->hDS; cd4++)
            for (int dur = 0; dur < s->hTS; dur++) {
              art_ii   += artpop.data_curr(dur, cd4, agr, sex);
              n_art_lo += artpop.data_curr(dur, cd4, h_lo, sex);
              n_art_up += artpop.data_curr(dur, cd4, s->hAG_1549, sex);
            }
      }
      for (int cd4 = 0; cd4 < s->hDS; cd4++) {
        n_hiv_lo += hivpop.data_curr(cd4, h_lo, sex);
        n_hiv_up += hivpop.data_curr(cd4, s->hAG_1549, sex);
      }
    }
    
    double hivp_inactive = n_pos_inactive - n_pos_inactive_lo * dt_ii;
    double hivn_both = n_neg_mf - n_neg_lo * dt_ii + n_neg_up * dt_ii;
    double hivp_active = n_pos_mf - n_pos_lo * dt_ii + n_pos_up * dt_ii;
    
    if (s->year >= s->tARTstart-1) {
      if (n_hiv_lo + n_art_lo > 0) {
        for (int sex = 0; sex < s->NG; sex++) {
          double art_trans = 0, hiv_trans = 0;
          for (int cd4 = 0; cd4 < s->hDS; cd4++) {
            for (int dur = 0; dur < s->hTS; dur++)
              art_trans += artpop.data_curr(dur, cd4, h_lo, sex);
            hiv_trans += hivpop.data_curr(cd4, h_lo, sex);
          }
          art_ii -= ( data_active(p_lo, sex, s->P) * 
                      art_trans / (hiv_trans + art_trans) ) * dt_ii;
        }
      }
      if (n_hiv_up + n_art_up > 0) {
        for (int sex = 0; sex < s->NG; sex++) {
          double art_trans = 0, hiv_trans = 0;
          for (int cd4 = 0; cd4 < s->hDS; cd4++) {
            for (int dur = 0; dur < s->hTS; dur++)
              art_trans += artpop.data_curr(dur, cd4, s->hAG_1549, sex);
            hiv_trans += hivpop.data_curr(cd4, s->hAG_1549, sex);
          }
          art_ii += ( data_active(s->pAG_1549, sex, s->P) * 
                      art_trans / (hiv_trans + art_trans) ) * dt_ii;
        }
      }
    }
    
    // Prob of contacting a sexual active, H+ in total pop
    double
    transm_prev = (hivp_active - art_ii * (1 - p->ic.relinfectART)) / 
                  (hivn_both + hivp_active + hivp_inactive);

    double w = (p->ic.proj_steps[ts] == p->ic.tsEpidemicStart) ? p->ic.iota : 0.0;
    double inc_rate = rvec(ts) * transm_prev + w;

    // Incidence: male = negative / female negative * sexratio + male negative; 
    //          female = male * sexratio
    double n_neg_m = 0, n_neg_f = 0;
    for (int age = p_lo; age < s->pAG_1549; age++) {
      n_neg_m += data_active(age, s->M, s->N);
      n_neg_f += data_active(age, s->F, s->N);
    }
    double adj_sex = (n_neg_m + n_neg_f) / (n_neg_m + n_neg_f * p->ic.incrr_sex[s->year]);
    double sex_inc[2] = {inc_rate * adj_sex, inc_rate * adj_sex * p->ic.incrr_sex[s->year]};
    // New infections distributed by age: ratio age_i/ 25-29 age
    for (int sex = 0; sex < s->NG; sex++) {
      double n_neg = 0, n_neg_rr = 0, adj_age;
      for (int age = p_lo; age < s->pAG_1549; age++) {
        n_neg += data_active(age, sex, s->N); 
        n_neg_rr += p->ic.incrr_age(age, sex, s->year) * data_active(age, sex, s->N);
      }
      adj_age = sex_inc[sex] / ( n_neg_rr / n_neg );
      for (int age = 0; age < s->pAG; age++) {
        if (sex == s->M) // age-specific incidence among circumcised men
          adj_age *= (1 - p->ic.circ_incid_rr * p->ic.circ_prop(age, s->year));
        infections_(age, sex) = p->ic.incrr_age(age, sex, s->year) * adj_age * 
          data_active(age, sex, s->N);
      }
    }

    // saving
    incrate15to49_ts(ts) = inc_rate;
    prev_last = (hivp_active + hivp_inactive) / (hivn_both + hivp_active + hivp_inactive);
    prev15to49_ts(ts) = prev_last;
  }

  void infect_mix (hivC<T>& hivpop, artC<T>& artpop, int ii) 
  {
    update_active_pop_to(s->year);
    
    epp::matrix<T> 
      N1 = hivpop.data_curr.sum(epp::by_rows) + artpop.data_curr.sum(epp::by_bay),
      N2 = hivpop.stage0.chip(s->year, 2) + (N1 - hivpop.stage0.chip(s->year, 2)) * p->ic.rel_vl[1],
      PP = data_active.chip(s->P, 2) / expand_age_group<T>(N1, s->h_ag_span),
      hiv_cd4_adj = PP * expand_age_group<T>(N2, s->h_ag_span);
    replace_na_with(hiv_cd4_adj);
    
    epp::index<1> ages = {s->pAG}, agexsex = {s->pAG * s->NG};
    data_active *= p->ic.est_senesence.reshape(agexsex).broadcast(ages);

    // balancing number of sex acts
    epp::matrix<T> 
      /* number of partnerships */
      nc_m = p->ic.mixmat.chip(s->M, 2) * p->ic.est_pcr.chip(s->M, 1).broadcast(ages),
      nc_f = p->ic.mixmat.chip(s->F, 2) * p->ic.est_pcr.chip(s->F, 1).broadcast(ages),
      /* number of pns formed by total pop */
      p_activ    = data_active.sum(epp::by_rack),
      nc_m_total = nc_m * p_activ.chip(s->M, 1).broadcast(ages),
      nc_f_total = nc_f * p_activ.chip(s->F, 1).broadcast(ages),
      /* balancing ratio */
      ratio_mf = nc_m_total / nc_f_total.shuffle(epp::transpose),
      nc_m_adj = nc_m * ratio_mf.pow(p->ic.balancing - 1),
      nc_f_adj = nc_f * ratio_mf.shuffle(epp::transpose).pow(p->ic.balancing),
      /* on negative only  */
      neg_prop = data_active.chip(s->N, 2) / p_activ;

    nc_m_adj *= (neg_prop.chip(s->M, 1).broadcast(ages));
    nc_f_adj *= (neg_prop.chip(s->F, 1).broadcast(ages));
    
    epp::matrix<T> art_cov(s->pAG, s->NG);
    art_cov.setZero();
    if (s->year >= s->tARTstart-1)
      art_cov = age_sex_cov(hivpop, artpop);

    epp::matrix<T>
      transm_prev = hiv_cd4_adj * (1 - art_cov * p->ic.relinfectART) / p_activ;
    // only two modes: "eppspectrum"=0 or "transm"=1
    T RR = (p->ic.incidmod == 0) ? p->ic.incrr_sex[s->year] : p->ic.mf_transm_rr[s->year];
    int ts = (s->year-1)/s->DT + ii;
    if (p->ic.proj_steps[ts] == p->ic.tsEpidemicStart)
      transm_prev = p->ic.leading_ev * p->ic.iota;
    transm_prev *= transm_prev.constant(rvec(ts));
    transm_prev.chip(0,s->M+1) = transm_prev.chip(0,s->M+1) * RR;
    
    epp::index<2> c2r = {1, s->pAG};
    auto m_neg = data_active.chip(s->N, 2).chip(s->M, 1).broadcast(ages);
    auto f_neg = data_active.chip(s->N, 2).chip(s->F, 1).broadcast(ages);
    auto m_pre = transm_prev.chip(s->M, 1).reshape(c2r).broadcast(ages);
    auto f_pre = transm_prev.chip(s->F, 1).reshape(c2r).broadcast(ages);
    auto m_cod = p->ic.est_condom.chip(s->year, 2).chip(s->M, 1).broadcast(ages);
    auto f_cod = p->ic.est_condom.chip(s->year, 2).chip(s->M, 1).reshape(c2r).broadcast(ages);
    epp::matrix<T>
      inc_m = nc_m_adj * m_neg * f_pre * (1. - m_cod),
      inc_f = nc_f_adj * f_neg * m_pre * (1. - f_cod);

    infections_.chip(s->M, 1) = inc_m.sum(epp::by_cols);
    infections_.chip(s->F, 1) = inc_f.sum(epp::by_cols);

    // prev15to49_ts_m should use this one! now just store as below
    epp::matrix<T> n_pos = data_curr.chip(s->P, 2);
    prev15to49_ts(ts) = sumArray(n_pos) / sumArray(data_curr);
    prev_last = prev15to49_ts(ts);
  }

  void epp_disease_model_direct(hivC<T>& hivpop, artC<T>& artpop) 
  {
    int a_l, a_r;
    if (p->ic.incidpopage) { // incidence for 15+ population
      a_l = s->p_age15plus_(0) - 1;
      a_r = s->pAG_15plus;
    } else { // incidence for 15 -49 population
      a_l = s->p_age15to49_(0) - 1;
      a_r = s->pAG_1549;
    }
    update_active_last_year();
    double n_m = 0, n_f = 0;
    for (int age = a_l; age < a_r; ++age) {
      n_m += active_last_year_(age, s->M, s->N);
      n_f += active_last_year_(age, s->F, s->N);
    }
    dvec sex_inc(s->NG);
    sex_inc[s->M] = (n_m + n_f) * p->ic.incidinput[s->year] / (n_m + n_f  * p->ic.incrr_sex[s->year]);
    sex_inc[s->F] = (n_m + n_f) * p->ic.incidinput[s->year] * p->ic.incrr_sex[s->year] / (n_m + n_f  * p->ic.incrr_sex[s->year]);
    dvec ageinc(s->NG);
    for (int sex = 0; sex < s->NG; sex++) {
      double neg_sa = 0, inc_sa = 0;
      for (int age = a_l; age < a_r; age++) {
        neg_sa += active_last_year_(age, sex, s->N);
        inc_sa += active_last_year_(age, sex, s->N) * p->ic.incrr_age(age, sex, s->year);
      }
      ageinc[sex] = inc_sa / neg_sa;
    }
    double new_infect = 0;
    for (int sex = 0; sex < s->NG; sex++)
      for (int age = 0; age < s->pAG; age++) {
        new_infect = p->ic.incrr_age(age, sex, s->year) * ( sex_inc[sex] / ageinc[sex]) * active_last_year_(age, sex, s->N);
        infections(age, sex, s->year) = new_infect;
        data_curr(age, sex, s->N) -= new_infect;
        data_curr(age, sex, s->P) += new_infect;
      }
    epp::matrix<double> infect_agrp = sumByAG2D<double>(infections.chip(s->year, 3), s->ag_, s->hAG);
    for (int sex = 0; sex < s->NG; sex++)
      for (int agr = 0; agr < s->hAG; agr++)
        for (int cd4 = 0; cd4 < s->hDS; cd4++)
          hivpop.data_curr(cd4, agr, sex) += p->nh.cd4_initdist(cd4, agr, sex) * infect_agrp(agr, sex);
    for (int sex = 0; sex < s->NG; sex++)
      for (int age = s->p_age15to49_(0) - 1; age < s->pAG_1549; age++)
        incid15to49(s->year) += infections(age, sex, s->year);
  }

};
