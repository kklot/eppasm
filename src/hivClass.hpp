#include "fpClass.hpp"
#pragma once

template <typename T> class artC;

// HIV class
template<typename T>
class hivC {
public: // fields
  StateSpace* const s;
  Parameters* const p;
  int N;
  Eigen::array<long, 3> dims;
  epp::array4D<T> data;
  epp::array4D<T> data_db; // debut only population
  epp::map_of_cube<T> data_curr;
  epp::map_of_cube<T> data_prev;
  epp::cube<T> grad;
  epp::cube<T> grad_db;
  epp::cube<T> death_; // death in this year
  epp::cube<T> death_db_; // death in this year
  epp::cube<T> cd4_mort_;
  epp::matrix<T> infect_by_agrp_;
  epp::cube<T> stage0;
public: // inits
  hivC(StateSpace* const s, Parameters* const p) :
    s              (s),
    p              (p),
    N              (s->hDS * s->hAG * s->NG),
    dims           ({s->hDS, s->hAG, s->NG}),
    data           (s->hDS, s->hAG, s->NG, s->PROJ_YEARS),
    data_db        (s->hDS, s->hAG, s->NG, s->PROJ_YEARS),
    data_curr      (&data(0), s->hDS, s->hAG, s->NG),
    data_prev      (&data(0), s->hDS, s->hAG, s->NG),
    grad           (dims),
    grad_db        (dims),
    death_         (dims),
    death_db_      (dims),
    cd4_mort_      (dims),
    infect_by_agrp_(s->hAG, s->NG),
    stage0         (s->hAG, s->NG, s->PROJ_YEARS)
  {
    data.setZero();
    data_db.setZero();
    stage0.setZero();
    cd4_mort_ = p->nh.cd4_mort; // copy bc. this can be scaled
    if (s->MODEL==2) {
      epp::cube<T> stage0_0_xp = expand_2to3<T>(p->ic.stage0_0, s->hDS);
      data_curr = p->ic.stages_0 + stage0_0_xp * p->nh.cd4_initdist;
      stage0.chip(0, 2) = p->ic.stage0_0;
    }
  };
// methods
  void update_views(int t) { // #TutorialMapPlacementNew
    new (&data_curr) epp::map_of_cube<T>(&data(0) + N * t, dims);
    new (&data_prev) epp::map_of_cube<T>(&data(0) + N * (t - 1), dims);
  }

  void aging(const epp::matrix<T>& ag_prob) 
  {
    data_curr = data_prev;
    epp::cube<T> nup = data_prev * expand_2to3<T>(ag_prob, s->hDS);
    epp::index<3> origin_lo = {0,0,0}, origin_up = {0,1,0},  xtents = {s->hDS, s->hAG-1, s->NG};
    data_curr.slice(origin_lo, xtents)   -= nup.slice(origin_lo, xtents);
    data_curr.slice(origin_up, xtents)   += nup.slice(origin_lo, xtents);

    if (s->MODEL == 2) 
    {
      data_db.chip(s->year, 3) = data_db.chip(s->year-1, 3);
      nup = data_db.chip(s->year-1, 3) * expand_2to3<T>(ag_prob, s->hDS);
      data_db.chip(s->year, 3).slice(origin_lo, xtents) -= nup.slice(origin_lo, xtents);
      data_db.chip(s->year, 3).slice(origin_up, xtents) += nup.slice(origin_lo, xtents);

      // moving the newly infected remained after the 10th time step last year
      stage0.chip(s->year, 2) = stage0.chip(s->year-1, 2);
      epp::matrix<T> s0up = stage0.chip(s->year-1, 2) * ag_prob;
      epp::index<2> originx = {0,0}, originy = {1,0}, xtentsx = {s->hAG-1, s->NG};
      stage0.chip(s->year, 2).slice(originx, xtentsx) -= s0up.slice(originx, xtentsx);
      stage0.chip(s->year, 2).slice(originy, xtentsx) += s0up.slice(originx, xtentsx);
    }
  }

  void add_entrants(const epp::vector<T>& artYesNo) 
  {
    epp::matrix<T> add(s->hDS, s->NG);
		for(int sex = 0; sex < s->NG; sex++)
			for(int cd4 = 0; cd4 < s->hDS; cd4++)
				add(cd4, sex) = p->ph.paedsurv_cd4dist(cd4, sex, s->year) * artYesNo(sex+2);
    if (s->MODEL == 1)
      data_curr.chip(0, 1) += add;
    if (s->MODEL == 2) // add to virgin then debut
      data_db.chip(s->year, 3).chip(0, 1) += add;
  }

  void sexual_debut() 
  {
    epp::index<3> origin = {0,0,0}, xtents = {s->hDS,s->hDB,s->NG};
    epp::cube<T> 
			dbxp = expand_2to3<T>(p->ic.db_rate.chip(s->year, 2), s->hDS), 
			ndb = data_db.chip(s->year, 3).slice(origin, xtents) * dbxp;
    data_curr.slice(origin, xtents)               += ndb;
    data_db.chip(s->year, 3).slice(origin, xtents) -= ndb;
  }

  void deaths (const epp::matrix<T>& survival_pr) 
  {
    epp::cube<T> sv = expand_2to3(survival_pr, s->hDS);
    data_curr *= sv;
    if (s->MODEL == 2)
      data_db.chip(s->year, 3) *= sv;
    stage0.chip(s->year,2) = stage0.chip(s->year,2) * survival_pr;
  }

  void migration (const epp::matrix<T>& migration_pr)
  {
    epp::cube<T> mg = expand_2to3(migration_pr, s->hDS);
    data_curr *= mg;
    if (s->MODEL == 2)
      data_db.chip(s->year, 3) *= mg;
    stage0.chip(s->year,2) = stage0.chip(s->year,2) * migration_pr;
  }

  void update_infection (const epp::matrix<T>& new_infect) 
  {
    infect_by_agrp_ = sumByAG2D<T>(new_infect, s->ag_, s->hAG);
    if (s->MODEL == 2) // tracking the newly infected at this time step
      stage0.chip(s->year, 2) += infect_by_agrp_ * infect_by_agrp_.constant(s->DT);
    grad = p->nh.cd4_initdist * expand_2to3<T>(infect_by_agrp_, s->hDS);
  }

  void scale_cd4_mort(artC<T>& artpop) 
  {
    epp::index<1> over_row({0});
    epp::cube<T> num = data_curr, 
                 den = artpop.data_curr.sum(over_row);
    if (s->MODEL==2) {
      num += data_db.chip(s->year, 3);
      den += artpop.data_db.chip(s->year, 4).sum(over_row);
    }
    num /= (num + den);
    replace_na_with(num, 1.);
    cd4_mort_ = num * p->nh.cd4_mort;
  }

  void grad_progress() 
  {
    epp::index<3> 
      origin_lo = {0,0,0}, origin_up = {1,0,0}, xtents = {s->hDS-1, s->hAG, s->NG};
    death_  = data_curr * cd4_mort_;
    grad   -= death_;
    epp::cube<T> nup = data_curr.slice(origin_lo, xtents) * p->nh.cd4_prog;
    grad.slice(origin_lo, xtents) -= nup;
    grad.slice(origin_up, xtents) += nup;
    if (s->MODEL == 2) {
      stage0.chip(s->year, 2) = stage0.chip(s->year, 2) * (1 - exp(-1/(p->ic.stage0_time * 12 * s->DT)));
      death_db_  = data_db.chip(s->year, 3) * cd4_mort_;
      grad_db    = -death_db_;
      nup = data_db.chip(s->year, 3).slice(origin_lo, xtents) * p->nh.cd4_prog;
      grad_db.slice(origin_lo, xtents) -= nup;
      grad_db.slice(origin_up, xtents) += nup;
    }
  }

  void distribute_artinit (epp::cube<T>& artinit, artC<T>& artpop) 
  {
    epp::cube<T> 
      debut_now = data_db.chip(s->year, 3) + grad_db * grad_db.constant(s->DT),
      all_hivpop = data_curr + grad * grad.constant(s->DT) + debut_now;
    epp::cube<bool> ifs = artinit > all_hivpop;
      artinit = ifs.select(all_hivpop, artinit);
    epp::cube<T> 
      pr_weight_db = debut_now / all_hivpop, // this can be NaN
      n_artinit_db = artinit * pr_weight_db;
      replace_na_with(n_artinit_db);
    artinit -= n_artinit_db;
    grad_db -= n_artinit_db / n_artinit_db.constant(s->DT);
    artpop.grad_db_init(n_artinit_db, s);
  }

  void add_grad_to_pop () {
    data_curr += s->DT * grad;
    if (s->MODEL == 2)
      data_db.chip(s->year, 3) += s->DT * grad_db;
  }

  void adjust_pop(const epp::matrix<T>& adj_prob) {
    epp::cube<T> adj = expand_2to3(adj_prob, s->hDS);
    data_curr *= adj;
    if (s->MODEL == 2)
      data_db.chip(s->year, 3) *= adj;
  }

  epp::matrix<T> n_by_agr() {
    epp::matrix<T> out = data_curr.sum(epp::index<1>({0}));
    if (s->MODEL == 2)
      out += data_db.chip(s->year, 3).sum(epp::index<1>({0}));
    return out;
  }
};
