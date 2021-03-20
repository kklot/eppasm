#include "fpClass.hpp"
#pragma once

// ART class
template<typename T>
class artC {
public: // fields
  StateSpace* const ss;
  Parameters* const ps;
  int                   N;
  Eigen::array<long, 4> dims;
  epp::array5D<T> data;
  epp::array5D<T> data_db;
  epp::map_of_4D<T> data_curr;
  epp::map_of_4D<T> data_prev;
  epp::array4D<T> gradART;
  epp::array4D<T> gradART_db;
  epp::array4D<T> death_;
  epp::array4D<T> death_db_;
  epp::vector<T> art_by_sex_;
public: // Inits
  artC(StateSpace* const s, Parameters* const p) :
    ss          (s),
    ps          (p),
    N           (ss->hTS * ss->hDS * ss->hAG * ss->NG),
    dims        ({ss->hTS, ss->hDS, ss->hAG, ss->NG}),
    data        (ss->hTS, ss->hDS, ss->hAG, ss->NG, ss->PROJ_YEARS),
    data_db     (ss->hTS, ss->hDS, ss->hAG, ss->NG, ss->PROJ_YEARS),
    data_curr   (&data(0), ss->hTS, ss->hDS, ss->hAG, ss->NG),
    data_prev   (&data(0), ss->hTS, ss->hDS, ss->hAG, ss->NG),
    gradART     (dims),
    gradART_db  (dims),
    death_      (dims),
    death_db_   (dims),
    art_by_sex_ (ss->NG)
  {
    data.setZero();
    data_db.setZero();
    death_.setZero();
    death_db_.setZero();
  }
// Methods
  void update_views(int t) { // #TutorialMapPlacementNew
    new (&data_curr) epp::map_of_4D<T>(&data(0) + N * t, dims);
    new (&data_prev) epp::map_of_4D<T>(&data(0) + N * (t - 1), dims);
  }

  void aging(const epp::matrix<T>& ag_prob) 
  {
    data_curr = data_prev;
    epp::array4D<T> 
      agxp = expand_2to4<T>(ag_prob, ss->hTS * ss->hDS).reshape(dims),
      nup = data_prev * agxp;
    epp::index<4> origin_lo = {0,0,0,0}, xtents = {ss->hTS, ss->hDS, ss->hAG - 1, ss->NG},
                  origin_up = {0,0,1,0};
    data_curr.slice(origin_lo, xtents) -= nup.slice(origin_lo, xtents);
    data_curr.slice(origin_up, xtents) += nup.slice(origin_lo, xtents);

    if (ss->MODEL == 2) { // using the same dimension in db!
      data_db.chip(ss->year, 4) = data_db.chip(ss->year-1, 4);
      nup = data_db.chip(ss->year-1, 4) * agxp;
      data_db.chip(ss->year, 4).slice(origin_lo, xtents) -= nup.slice(origin_lo, xtents);
      data_db.chip(ss->year, 4).slice(origin_up, xtents) += nup.slice(origin_lo, xtents);
    }
  }

  void add_entrants(const Eigen::Tensor<T, 1>& artYesNo) 
  {
    epp::cube<T> add(ss->hTS, ss->hDS, ss->NG);
    for(int sex = 0; sex < ss->NG; sex++)
      for(int cd4 = 0; cd4 < ss->hDS; cd4++)
        for(int dur = 0; dur < ss->hTS; dur++)
          add(dur, cd4, sex) = ps->ph.paedsurv_artcd4dist(dur, cd4, sex, ss->year) * artYesNo(sex);
    if (ss->MODEL == 1)
      data_curr.chip(0, 2) += add;
    if (ss->MODEL == 2) // add to virgin then debut
      data_db.chip(ss->year, 4).chip(0, 2) += add;
  }

  void sexual_debut()
  {
    epp::index<4> origin = {0,0,0,0}, xtents = {ss->hTS, ss->hDS, ss->hDB, ss->NG};
    epp::array4D<T> 
      dbxp = expand_2to4<T>(ps->ic.db_rate.chip(ss->year, 2), ss->hTS * ss->hDS).reshape(xtents),
      ndb = data_db.chip(ss->year, 4).slice(origin, xtents) * dbxp;
    data_curr.slice(origin, xtents)                 += ndb;
    data_db.chip(ss->year, 4).slice(origin, xtents) -= ndb;
  }

  void deaths(const epp::matrix<T>& survival_pr) 
  {
    epp::array4D<T> sv = expand_2to4<T>(survival_pr, ss->hDS * ss->hTS).reshape(dims);
    data_curr *= sv;
    if (ss->MODEL == 2)
      data_db.chip(ss->year, 4) *= sv;
  }

  void migration(const epp::matrix<T>& migration_pr) 
  {
    epp::array4D<T> mg = expand_2to4<T>(migration_pr, ss->hDS * ss->hTS).reshape(dims);
    data_curr *= mg;
    if (ss->MODEL == 2)
      data_db.chip(ss->year, 4) *= mg;
  }

  void grad_progress() 
  {
    gradART = -death_;
    epp::index<4> 
      origin_lo = {0,0,0,0}, extent = {ss->hTS-1,ss->hDS,ss->hAG,ss->NG}, 
      origin_up = {1,0,0,0};
    epp::array4D<T> art_up = 2.0 * data_curr.slice(origin_lo, extent);
    gradART.slice(origin_lo, extent) -= art_up;
    gradART.slice(origin_up, extent) += art_up;
    
    if (ss->MODEL == 2) {
      gradART_db = -death_db_;
      art_up = 2.0 * data_db.chip(ss->year, 4).slice(origin_lo, extent);
      gradART_db.slice(origin_lo, extent) -= art_up;
      gradART_db.slice(origin_up, extent) += art_up;
    }
  }

  void art_dropout(hivC<T>& hivpop) 
  {
    epp::array4D<T> n_dropout = data_curr * ps->ad.art_dropout[ss->year];
    hivpop.grad += n_dropout.sum(epp::index<1>({0}));
    gradART     -= n_dropout;
    if (ss->MODEL == 2) {
      n_dropout       = data_db.chip(ss->year, 4) * ps->ad.art_dropout[ss->year];
      hivpop.grad_db += n_dropout.sum(epp::index<1>({0}));
      gradART_db     -= n_dropout;    
    }
  }

  void update_current_on_art() 
  {
    epp::array4D<T> art_ts = data_curr + gradART * ss->DT;
    if (ss->MODEL==2)
        art_ts += data_db.chip(ss->year, 4) + gradART_db * ss->DT;
    art_by_sex_ = art_ts.sum(epp::index<3>({0,1,2}));
  }

  void grad_init(const epp::cube<T>& artinit) {
    gradART.chip(0, 0) += artinit / ss->DT;
    data_curr += ss->DT * gradART;
  }

  void grad_db_init(const epp::cube<T>& artinit_db, StateSpace* const s) {
    gradART_db.chip(0, 0) += artinit_db / ss->DT;
    data_db.chip(ss->year, 4) += ss->DT * gradART_db;
  }

  void adjust_pop(const epp::matrix<T>& adj_prob) {
    epp::array4D<T> adj_xp = expand_2to4<T>(adj_prob, ss->hTS * ss->hDS).reshape(dims);
    data_curr *= adj_xp;
    if (ss->MODEL == 2)
      data_db.chip(ss->year, 4) *= adj_xp;
  }

  void count_death() {
    epp::array4D<T>
      mxxp = ps->nh.artmx_timerr.chip(ss->year, 1)
        .broadcast(epp::index<1>({ss->hDS * ss->hAG * ss->NG})).reshape(dims),
      artmx = ps->nh.art_mort * mxxp;
    death_ = data_curr * artmx;
    if (ss->MODEL==2)
      death_db_ = data_db.chip(ss->year, 4) * artmx;
  }

  epp::matrix<T> n_by_agr() {
    epp::matrix<T> out = data_curr.sum(epp::by_bay);
    if (ss->MODEL == 2)
      out += data_db.chip(ss->year, 4).sum(epp::by_bay);
    return out;
  }
};
